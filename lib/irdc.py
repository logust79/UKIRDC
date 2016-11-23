'''
irdc analysis
'''
import os
import sys
import logging
import sqlite3
sys.path.append('../BioTools')
from CommonFuncs import *

# store useful info into irdc.db
irdc_conn = sqlite3.connect('irdc.db')
irdc_c = irdc_conn.cursor()

# create tables if not exist
irdc_c.execute('''CREATE TABLE IF NOT EXISTS variants
    (id text, exac text, kaviar_af real, cadd_phred real)''')
#irdc_c.execute('''CREATE TABLE IF NOT EXISTS patients
#    (id text, hpos text, exac text, kaviar_af real, cadd_phred real)''')
#irdc_c.execute('''CREATE TABLE IF NOT EXISTS patient_variant
#    (id text, hpos text, exac text, kaviar_af real, cadd_phred real)''')
irdc_conn.commit()

def _update_db(table,field,value_dict):
    # update database
    remain_fields = list(set(['exac','kaviar_af','cadd_phred'])-set([field]))
    for k,v in value_dict.iteritems():
        irdc_c.execute('''INSERT OR REPLACE INTO %(table)s (id, %(field)s, %(f1)s, %(f2)s)
            VALUES (  ?,
            ?,
            (SELECT %(f1)s FROM variants WHERE id = ?),
            (SELECT %(f2)s FROM variants WHERE id = ?)
            )''' % {'table':table,'field':field,'f1':remain_fields[0],'f2':remain_fields[1]}, (k,v,k,k))
    irdc_conn.commit()


def _dict_factory(cursor, row):
    # convert from sqlite tuple to dictionary
    # if row is None, return None
    if row == None: return None
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

class Variant(object):
    def __init__(self, variant_id=None):
        '''
        attributes:
            variant_id
            cleaned_variant_id
        '''
        self.variant_id = variant_id if variant_id else logging.warning('Need variant_id')
        self.cleaned_variant_id = clean_variant(variant_id)
    
    @property
    def exac(self):
        # exac annotation for one variant.
        # check local database first. If not,
        #   use CommonFuncs to annotate exac, then store in database
        if getattr(self, '_exac', None) is None:
            irdc_c.execute('SELECT * FROM variants WHERE id=?',(self.cleaned_variant_id,))
            db_var = _dict_factory(irdc_c,irdc_c.fetchone())
            if db_var == None or db_var['exac'] == None:
                # query web
                print 'querying the web for exac'
                exac = anno_exac(self.cleaned_variant_id)
                # insert into database
                _update_db(
                           'variants',
                           'exac',
                           {self.cleaned_variant_id:json.dumps(exac,indent=4)}
                           )
            else:
                exac = json.loads(db_var['exac'])
            self._exac = exac
        return self._exac
    @property
    def kaviar_af(self):
        self._kaviar_af = 0.
        return self._kaviar_af

class Variants:
    def __init__(self, vars=[]):
        self.variants = vars
        self.cleaned_variants = {i:clean_variant(i) for i in vars}
    @property
    def exac(self):
        # Check local database first. If not,
        #   use CommonFuncs to annotate exac, then store in database
        if getattr(self, '_exac', None) is None:
            # check database
            sql="select * from variants where id in ({seq})".format(
                    seq=','.join(['?']*len(self.variants)))
            result = irdc_c.execute(sql,self.cleaned_variants.values())
            data = {}
            new_vars = {}
            exac = {}
            for i in result:
                temp = _dict_factory(irdc_c,i)
                data[temp['id']] = json.loads(temp['exac'])
            for k,v in self.cleaned_variants.iteritems():
                if v in data:
                    exac[k] = data[v]
                else:
                    # not in database, push to array for later query
                    new_vars[k] = v
            if new_vars:
                print 'querying web'
                new_result = anno_exac_bulk(new_vars.values())
                # update database
                _update_db(
                           'variants',
                           'exac',
                           {k:json.dumps(v,indent=4) for k,v in new_result.iteritems()}
                           )
                # populate exac
                for k,v in new_vars.iteritems():
                    exac[k] = new_result[v]
                irdc_conn.commit()
        self._exac = exac
        return self._exac
