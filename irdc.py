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
                if db_var == None:
                    irdc_c.execute('INSERT INTO variants VALUES (?,?,?,?)',(self.cleaned_variant_id,json.dumps(exac, indent=4),None,None))
                else:
                    irdc_c.excute('UPDATE variants SET exac=? where id=?',(json.dumps(exac, indent=4)), self.variant_id)
                irdc_conn.commit()
            else:
                exac = json.loads(db_var['exac'])
            self._exac = exac
        return self._exac

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
                for k,v in new_vars.iteritems():
                    irdc_c.execute('''INSERT OR REPLACE INTO variants (id, exac, kaviar_af, cadd_phred)
                        VALUES (  ?,
                        ?,
                        (SELECT kaviar_af FROM variants WHERE id = ?),
                        (SELECT cadd_phred FROM variants WHERE id = ?)
                        )''', (v,json.dumps(new_result[v],indent=4),v,v))
                    exac[k] = new_result[v]
                irdc_conn.commit()
        self._exac = exac
        return self._exac
