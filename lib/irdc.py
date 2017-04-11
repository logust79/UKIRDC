'''
irdc analysis
'''
import os
import sys
import logging
import sqlite3
import json
import csv
sys.path.append('../../BioTools')
import Variants
import Genes
import HPO
import Compare
import mysql.connector
import pymongo
import sqlite_utils
import re
import copy
import time
from collections import defaultdict
from itertools import combinations
from urllib2 import HTTPError, URLError
from Bio import Entrez
import pandas as pd
import numpy as np
import xlsxwriter as xw
from fields_update_methods import field_helpers

'''
globals
'''
EARLY_ONSET_HPOS = [
    'HP:0003623', #Neonatal onset
    'HP:0003577', #Congenital onset
    'HP:0011463', #Childhood onset
    'HP:0003621', #Juvenile onset
    'HP:0003593', #Infantile onset
    'HP:0030674', #Antenatal onset
    'HP:0011460', #Embryonal onset
    'HP:0011461', #Fetal onset
    ]

'''
given an array and a header, give back a dict
'''
def _make_dict(data,header):
    result = {}
    for i,d in enumerate(header):
        result[d] = data[i]
    return result

'''
check if two variants are close together to warrant a igv check
'''
def gap_check(l1,l2,gap_size):
    # check if l1 and l2 are close given gap_size
    if ((l1[0] <= l2[0] and l1[1] >= l2[0]) or # overlap
            (l1[0] >= l2[0] and l1[0] <= l2[1]) or # overlap
            (abs(l1[0] - l2[1]) <= gap_size) or  # gap
            (abs(l2[0] - l1[1]) <= gap_size)):    # gap
        return 'Y'
    else:
        return None

'''
check if i is unique among given files. Files have to be first filtered against proband and relatives
'''
def _get_exome_unique_cnv(i,files):
    bad = 0
    for f in files:
        inf = open(f,'r')
        header = []
        for row in inf:
            row = row.rstrip().split('\t')
            if not header:
                header = row
                continue
            row = _make_dict(row,header)
            # find chrom
            if row['chromosome'] != i['chrom']: continue
            # find overlap
            this_start = int(row['start'])
            this_end = int(row['end'])
            
            if max(this_end,i['end']) - min(this_start,i['start']) < sum([this_end - this_start,i['end'] - i['start']]) - 1:
                # overlap
                bad = 1
                break
        if bad == 1: break
    return bad

'''
check if given hpo (rod-cone or cone-rod) is in hpos's ancestors.
if yes, add to hpos and return
'''
def _check_add_hpo(hpo, hpos):
    hpo_ids = [i['id'] for i in hpos]
    if hpo['id'] not in hpo_ids:
        add_flag = 0
        for i in hpo_ids:
            H = HPO.Hpo(sqlite3.connect('irdc.db'),i)
            if hpo['id'] in H.ancestors:
                add_flag = 1
                break
        if add_flag:
            hpos.append(hpo)
    return hpos

'''
parse exac info
'''
def parse_exac(exac, hom=False):
    if 'allele_freq' in exac['variant']:
        # is covered and exists on exac
        if not hom:
            return exac['variant']['allele_freq']
        else:
            return float(exac['variant']['hom_count']) / (sum(exac['variant']['pop_ans'].values()) / 2 or 1)
    if exac['base_coverage'] and exac['base_coverage'][0]['has_coverage']:
        # is covered but not exsits on exac
        return 0
    # not covered
    return None

'''
parse gnomad info. merge exome and genome
'''
def parse_gnomad(gnomad, hom=False):
    result = None
    total_an = 0
    for mode in ['genomes','exomes']:
        if gnomad[mode]['any_covered']:
            result = result or 0
            if not gnomad[mode]['data']:
                continue
            total_an += gnomad[mode]['data'][0]['an']
            if hom:
                result += gnomad[mode]['data'][0]['hom'] * 2
            else:
                result += gnomad[mode]['data'][0]['ac']
    # not covered
    if not result:
        return result
    else:
        return float(result) / total_an

'''
parse protein atlas, given tissue
result = {
    ENSG0000:[{
        'cell_type':
        'level':
        'reliability':
    }]
}
'''
def parse_protein_atlas(file,tissue):
    result = {}
    csvreader = csv.reader(open(file,'r'), delimiter=',', quotechar='"')
    header = []
    for row in csvreader:
        if not header:
            header = row
            continue
        # tissue?
        if row[header.index('Tissue')] != tissue: continue

        g = row[header.index('Gene')]
        result[g] = result.get(g,[])
        result[g].append({
            'cell_type':row[header.index('Cell type')],
            'level':row[header.index('Level')],
            'reliability':row[header.index('Reliability')],
        })
    return result
'''
find start and end given a variant or cnv id
'''
def parse_pos(id,type):
    if type == 'snp':
        start = int(id.split('-')[1])
        end = start + len(id.split('-')[2]) - 1
    else:
        start = int(id.split('-')[0].split(':')[1])
        end = int(id.split('-')[1])
    return (start,end)

'''
format excel
'''
def format_excel(a,writer,changes,rules):
    """ Add Excel specific formatting to the workbook
    """
    # Get the workbook and the sheets so we can add the formatting
    workbook = writer.book
    # format wrap
    wrap_fmt = workbook.add_format({'text_wrap': True})
    rule_fmt = {k:workbook.add_format(v) for k,v in rules.items()}
    
    # format each sheet
    for sheet in a:
        header = list(a[sheet])
        worksheet = writer.sheets[sheet]
        shape = a[sheet].shape
        # highlight changes, wrap variants, pubmed
        if sheet in ['recessive','dominant','X']:
            
            # changes highlight, highlight first cell and changed cell
            if changes:
                for k1,v1 in changes[sheet].items():
                    if k1 == '<>':
                        for k2,v2 in v1.items():
                            row_ind = np.where(a[sheet]['gene_id'] == k2)[0][0]
                            cols = v2['change'].keys()
                            col_inds = [0] + [header.index(i) for i in cols]
                            for i in col_inds:
                                cell_value = a[sheet].get_value(row_ind,header[i])
                                worksheet.write(row_ind+1, i, cell_value, rule_fmt[k1])
                    else:
                        if not v1: continue
                        row_inds = a[sheet][a[sheet]['gene_id'].isin(v1)].index.tolist()
                        for i in row_inds:
                            cell_value = a[sheet].get_value(i,header[0])
                            worksheet.write(i+1, 0, cell_value, rule_fmt[k1])
        
            # columns highlight
            variant_index = header.index('variants')
            v_col = xw.utility.xl_col_to_name(variant_index)
            pubmed_index = header.index('pubmed')
            p_col = xw.utility.xl_col_to_name(pubmed_index)
            protein_index = header.index('protein_atlas')
            pr_col = xw.utility.xl_col_to_name(protein_index)
            worksheet.set_column(':'.join([v_col,v_col]),20, wrap_fmt)
            worksheet.set_column(':'.join([p_col,p_col]),10, wrap_fmt)
            worksheet.set_column(':'.join([pr_col,pr_col]),10, wrap_fmt)
            
            # add table. uncomment if necessary
            '''
            this_header = [{'header': di} for di in a[sheet].columns.tolist()]
            cell_range = xw.utility.xl_range(0,0,shape[0],shape[1])
            worksheet.add_table(cell_range,{'header_row': True,'columns':this_header})
            '''
        elif sheet == 'notes':
            # notes, wider and wrap
            worksheet.set_column('A:A',20,wrap_fmt)
            worksheet.set_column('B:B',40,wrap_fmt)
'''
write data to excel
'''
def write_to_excel(a,f,ped_loc,changes,rules):
    #f = '../analysis_data/test/result.xlsx'
    writer = pd.ExcelWriter(f,engine='xlsxwriter')
    for k in ['recessive','dominant','X','relatives','history']:
        if k == 'history':
            if k in a:
                print(a[k])
                a[k].to_excel(writer,sheet_name=k)
        else:
            a[k].to_excel(writer,sheet_name=k,index=False)
    # write notes. first transpose, then add pedigree figure if available
    a['notes'].transpose().to_excel(writer,sheet_name='notes',header=False)
    
    # figure
    if os.path.isfile(ped_loc):
        worksheet = writer.sheets['notes']
        worksheet.insert_image('D1',ped_loc)
    
    format_excel(a,writer,changes,rules)
    writer.save()

'''
make pandas dataframe
'''
def make_df(a):
    result = {'recessive':{},'dominant':{},'X':{},'relatives':{}}
    H = {}
    for k,v in a.iteritems():
        # make header
        if k in ['known','hpos']: continue
        if k in ['recessive','dominant','X']:
            header = ['symbol','retnet','pubmed_score','pubmed','protein_atlas','variants','consequence','filter','gnomad_af','kaviar_af','cadd_phred','cnvs','igv_check','gene_id','original_gene_id']
            if k == 'recessive':
                header = header[:7]+['gnomad_hom_af']+header[7:]
            else:
                header = header[:4]+['pLI']+header[4:]
            # add knowns
            header = header[:2] + a['known'] + header[2:]
            # add hpos
            header = header + ['%(name)s(%(id)s)' % i for i in a['hpos']]
        elif k == 'relatives':
            header = ['Name','project_associated_id','Affected','HPO','role']
    
        H[k] = header
        # make arrays
        for h in header:
            result[k][h] = make_array(h,v)

    # convert to dataframe
    for k in result:
        result[k] = pd.DataFrame(result[k])
        result[k] = result[k][H[k]]
    return result

'''
make array for making pandas df
'''
def make_array(h,v):
    # h is the header, and v is an array ( or a dict in case of relatives )
    ary = []
    if h == 'retnet':
        for i in v:
            if i['retnet']:
                ary.append( 'mode:%(mode)s, disease:%(disease)s' % i[h] )
            else:
                ary.append(None)
    elif h == 'pubmed':
        ary = [', '.join(i['pubmed']) for i in v]
    elif h == 'protein_atlas':
        for i in v:
            if i[h]: ary.append('; '.join(['[cell type: %(cell_type)s, level:%(level)s, reliability: %(reliability)s]' % j for j in i[h]]))
            else: ary.append(None)
    elif h == 'variants':
        # v_id1:cleaned_id:genotype:AAChange,v_id1..
        for i in v:
            this = ['%(cleaned_id)s:%(genotype)s:%(AAChange)s' % j for j in i[h]]
            ary.append(', '.join(this))
    elif h in ['filter','consequence','gnomad_af','kaviar_af','gnomad_hom_af','cadd_phred']:
        for i in v:
            this = [str(j[h]) for j in i['variants']]
            ary.append(', '.join(this))
    elif h == 'cnvs':
        for i in v:
            this = ['%(id)s:%(type)s:reads_observed:%(reads_observed)s:reads_expected:%(reads_expected)s:ratio:%(ratio)s:symbols:%(symbols)s' % j for j in i['cnv']]
            ary.append(', '.join(this))
    elif h in ['Name','project_associated_id','Affected','role']:
        for k,v in v.items():
            ary.append(v[h])
    elif h == 'HPO':
        for k,v in v.items():
            this = ['%(name)s(%(id)s)' % j for j in v[h]]
            ary.append(', '.join(this))
    elif '(' in h:
        # phenogenon
        id = h.split('(')[1][:-1]
        for i in v:
            this = [j['p_value'] for j in i['hpos'] if j['id'] == id][0]
            ary.append(this)
    else:
        ary = [i[h] for i in v]
    return ary

'''
pubmed batch
'''

'''
find the freaking PID, Title or Abstract no matter what!
'''
def find_item(obj, key):
    if key in obj:
        return obj[key]
    if isinstance(obj, dict):
        for k in obj:
            if isinstance(obj[k], dict):
                item = find_item(obj[k], key)
                if item is not None:
                    return item
            elif isinstance(obj[k], list):
                for i in obj[k]:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item
    elif isinstance(obj, list):
        for k in obj:
            if isinstance(k, dict):
                item = find_item(k, key)
                if item is not None:
                    return item
            elif isinstance(k, list):
                for i in k:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item

def pubmed_query(gene,keywords,lag=None,email='logust@yahoo.com'):
    Entrez.email = email
    reg = '\\b|\\b'.join(keywords)
    reg = '\\b' + reg + '\\b'
    reg = re.compile(reg, re.IGNORECASE)
    term = gene + ' AND (' + ' OR '.join(['"' + t + '"' + '[Title/Abstract]' for t in keywords]) + ')'

    if lag:
        lag = int(lag/3600/24) # convert it to days, at least 1 day
        lag = 1 if lag < 1 else lag
            # need to update
        attempt = 1
        while attempt <= 10:
            try:
                search_results = Entrez.read(Entrez.esearch(db='pubmed', term=term, reldate=lag, datetype='pdat', usehistory='y'))
                break
            except URLError as err:
                print ('!!URLError %s' % err)
                time.sleep(2)
                attempt += 1
    else:
        # just search
        attempt = 1
        while attempt <= 10:
            try:
                search_results = Entrez.read(Entrez.esearch(db='pubmed',retmax=50, term=term, usehistory='y'))
                break
            except URLError as err:
                print ('URLError: %s at line 249' % err)
                time.sleep(2)
                attempt += 1
            except RuntimeError as err:
                print ('Runtime error: %s at line 276' % err)
                time.sleep(2)
                attempt += 1
    # now done the search. let's get results
    count = int(search_results["Count"])
    print count
    results = {'results':[], 'total_score':0}
    # get search content
    if count:
        attempt = 1
        while attempt <= 10:
            try:
                handle = Entrez.efetch("pubmed",
                                       restart=0,
                                       retmax=50,
                                       retmode="xml",
                                       webenv=search_results['WebEnv'],
                                       query_key=search_results['QueryKey']
                                       )
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print('Received error from server %s' % err)
                else:
                    print('Something is wrong while efetch..')
                print('Attempt %i of 10' % attempt)
                attempt += 1
                time.sleep(5)
            except SocketError as err:
                print('Socket error')
                time.sleep(2)
            except URLError as err:
                print ('URLError')
                time.sleep(2)
        record = Entrez.read(handle)
        if record:
            # got something. let's do some calculation
            for r in record['PubmedArticle']:
                # calculate score
                score = 0
                pid = str(find_item(r, 'PMID'))
                abstract_list = find_item(r, 'AbstractText')
                # parse abstract
                abstract = ''
                if abstract_list:
                    for a in abstract_list:
                        if hasattr(a, 'attributes') and 'Label' in a.attributes:
                            abstract = abstract + '<b>' + a.attributes['Label'] + ': </b>'
                            abstract = abstract + a + '<br/>'
                        else:
                            abstract = abstract + a

                title = find_item(r, 'ArticleTitle')
                if title:
                    score = score + len(reg.findall(title))
                if abstract:
                    score = score + len(reg.findall(abstract))

                # add result to genes[gene_name]
                if score:
                    results['results'].append({
                        'id': pid,
                        'title': title,
                        'abstract': abstract,
                        'score': score
                    })
                    results['total_score'] = results['total_score'] + score
        results['results'] = sorted(results['results'], key=lambda k: k['score'], reverse=True)
    return results

class IRDC_variant(Variants.Variant):
    '''
    extra attributes. Store facts(annotations) in db by Variants.Variant
    '''
    def __init__(self,variant_id,db_conn=sqlite3.connect('irdc.db')):
        db_conn.text_factory = str
        Variants.Variant.__init__(self,db_conn,variant_id)
        self.filter = None
        self.patients = []
        self.genes = []
        self.most_severe_consequence = None

class IRDC_variants(Variants.Variants):
    def __init__(self,variant_ids,db_conn=sqlite3.connect('irdc.db'),varsome_key=None):
        db_conn.text_factory = str
        Variants.Variants.__init__(self,db_conn,variant_ids,varsome_key=varsome_key)

class IRDC_gene(Genes.Gene):
    '''
    attributes:
        gene_id
        pLI:
        mis_z:
        retnet:{
            in_retnet:1,
            mode:['r'],
            description:''
        }
        patients:[],
        pubmeds:[pubmed_ids]
    '''
    def __init__(self,gene_id,db_conn=sqlite3.connect('irdc.db')):
        db_conn.text_factory = str
        Genes.Gene.__init__(self,db_conn,gene_id)

class IRDC_genes(Genes.Genes):
    # a simple wrapper of Genes.Genes
    def __init__(self,gene_ids=[],db_conn=sqlite3.connect('irdc.db')):
        db_conn.text_factory = str
        Genes.Genes.__init__(self,db_conn,gene_ids)

class Patient:
    '''
    
    '''
    def __init__(self,id,options,G=None,known=None,retnet=None,protein_atlas=None):
        self.id = id
        self.options = options
        self.relatives = {}
        self.exome_rare_file_path = os.path.join(options['irdc_exome_folder'],'combinedAnalysis/point_mutations/latest/rare_variants')
        self.exome_cnv_file_path = options['irdc_exome_cnv_folder']
        self.wgs_rare_file_path = options['irdc_wgs_folder']
        self.G = G # a shared IRDC_genes object, to help translate symbols
        self.known = known
        self.retnet = retnet
        # mongo and mysql
        self.mongo_client = pymongo.MongoClient(host=options['mongo']['host'], port=options['mongo']['port'])
        self.mongo_db = self.mongo_client[options['mongo']['database']]
        self.mongo = self.mongo_db[options['mongo']['table']]
        mysqldb = mysql.connector.connect(**options['mysql'])
        cursor = mysqldb.cursor()
        # query mysql
        sql = '''SELECT * FROM patient JOIN project_relation ON patient.ID = linked_id WHERE Name = %s AND db_name = 'patient' AND project_id = 1'''
        cursor.execute(sql,(id,))
        self.mysql = sqlite_utils.dict_factory(cursor,cursor.fetchone())
        self.irdc_id = self.mysql['project_associated_id']
        # find rare and cnv files locations
        self.exome_rare_file, self.exome_cnv_file, self.exome_cnv_X_file, self.wgs_rare_file = self._populate_file_locations(self.mysql)
        '''
        find relatives and their Affected / HPO
        '''
        sql = '''SELECT Name, patient_id, role, Notes, Affected, HPO, SEX, project_associated_id FROM patient_family_relation
                JOIN patient ON patient_id = patient.ID JOIN project_relation ON patient_id = linked_id
                WHERE family_id = (SELECT family_id FROM patient_family_relation WHERE patient_id = %s) AND patient_id > 0
				AND db_name = 'patient' AND project_id = 1'''
        cursor.execute(sql,(self.mysql['linked_id'],))
        for i in cursor:
            temp = sqlite_utils.dict_factory(cursor,i)
            # find only irdc patient
            '''
            need to change it for generic scenarios
            '''
            if re.match('3\d+',temp['Name']):
                # turn Hpo,Notes back to data from json
                temp['HPO'] = json.loads(temp['HPO'])
                temp['Notes'] = json.loads(temp['Notes'])
                self.relatives[temp['Name']] = temp
    
        # find rare files and cnv files
        
        for k,v in self.relatives.iteritems():
            # get batch folder
            v['exome_rare_file'], v['exome_cnv_file'], v['exome_cnv_X_file'], v['wgs_rare_file'] = self._populate_file_locations(v)

    '''
    populate file locations. One for rare file, two for cnv (cnv and cnv_X) files
    '''
    def _populate_file_locations(self,v):
        batch = 'IRDC_'+v['project_associated_id'].split('_')[1]
        cnv_name = '_'.join(v['project_associated_id'].split('_')[2:])
        exome_rare = os.path.join(self.exome_rare_file_path, v['project_associated_id']+'.csv')
        exome_cnv = os.path.join(self.exome_cnv_file_path,batch,'multi_exons',cnv_name+'_sorted_unique.bam.cnv')
        exome_cnv_x = os.path.join(self.exome_cnv_file_path,batch,'multi_exons',cnv_name+'_sorted_unique.bam_X.cnv')
        wgs_rare = os.path.join(self.wgs_rare_file_path,'rare-VEP_OXF_%s-annotations.csv' % v['Name'])
        return (
            exome_rare,
            exome_cnv if os.path.isfile(exome_cnv) else None,
            exome_cnv_x if os.path.isfile(exome_cnv_x) else None,
            wgs_rare if os.path.isfile(wgs_rare) else None,
        )
    
    '''
    read csv 
    return {
        ENSG0000: {
            original_symbol:
            variants: [
                id
                genotype
                filter
                consequence
                gnomad
                kaviar
                cadd
                het_p:[who,who,who]
                hom_p:[who,who,who]
            ]
        }
    }

    '''
    def _get_wgs_rare_snp_genes(self, file):
        result = {}
        csvreader = csv.reader(open(file,'r'), delimiter=',', quotechar='"')
        header = []
        
        return
    '''
    read csv.
    return {
        ENSG0000: {
            original_symbol:
            variants: [
                id
                genotype
                filter
                consequence
            ]
        }
    }
    '''
    def _get_exome_rare_snp_genes(self, file):
        result = {}
        csvreader = csv.reader(open(file,'r'), delimiter=',', quotechar='"')
        header = []
        variants = []
        for row in csvreader:
            if not header:
                # read header
                header = row
                continue
            # convert row from array to dict with header
            row = _make_dict(row,header)
            # get gene_ids
            # separate genes on non words
            gene_ids = re.split(r'[^\w]\s*(?![^()]*\))', row['ensemblID'])
            v_id = '-'.join(row['signature'].split('_'))
            variants.append(v_id)
            genotype = row[self.irdc_id]
            # some rows have retired ensembl ids. need original_symbol to work it out
            original_symbol = row['HUGO.no.splice.info']
            gene_ids = set(gene_ids)
            for g in gene_ids:
                result[g] = result.get(g,{'variants':[],'original_symbol':original_symbol})
                result[g]['variants'].append({
                    'variant_id':v_id,
                    'genotype':'hom' if genotype[0] == '1' and genotype[2] == '1' else 'het',
                    'filter':row['FILTER'],
                    'consequence':row['ExonicFunc'] or row['Func'],
                    'AAChange':row['AAChange'],
                })
        # annotate variants
        V = IRDC_variants(variants,varsome_key=self.options['varsome_API_key'])
        V.cadd_file = self.options['cadd_file']
        kaviars = V.kaviar_af
        gnomads = V.gnomad
        cadd = V.cadd_phred
        cleaned = V.cleaned_variants
        for k,v in result.iteritems():
            for i in v['variants']:
                i['cleaned_id'] = cleaned[i['variant_id']]
                i['id'] = i['cleaned_id']
                i['gnomad'] = gnomads[i['variant_id']]
                i['cadd_phred'] = cadd[i['variant_id']]
                i['kaviar_af'] = kaviars[i['variant_id']]
        return result
    
    '''
    read cnv.
    return [
        chrom:
        start:
        end:
        reads_observed:
        reads_expected:
        ratio:
        type:
        genes:[]
    ]
    '''
    def _get_exome_cnv(self):
        # find relatives (it doesn't matter if they are affected)
        relatives = []
        if self.relatives:
            relatives = [v['project_associated_id'] for k,v in self.relatives.items()]
        # add itself
        relatives.append(self.irdc_id)
        relatives = set(relatives)
        # find batch number
        batch_no = self.irdc_id.split('_')[1]
        proband_df = pd.DataFrame()
        all_df = pd.DataFrame()
        # merge everything into a big df
        for folder in [i for i in os.listdir(self.exome_cnv_file_path) if 'batch' in i]:
            # if in the right folder? if not, append to all_df. if yes, store proband related cnvs in df, else in all_df
            if folder == batch_no:
                bingo = True
            else:
                bingo = False
            for file in os.listdir(os.path.join(self.exome_cnv_file_path,folder)):
                if not file.endswith('.csv'): continue
                f = os.path.join(self.exome_cnv_file_path,folder,file)
                # batch1 has strange ids. (not start with IRDC) remove them.
                this_df = pd.read_csv(f)
                this_df = this_df[this_df.apply(lambda x:x['sample'].startswith('IRDC'), axis=1)]
                # normalise this_df's sample field
                this_df['sample'] = this_df['sample'].apply(lambda x: x[:-18])
                g = this_df.groupby('id')['sample'].apply(list)
                
                # don't change this iteritems!!! this is a pandas series
                # all_df
                good_key = []
                for k,v in g.iteritems():
                    if len(set(v) - relatives) == 1:
                        good_key.append(k)
            
                all_df = all_df.append( this_df[this_df['id'].isin(good_key)] )
                # remove rows with not-unique ids (if not only duplicated in relatives), or not belong to proband
                if bingo:
                    good_key = []
                    for k,v in g.iteritems():
                        if not (set(v) - relatives):
                            good_key.append(k)
                    print('good_key',good_key)
                    this_df = this_df[this_df['id'].isin(good_key)]
                    proband_df = proband_df.append(this_df)
        all_df = all_df.reset_index(drop=True)
        proband_df = proband_df.reset_index(drop=True)
        return proband_df, all_df
    
    '''
    NOT FINISHED!!!!
    get unique cnv (excluding relatives)
    unique is 'batch' unique. 
    will also look for patients with overlapping cnvs
    '''
    @property
    def exome_unique_cnv(self):
        if getattr(self, '_exome_unique_cnv', None) is None:
            result = []
            proband_df, all_df = self._get_exome_cnv()
            sys.exit()
                    
        
            for i in self.exome_cnv:
                # search all the files
                bad = None
                if i['chrom'] != 'X':
                    bad = _get_exome_unique_cnv(i,cnv_files)
                else:
                    bad = _get_exome_unique_cnv(i,cnv_X_files)
                if not bad:
                    result.append(i)
            # convert the data into a dict of genes
            result2 = {}
            for i in result:
                for g in i['genes']:
                    result2[g] = result2.get(g,[])
                    result2[g].append(i)
            self._exome_unique_cnv = result2
        return self._exome_unique_cnv

    @property
    def exome_cnv(self):
        if getattr(self, '_exome_cnv', None) is None:
            if self.exome_cnv_file:
                self._exome_cnv = self._get_exome_cnv(self.exome_cnv_file) + self._get_exome_cnv(self.exome_cnv_X_file)
            else:
                self._exome_cnv = []
        return self._exome_cnv
    
    '''
    extract genes from the rare_variant file
    '''
    @property
    def exome_rare_snp_genes(self):
        if getattr(self, '_exome_rare_snp_genes', None) is None:
            rare_snp_genes = self._get_exome_rare_snp_genes(self.exome_rare_file)
            self._exome_rare_snp_genes = rare_snp_genes
        return self._exome_rare_snp_genes

class report:
    def __init__(self,options):
        self.options = options
        G = IRDC_genes()
        # shared bad_genes
        G._bad_genes = open(options['bad_genes'],'r').readline().split()
        self.G = G
        # relatives are done together with probands. push them here to avoid repetitive calculation
        self.done = []
    
    '''
    fields update rules
    '''
    @property
    def field_rules(self):
        if getattr(self, '_field_rules', None) is None:
            fields_to_check = self.options['fields_to_check']['fields']
            # fields_update_methods have some predefined helpers
            
            # populate missing fields in field_helpers with method None
            field_rules = {}
            for i in fields_to_check:
                field_rules[i] = field_helpers.get(i, None)

            # make gnomad / kaviar / pubmedscore / pLI methods
            for i in ['gnomad_af', 'gnomad_hom_af', 'kaviar_af', 'pubmed_score', 'pLI']:
                if i not in fields_to_check: continue
                field_rules[i] = field_helpers[i](fields_to_check[i])
            self._field_rules = field_rules
        return self._field_rules
        
    '''
    known lists
    translate symbols to ensembl ids
    '''
    @property
    def known(self):
        if getattr(self, '_known', None) is None:
            known = {}
            for f in self.options['known_genes_lists']:
                symbols = open(f,'r').readline().split()
                known[f] = self.G.symbols_to_ensemblIds(symbols).values()
            self._known = known
        return self._known
    
    '''
    retnet
    '''
    @property
    def retnet(self):
        if getattr(self, '_retnet', None) is None:
            retnet = json.load(open(self.options['retnet'],'r'))
            myd = self.G.symbols_to_ensemblIds(retnet.keys())
            for k in retnet.keys():
                # change symbols to ensembl ids
                if k not in myd: continue
                retnet[myd[k]] = retnet.pop(k)
            self._retnet = retnet
        return self._retnet
    
    '''
    protein_atlas
    '''
    @property
    def protein_atlas(self):
        if getattr(self, '_protein_atlas', None) is None:
            self._protein_atlas = parse_protein_atlas(self.options['protein_atlas']['file'],'retina')
        return self._protein_atlas
    '''
    run analysis for all the patients
    '''
    def run(self):
        for p in self.options['patients']:
            print '----doing----'
            print p
            if p in self.done:
                # already done, pass
                continue
            # analysis
            result = self.analysis(p)
            doing = [p]
            if result['relatives']:
                # has relative,find similarly affected p
                this_affected = [v['Affected'] for k,v in result['relatives'].iteritems() if k == p][0]
                doing = [k for k,v in result['relatives'].iteritems() if v['Affected'] == this_affected]
            # make df for analysis
            export = make_df(result)
            # write to excels for this and relatives
            self.done.extend(doing)
            for r in doing:
                # add Notes to export, and find the pedigree figure
                ped_loc = self.options['photo_folder']
                if r == p:
                    this_p = Patient(r,self.options)
                    # relatives:
                    ped_loc = os.path.join(ped_loc, 'pt%s' % this_p.mysql['linked_id'], '_pedigree.png')
                    notes = json.loads(this_p.mysql['Notes'])
                else:
                    this_p = [v for k,v in result['relatives'].iteritems() if k == r][0]
                    notes = this_p['Notes']
                    ped_loc = os.path.join(ped_loc, 'pt%s' % this_p['patient_id'], '_pedigree.png')
                # make notes a dict of arrays
                for k in notes:
                    notes[k] = [notes[k]]
                export['notes'] = pd.DataFrame(notes)
                # write to excel
                f = os.path.join(self.options['output_dir'],r+'.xlsx')
                # now need to check if target excel file exists.
                # if yes, compare and highlight.
                # if no, write
                changes = {}
                if os.path.isfile(f):
                    # add change message
                    msg = 'comparing with old file to highlight changes'
                    print(msg)
                    excel_data = pd.ExcelFile(f)
                    # multiple index for time and msg
                    now = time.strftime('%Y-%m-%d %H:%M')
                    for sheet in ['recessive','dominant','X']:
                        old_df = excel_data.parse(sheet)
                        changes[sheet] = Compare.compare_dfs(
                            df1 = old_df,
                            df2 = export[sheet],
                            key = 'gene_id',
                            fields = self.field_rules,
                        )
                        # add deleted rows back to new df
                        headers = list(export[sheet])
                        export[sheet] = export[sheet].append(old_df[old_df['gene_id'].isin(changes[sheet]['-'])])[headers].reset_index(drop=True)
                    
                        # and add non-overlapping column back. this is clumsy. should have a better way to do this
                        added_columns = set(list(old_df)) - set(list(export[sheet]))
                        for a in added_columns:
                            this_list = []
                            for g in list(export[sheet]['gene_id']):
                                if g not in list(old_df['gene_id']):
                                    this_list.append(None)
                                else:
                                    this_list.append(list(old_df[old_df['gene_id'] == g][a])[0])
                            export[sheet][a] = this_list
                        
                    # history df. first make dict, then convert to df, then merge with old history
                    history_df = {}
                    for k1,v1 in changes.items():
                        history_df[k1] = {}
                        for k2,v2 in v1.items():
                            if type(v2) is set: v2 = list(v2)
                            history_df[k1][(now,self.options['message'],k2)] = json.dumps(v2,indent=4)
                    history_df = pd.DataFrame.from_dict(history_df)
                    history_df.index.names = ['time','message','change']
                    if 'history' in excel_data.sheet_names:
                        old_history = excel_data.parse('history',index_col=[0,1,2])
                        history_df = history_df.append(old_history)
                    export['history'] = history_df
                    
                write_to_excel(export,f,ped_loc,changes,self.options['highlight_rules'])
        # write bad genes
        outf = open(self.options['bad_genes'],'w')
        outf.write(' '.join(set(self.G._bad_genes)))
        print 'All done'
    '''
    analysis of a patient

    this is the main analysis. will seperate into dominant/recessive/X
    filter variants based on provided filters in self.options
    group by genes
    segregate if relatives available, and mark the relatives as done (not repeat same analysis)
    use cnv info if available. only consider cnv not overlapped by not(affected sibs)
    hpo headers are a union of relative hpos

    result = {
        //family_history: this will be added when writing to excel
        relatives:{}
        dominant:[{recessive + pLI - gnomad_hom_af}]
        X:[{recessive + pLI - gnomad_hom_af}]
        recessive:[{
            gene_id:
            symbol:
            retnet:
            knowns:[names]
            known1:
            known2:
            pubmed:
            protein_atlas:
            hpos:[{id,name,p_value}]
            igv_check:
            variants:[{
                variant_id:
                AAChange:
                consequence:
                filter:
                genotype:
                gnomad_hom_af:
                gnomad_af:
                kaviar_af:
                cadd_phred:
            }]
            cnv:[{
                chrom:
                start:
                end:
                reads_observed:
                reads_expected
                ratio
                type
                genes
            }]
        }]
    }


    '''
    def analysis(self, patient_id):
        P = Patient(patient_id, self.options, self.G, self.known, self.retnet, self.protein_atlas)
        # make P of relatives
        relatives_P = [Patient(i, self.options, self.G) for i in P.relatives]
        # extract all genes and push them to GENES. Note that G is only for translation
        if self.options['cnv']['switch']:
            exome_unique_cnv = P.exome_unique_cnv
        else:
            exome_unique_cnv = {}
        GENES = IRDC_genes(P.exome_rare_snp_genes.keys() + exome_unique_cnv.keys())
        GENES._bad_genes = self.G._bad_genes
        #G._bad_genes.extend(GENES._bad_genes)
        # hpos
        # intersection of hpos of affected relatives
        hpos = set([i['id'] for i in json.loads(P.mysql['HPO'])])
        for k2,v2 in P.relatives.iteritems():
            if v2['Affected'] == 1:
                hpos = hpos & set([i['id'] for i in v2['HPO']])
        # turn hpos to an array of dict
        hpos = [i for i in json.loads(P.mysql['HPO']) if i['id'] in hpos]
        # since some terminal hpos in RD are not well calculated in phenogenon, add RD and rod-cone or code-rod to hpos if there is not one already
        hpo_ids = [i['id'] for i in hpos]
        # rod-cone
        hpos = _check_add_hpo({'id':'HP:0000510','name':'Rod-cone dystrophy'}, hpos)
        # cone-rod
        hpos = _check_add_hpo({'id':'HP:0000548','name':'Cone/cone-rod dystrophy'}, hpos)
        # rd
        if 'HP:0000556' not in hpo_ids:
            hpos.append({'id':'HP:0000556','name':'Retinal dystrophy'})
        # load the hpo data
        for h in hpos:
            file_name = os.path.join('../rd/hpo_gene',h['id']+'.json')
            if os.path.isfile(file_name):
                h['data'] = json.load(open(os.path.join('../rd/hpo_gene',h['id']+'.json'),'r'))
            else:
                h['data'] = {'data':{'unrelated':{'recessive':[],'dominant':[]}}}
        # early onset?
        onset = 'adult_onset'
        if set(hpo_ids) & set(EARLY_ONSET_HPOS):
            onset = 'early_onset'
        # initialise with relatives
        result = {
            'hpos':hpos,
            'known':P.known.keys(),
            'relatives':P.relatives,
            'dominant':[],
            'recessive':[],
            'X':[],
        }

        for k1,v1 in P.exome_rare_snp_genes.iteritems():
            # make a temp entry
            original_id = k1
            symbol = GENES.symbol.get(k1,None)
            pLI = GENES.pLI.get(k1,None)
            if not symbol:
                # highly likely a retired ensemblID
                # use the original symbol to find gene_id
                temp =  self.G.symbols_to_ensemblIds([v1['original_symbol']])
                if temp:
                    k1 = temp.values()[0]
                    temp_G = IRDC_gene(k1)
                    symbol = temp_G.symbol
                    pLI = temp_G.pLI
                else:
                    # can't do anything with this, such as ENSG00000211940: IGHV3-9
                    symbol = v1['original_symbol']
        
            this = {
                'gene_id':k1,
                'original_gene_id': original_id,
                'symbol': symbol,
                'protein_atlas': self.protein_atlas.get(k1,[]),
                'retnet': P.retnet.get(k1,{}),
            }
            # pubmed
            key = this['symbol']+P.options['mongo']['key'] if this['symbol'] else ''
            pubmed = P.mongo.find_one({'key':key}) if this['symbol'] else {'result':{'results':[],'total_score':0}}
            if not pubmed:
                # use scrutinise to search
                keywords = P.options['mongo']['key'][1:-1].split(',')
                print 'search pubmed for '+this['symbol']
                pubmed_result = pubmed_query(this['symbol'],keywords)
                pubmed = {'result':pubmed_result}
                # update database
                now = time.mktime(time.localtime())
                P.mongo.update({'key': key}, {'$set': {'result': pubmed_result, 'date': now}},upsert=True)
        
            this['pubmed_score'] = pubmed['result']['total_score']
            this['pubmed'] = [i['id'] for i in pubmed['result']['results']]
            # known
            for k2,v2 in P.known.iteritems():
                this[k2] = 'Y' if k1 in v2 else None
            # retnet
            this['retnet'] = P.retnet.get(k1,{})

            # check which group it belongs
            cnv = exome_unique_cnv.get(k1,[])

            # dominant?
            rare_variants = []
            for v in v1['variants']:
                bad = 0
                for p in relatives_P:
                    this_v = [z['cleaned_id'] for z in p.exome_rare_snp_genes.get(k1,{'variants':[]})['variants']]
                    # affected/unaffected relatives have it?
                    if p.mysql['Affected'] != P.mysql['Affected'] and v['cleaned_id'] in this_v:
                        bad = 1
                        break
                    if p.mysql['Affected'] == P.mysql['Affected'] and v['cleaned_id'] not in this_v:
                        bad = 1
                        break
                if bad: continue
                # gnomad?
                v['gnomad_af'] = parse_gnomad(v['gnomad'])
                if v['gnomad_af'] != None and v['gnomad_af'] <= self.options['cut_offs'][onset]['gnomad']:
                    rare_variants.append(v)
                elif v['gnomad_af'] == None and v['kaviar_af'] <= self.options['cut_offs'][onset]['kaviar']:
                    rare_variants.append(v)
            rare_cnv = []
            for c in cnv:
                bad = 0
                for p in relatives_P:
                    this_c = [z['id'] for z in exome_unique_cnv.get(k1,[])]
                    if p.mysql['Affected'] != P.mysql['Affected'] and z['id'] in this_c:
                        bad = 1
                        break
                    if p.mysql['Affected'] == P.mysql['Affected'] and z['id'] not in this_c:
                        bad = 1
                        break
                if bad: continue
                rare_cnv.append(c)

            if len(rare_variants) + len(rare_cnv):
                # add to dominant
                d_this = copy.copy(this)
                # add pLI
                d_this['pLI'] = pLI
                # add variants
                d_this['variants'] = rare_variants
                d_this['cnv'] = rare_cnv
                # add hpos
                d_this['hpos'] = []
                for h in hpos:
                    d_this['hpos'].append({
                        'id':h['id'],
                        'name':h['name'],
                        'p_value':([i['pass_p_val'] for i in h['data']['data']['unrelated']['dominant'] if i['gene_id'] == k1]+[None])[0] # add None to prevent indexing an empty array
                    })
                
                # any of the variants are close together?
                rare_ids = set([i['cleaned_id'] for i in rare_variants] + [i['id'] for i in cnv])
                igv_check = None
                for c in combinations(rare_ids,2):
                    start0 = start1 = end0 = end1 = None
                    if c[0] in [i['cleaned_id'] for i in rare_variants]:
                        start0, end0 = parse_pos(c[0],'snp')
                    else:
                        start0, end0 = parse_pos(c[0],'cnv')
                    if c[1] in [i['cleaned_id'] for i in rare_variants]:
                        start1, end1 = parse_pos(c[1],'snp')
                    else:
                        start1, end1 = parse_pos(c[1],'cnv')
                    igv_check = gap_check((start0,end0),(start1,end1),self.options['igv_check_size'])
                    if igv_check: break
                d_this['igv_check'] = igv_check

                result['dominant'].append(d_this)
                
                # is this an X-linked gene?
                if GENES.genomic_pos.get(k1,{'chr':None})['chr'] == 'X':
                    result['X'].append(d_this)

            # recessive?
            rare_variants = []
            # rare_cnv remains the same
            for v in v1['variants']:
                # gnomad hom? not check kaviar since it doesnt have hom af
                v['gnomad_af'] = parse_gnomad(v['gnomad'])
                v['gnomad_hom_af'] = parse_gnomad(v['gnomad'],hom=True)
                if v['gnomad_hom_af'] == None or v['gnomad_hom_af'] <= self.options['cut_offs'][onset]['gnomad']:
                    rare_variants.append(v)
            rare_ids = set([i['cleaned_id'] for i in rare_variants] + [i['id'] for i in cnv])

            # intersect with affected
            for p in [i for i in relatives_P if i.mysql['Affected'] == P.mysql['Affected']]:
                this_rare_id = [j['cleaned_id'] for j in p.exome_rare_snp_genes.get(k1,{'variants':[]})['variants']] + [j['id'] for j in exome_unique_cnv.get(k1,[])]
                rare_ids = rare_ids & set(this_rare_id)
            # more than 2? hom * 2 + het * 1. count cnvs as hom based on options['cnv']
            count = (len([i for i in rare_variants if i['cleaned_id'] in rare_ids and i['genotype'] == 'hom']) + len([i for i in rare_cnv if i['id'] in rare_ids and i['genotype'] == 'hom'])) * 2 + len([i for i in rare_variants if i['cleaned_id'] in rare_ids and i['genotype'] == 'het']) + len([i for i in rare_cnv if i['id'] in rare_ids and i['genotype'] == 'het'])
            if count > 1:
                # great, let's check unaffected. unaffected can't exaust all possible combinations of 2 or more than 2 of rare_ids
                unaffected = [i for i in relatives_P if i.mysql['Affected'] != P.mysql['Affected']]
                # expand hom * 2
                extra = []
                for i in rare_variants+rare_cnv:
                    if i['id'] in rare_ids and i['genotype'] == 'hom':
                        extra.append(i['cleaned_id'])
                rare_ids = list(rare_ids) + extra

                good = None
                if unaffected:
                    # combo with 2
                    combo = combinations(sorted(rare_ids),2)
                    relatives_combo = []
                    for p in unaffected:
                        pv = p.exome_rare_snp_genes.get(k1,{'variants':[]})['variants']
                        pc = exome_unique_cnv.get(k1,[])
                        ids = []
                        for i in pv+pc:
                            if i['id'] in rare_ids:
                                if i['genotype'] == 'hom':
                                    ids = ids + [i['id'] + i['id']]
                                else:
                                    ids.append(i['id'])

                        relatives_combo.extend(combinations(sorted(ids),2))
                    good = set(combo) - set(relatives_combo)
            
                    if not good: continue
                    # flatten good
                    rare_ids = reduce((lambda x, y: x + y), good)
                r_this = copy.copy(this)
                # add variants, only the good ones
                r_this['variants'] = [i for i in rare_variants if i['cleaned_id'] in rare_ids]
                for i in r_this['variants']:
                    i['gnomad'] = None
                r_this['cnv'] = [ i for i in rare_cnv if i['id'] in rare_ids]
                
                # any of the variants are close together?
                igv_check = None
                for c in combinations(rare_ids,2):
                    start0 = start1 = end0 = end1 = None
                    if c[0] in [i['cleaned_id'] for i in rare_variants]:
                        start0, end0 = parse_pos(c[0],'snp')
                    else:
                        start0, end0 = parse_pos(c[0],'cnv')
                    if c[1] in [i['cleaned_id'] for i in rare_variants]:
                        start1, end1 = parse_pos(c[1],'snp')
                    else:
                        start1, end1 = parse_pos(c[1],'cnv')
                    igv_check = gap_check((start0,end0),(start1,end1),self.options['igv_check_size'])
                    if igv_check: break
                r_this['igv_check'] = igv_check
                
                # add hpos
                r_this['hpos'] = []
                for h in hpos:
                    r_this['hpos'].append({
                        'id':h['id'],
                        'name':h['name'],
                        'p_value':([i['pass_p_val'] for i in h['data']['data']['unrelated']['recessive'] if i['gene_id'] == k1]+[None])[0] # add None to prevent indexing an empty array
                    })
                result['recessive'].append(r_this)
        self.G._bad_genes.extend(GENES._bad_genes)
        return result
