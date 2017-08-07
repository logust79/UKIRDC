from __future__ import print_function, division
import sys
import tabix
import os

path = '/media/jing/SEAGATE/db/gnomAD'
release = 'release-170228'

'''
coverage
'''
def coverage(v,path_to_gnomad,mode='exome'):
    # pytabix does not support header yet. hard code it
    header = ['chrom','pos','mean','median',1,5,10,15,20,25,30,50,100,]
    chrom,pos,ref,alt = v.split('-')
    if mode == 'exome':
        file = os.path.join(path_to_gnomad,'exomes','coverage','exacv2.chr'+chrom+'.cov.txt.gz')
    elif mode == 'genome':
        file = os.path.join(path_to_gnomad,'genomes','coverage','gnomad.chr'+chrom+'.cov.txt.gz')
    else:
        msg = "mode only accepts 'exome' or 'genome'"
        raise ValueError(msg)
    tb = tabix.open(file)
    r = tb.query(chrom, int(pos)-1, int(pos))
    r = list(r)
    if not r:
        # not covered
        return None
    r = r[0]
    return {a:b for a,b in zip(header,r)}

'''
exome freqs
'''
def freqs(v,path_to_gnomad,mode='exome'):
    # pytabix does not support header yet. hard code it
    header = ['chrom','pos','id','ref','alt','quality','filter','info']
    chrom,pos,ref,alt = v.split('-')
    if mode == 'exome':
        file = os.path.join(path_to_gnomad,'exomes','vcf','gnomad.exomes.r2.0.1.sites.vcf.gz')
    elif mode == 'genome':
        file = os.path.join(path_to_gnomad,'genomes','vcf','gnomad.genomes.r2.0.1.sites.'+chrom+'.vcf.gz')
    tb = tabix.open(file)
    records = tb.query(chrom, int(pos)-1, int(pos))

    found = False
    for r in records:
        if not r: return None
        data = {a:b for a,b in zip(header,r)}
        # ref ok?
        if data['ref'] != ref:
            continue

        found = True
        # find alt
        g_alts = data['alt'].split(',')
        if alt not in g_alts:
            return None
        alt_ind = g_alts.index(alt)
        
        # parse info
        # no need for annotation?

        info = data['info'].split(';CSQ=A')[0] # 1 for annotation
        info = info.split(';')
        info_dict = {}
        for i in info:
            if not '=' in i: continue
            a,b = i.split('=')
            b = b.split(',')
            ind = min(alt_ind,len(b)-1)
            b = b[ind]
            # convert to number if possible
            try:
                if '.' in b:
                    b = float(b)
                else:
                    b = int(b)
            except ValueError:
                pass
            info_dict[a] = b
        info_dict['FILTER'] = info_dict['AS_FilterStatus']
        return info_dict
                

    if not found: return None

'''
genome coverage
'''
def genome_coverage(v):
    # pytabix does not support header yet. hard code it
    header = ['chrom','pos','mean','median',1,5,10,15,20,25,30,50,100,]
    chrom,pos,ref,alt = v.split('-')
    file = os.path.join(path,release,'genomes','coverage','gnomad.chr'+chrom+'.cov.txt.gz')
    tb = tabix.open(file)
    r = tb.query(chrom, int(pos)-1, int(pos))
    r = list(r)
    if not r:
        # not covered
        return None
    r = r[0]
    return {a:b for a,b in zip(header,r)}
