'''
helper functions to highlight changes
'''
from __future__ import print_function, division
import pandas as pd
import numpy as np

'''
convert np.nan to y
'''
def nan_to_(x,y):
    return [y if pd.isnull(i) else i for i in x]
'''
retnet, check mode only
'''
def retnet_field_helper(A,B):
    # convert np.nan to None
    A,B = nan_to_([A,B],None)
    # if one is empty, the other is not?
    if bool(A) ^ bool(B):
        return True
    if not A:
        return False
    def parse_mode(x):
        return x.split(', ')[0].split(':')[1]
    if parse_mode(A) != parse_mode(B):
        return True
    else:
        return False

'''
pubmed_score
pass in a fold to trigger highlighting
if fold = 2, new score has to be at least 2 times the old one to trigger highlighting
'''
def pubmed_score_field_helper_factory(fold):
    def inner_cb(A,B):
        # convert np.nan to None
        A,B = nan_to_([A,B],0)

        if B > A * fold:
            return True
        else:
            return False

    return inner_cb

'''
variants, simply check variant ids
'''
def variants_field_helper(A,B):
    def parse_variants(x):
        return {i.split(':')[0] for i in x.split(', ')}
    
    if parse_variants(A) == parse_variants(B):
        return False
    else:
        return True

'''
exac / kaviar
'''
def ek_field_helper_factory(cutoff):
    # there might be a ',' in the cell
    def parse(x):
        return [float(i) for i in x.split(', ')]
    def inner_cb(A,B):
        # convert to string
        A,B = str(A), str(B)
        # equal?
        if A == B:
            return False
        
        A,B = parse(A), parse(B)
        if len(A) != len(B):
            return True
        for a,b in zip(A,B):
            # one > thrd, one < thrd?
            if sorted([a,b,cutoff])[1] == cutoff: return True
            # both < thrd?
            if max(a,b,cutoff) == cutoff: return True
        return False

    return inner_cb

'''
cnvs, checks id and type
'''
def cnvs_field_helper(A,B):
    def parse_cnvs(x):
        return {i[0]+i[1] for i in x.split(', ')}
       
    A,B = nan_to_([A,B],None)
    if bool(A) ^ bool(B):
        return True
    if not A:
        return False
    if parse_cnvs(A) != parse_cnvs(B):
        return True
    else:
        return False

'''
pLI
'''
def pLI_field_helper_factory(change):
    def inner_cb(A,B):
        A,B = nan_to_([A,B],None)
        if bool(A) ^ bool(B):
            return True
        if not A:
            return False

        if abs(B-A) >= change:
            return True
        else:
            return False
            
    return inner_cb

field_helpers = {
    'retnet': retnet_field_helper,
    'pubmed_score': pubmed_score_field_helper_factory,
    'variants': variants_field_helper,
    'exac_hom_af': ek_field_helper_factory,
    'exac_af': ek_field_helper_factory,
    'kaviar_af': ek_field_helper_factory,
    'cnvs': cnvs_field_helper,
    'pLI': pLI_field_helper_factory,
};
