'''
Given a gene, find all individuals carrying rare variants (given the rare folder path)
Need to mount gnomad drive first
Then run the script without --cadd arg, to get a vcf for CADD scoring (I do it online).
After getting the cadd.tsv, feed it to --cadd and run the script again to produce the final report
'''
import os
import csv
import argparse
import sys
from collections import Counter
sys.path.append('../BioTools')
import CommonFuncs
import gnomad_utils


def main(params):
    # genotype dict
    genotype_dict = {
        1: 'het',
        2: 'hom'
    }
    # read patient info
    patient = {}
    patient_header = []
    if not params.cadd:
        # no cadd provided, write to vcf
        outfile = '{}.vcf'.format(params.gene)
    else:
        outfile = '{}.txt'.format(params.gene)
    with open(PATIENT_CSV, 'rt', encoding='utf-8-sig') as inf:
        csvreader = csv.reader(inf)
        for row in csvreader:
            row = row[:11]
            if not patient_header:
                patient_header = row
                continue
            record = dict(zip(patient_header, row))
            del record['IRDC ID']
            patient[row[0]] = record

    variants = set()
    report = {}
    for csvfile in os.listdir(PATH_TO_CSVS):
        if not csvfile.endswith('.csv'):
            continue
        header = []
        genotype_header = None
        with open(os.path.join(PATH_TO_CSVS, csvfile), 'rt') as inf:
            csvreader = csv.reader(inf)
            for row in csvreader:
                if not header:
                    header = row
                    continue
                record = dict(zip(header, row))
                # has gene?
                genes = record['HUGO.no.splice.info'].split(',')
                if params.gene in genes:
                    variant = CommonFuncs.find_leftmost_synonymous_variant(
                        CommonFuncs.clean_variant(record['clean.signature'].replace('_', '-')))
                    variants.add(variant)
                    sample = csvfile.split('.csv')[0]
                    genotype = genotype_dict.get(
                        Counter(record[sample].split(':')[0])['1'], 'unknown')
                    if variant not in report:
                        report[variant] = record
                        report[variant]['samples'] = [
                            {'id': sample, 'genotype': genotype}]
                    else:
                        report[variant]['samples'].append(
                            {'id': sample, 'genotype': genotype})
    if not params.cadd:
        # sort and write vcf
        with open(outfile, 'wt') as outf:
            # write vcf header
            outf.write('##VCF4.1\n')
            outf.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT']) + '\n')
            for variant in sorted(list(variants), key=lambda x: int(x.split('-')[1])):
                chrom, pos, ref, alt = variant.split('-')
                row = [chrom, pos, '.', ref, alt]
                outf.write('\t'.join(row) + '\n')
    else:
        # get gnomads
        gnomads = gnomad_utils.overall_freqs(list(variants), PATH_TO_GNOMAD)
        # get cadd
        cadds = {}
        with open(params.cadd, 'rt') as inf:
            for line in inf:
                if line.startswith('#'):
                    continue
                row = line.rstrip().split('\t')
                cadds['-'.join(row[:4])] = row[-1]
        # write report
        with open(outfile, 'wt') as outf:
            for variant in sorted(list(variants), key=lambda x: int(x.split('-')[1])):
                outf.write(variant + ':\n')
                outf.write('\tFilter: {}\n'.format(report[variant]['FILTER']))
                outf.write('\t{}\n'.format(report[variant]['AAChange']))
                outf.write('\tPolyphen: {}, SIFT: {}, MutationTaster: {}\n'.format(
                    report[variant]['LJB_PolyPhen2_Pred'], report[variant]['LJB_SIFT_Pred'], report[variant]['LJB_MutationTaster_Pred']))
                outf.write('\tgnomad_af:{}, gnomad_hom_f:{}, cadd:{}\n'.format(
                    gnomads[variant]['gnomad_af'], gnomads[variant]['gnomad_hom_f'], cadds[variant]))
                for sample in report[variant]['samples']:
                    outf.write('\t{} ({}):\n'.format(
                        sample['id'], sample['genotype']))
                    study_number = sample['id'].split('_')[3]
                    if study_number in patient:
                        for h in patient_header:
                            if h in patient[study_number]:
                                outf.write('\t\t{}: {}\n'.format(
                                    h, patient[study_number][h]))
                    outf.write('\n')


if __name__ == '__main__':
    PATH_TO_CSVS = '/Users/logust/Dropbox (UKIRDC)/UKIRDC Team Folder/UKIRDC_March_2018/rare_variants'
    # note that some of the patient info might be missing in the file
    PATIENT_CSV = 'data/UKIRDC_Analysis_Summary.csv'
    PATH_TO_GNOMAD = 'data/gnomad'
    parser = argparse.ArgumentParser(
        description='Given gene, get info. If cadd tsv is given, produce final report. otherwise produce vcf for getting cadd tsv')

    parser.add_argument('--gene', dest='gene',
                        help='gene symbol (ABCA4)')

    parser.add_argument('--cadd', dest='cadd',
                        help='cadd scores for all variants')
    args = parser.parse_args()
    main(args)
    print('done')
