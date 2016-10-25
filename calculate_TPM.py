# coding: utf-8

### LW 10-10-16 ##
## This code reads in a .counts file developed through the Rsubread R package (RNA-seq data) and calculates the TPM
## USAGE: in command line type: python calculate_TPM.py -gtf [path to input gtf file] -counts [path to input counts file (csv)
## -o [output file of genes/lengths]-O [output file of TPM calculations] -f [featurename in gtf file, ex: 'exon', 'gene']



import pandas as pd
import re
from collections import defaultdict
import csv
import sys
import argparse

def check_arg():
    parser = argparse.ArgumentParser(description='Process gff/gtf file, counts file to obtain TPM values per gene')
    parser.add_argument('-gtf', '--input_file', 
                        help = 'path to input gtf file',
                        required=True )
    parser.add_argument('-counts', '--counts_file',
                        help='path to input counts file, must contain a column \n with gene names and contain column header of "genes" and must be comma delimited',
                        required=True)
    parser.add_argument('-o', '--output1', help='path to output file containing genes and lengths', required=True)
    parser.add_argument('-O', '--output2', help='path to output file containing TPM calculations', required=True)
    parser.add_argument('-f', '--featureName', help='feature used to calculate gene length, default = "gene"', required=False, default='gene')
    args = parser.parse_args()
    featureName = args.featureName
    in_file = args.input_file
    gene_lengths_file = args.output1
    counts = args.counts_file
    TPM_file = args.output2
    return featureName, in_file, gene_lengths_file, counts, TPM_file

def get_gene_info(filename, featureName):
    x = []
    y = []
    with open(filename, 'r') as f:
        reader = f.readlines()[7:]
        for line in reader:
            gene_name = re.compile(r'Name=([a-z]{2,}[A-Z]{0,1})')
            rr_gene = re.compile(r'rr[l|s][a-zA-z]{1}')
            row = line.split() 
            if len(row) == 9:
                if row[2] == featureName:
                    gene_match = gene_name.match(row[8].split(';')[2])
                    if gene_match:
                        gene = gene_match.group(1)
                        ribo = rr_gene.match(gene)
                        start = float(row[3])
                        end = float(row[4])
                        length = end - start
                        if ribo:
                            pass
                        else:
                            x.append(gene)
                            y.append(length)
    return x, y
            



def format_dataframes(genes, lengths, counts):
    genes_lengths = {'genes' : genes, 'lengths':lengths}
    df = pd.DataFrame.from_dict(genes_lengths)
    counts_df = pd.read_table(counts, sep=',', header=0)
    df_sorted = df.sort_values('genes', ascending=True)
    df = df_sorted.copy()
    counts_df_sorted = counts_df.sort_values('Gene', ascending=True)
    counts = counts_df_sorted.copy()
    unique_counts = counts_df_sorted.drop_duplicates('Gene')
    unique_lengths = df.drop_duplicates('genes')
    df = df.reset_index(drop=True)
    unique_counts['lengths'] = unique_lengths['lengths'].values
    columns = list(unique_counts)
    return columns, unique_counts
    
def create_new_df(df, columns):
    new_counts = df.copy()
    for i in range(len(columns)): 
        if columns[i] == 'Gene':
            continue
        if columns[i] == 'lengths':
            continue
        else:
            col_RPK = columns[i]+'_RPK'
            col_RPKM = columns[i]+'_RPKM'
            col_TPM = columns[i]+'_TPM'
            sample =  list(new_counts[columns[i]].values)
            lengths = list(new_counts['lengths'].values)
            rpk = [float(s) / float(l) for s,l in zip(sample, lengths)] # Divide the counts for each gene by the length of each gene
            new_counts[col_RPK] = rpk
            rpk_column = list(new_counts[col_RPK].values)
            RPKM = sum(new_counts[col_RPK])/1000000 # Take the sum of the counts for each sample and divide by 1000000 to get the RPKM
            TPM = [r/RPKM for r in rpk_column] ## divide the rpk by the sum of the counts/1000000 to get the TPM
            new_counts[col_TPM] = TPM
    return new_counts








if __name__ =='__main__': #this code wil only be executed when run directly, not imported as a module
    featureName, in_file, gene_lengths_file, counts, TPM_file = check_arg()
    
    print '________________________________FILES_____________________________________'
    print '{} {}'.format('input file (gtf):', in_file)
    print '{} {}'.format('gene/length output file (csv):', gene_lengths_file)
    print '{} {}'.format('TPM output file (csv):', TPM_file)
    print '{} {}'.format('counts input file (csv):', counts)
    print '{} {}'.format('feature used for gene lengths:', featureName)
    print '______________________________HERE WE GO__________________________________'
    
    genes, lengths = get_gene_info(in_file, featureName)
    columns, df = format_dataframes(genes, lengths, counts)
    df.to_csv(gene_lengths_file)  
    c3 = create_new_df(df, columns)  
    c3.to_csv(TPM_file)
    
    print 'ALL DONE!'