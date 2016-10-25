
# coding: utf-8

## LW 9-2016 
## takes a gtf file, gets a list of gene names and corresponding start positions, returns dataframe#
## note: this will return duplicate gene names and positions for a gene, need to sort out##
## Usage: python assign_ref.py <gff_file> <file_to_add_ref_to> <out_file>

import pandas as pd
import numpy as np
import re
import sys



GFF = sys.argv[1]
COUNTS = sys.argv[2]
OUT = sys.argv[3]

print '{} {}'.format('GFF file:', GFF)
print '{} {}'.format('Counts file (csv):', COUNTS)
print '{} {}'.format('Outfile:', OUT)

def get_gene_info(filename):
    with open(filename, 'r') as f:
        reader = pd.read_table(f, sep='\t', usecols=[0,3,4,8], header=None)
        gene_name = reader[8]
        gene_pos = reader[3]
        genes = []
        for thinger in gene_name:
            item = thinger.split(";")[0]
            m = re.search(r'"(.*)"', item)
            try:
                genes.append(m.group(0).strip(r'""'))
            except AttributeError:
                genes.append('error')
        df1 = {'gene_name':genes, 'position':gene_pos}
        df2 = pd.DataFrame(data=df1)
        return df2
    

def generate_positions(df):
    positions =[] #set-up empty list of positions
    current_gene = df['gene_name'][0] #start current_ref at the beginning
    for key,value in enumerate(df['gene_name']):    
        if current_gene == value: #if current_gene is equal to value of RefPos
            positions.append(df['position'][key]) #append to list called positions                  
        else:    #if not equal, yield the current_ref and IPD list to function 
            yield current_gene, positions
            del positions[:] #delete contents of IPD lis
            positions.append(df['position'][key]) #but since we have moved to next line, attach that to new IPD
            current_gene = value  #and make sure current_ref value is next RefPos value
    yield current_gene, positions  #this yields last IPD list when eof is reached        
            


def create_gene_pos_df():
    gene = []
    pos = []
    for a,b in generate_positions(gene_info):
        gene.append(a)
        pos.append(min(b))
    df = {'gene': gene,
          'position': pos}
    df2 = pd.DataFrame(df)
    return df2



def write_positions(filename, df2):
    with open(filename, 'r') as f:
        table = pd.read_table(f, sep=',')
        gene_name = table['gene']
        table['position'] = 0
        for i, k in enumerate(gene_name):
            if 'ins' in k:
                pass        
            try:
                pos = df2[df2['gene'] == k]['position'].iloc[0]
                table.loc[i,'position']= pos
            except IndexError:
                pass
        pd.DataFrame.to_csv(table, OUT)

if __name__ == '__main__':
    gene_info = get_gene_info(GFF)  
    uni_genes = gene_info['gene_name'].unique()
    #list_positions = generate_positions(uni_genes)
    write_positions(COUNTS, create_gene_pos_df())
    