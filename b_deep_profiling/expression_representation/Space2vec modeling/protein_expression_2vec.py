#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 1/1/2019 6:01 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : protein_expression_2vec.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import codecs
import gene_expression as ge

# TPM was expressed in all samples of all tissues for each protein
def prepare_gene_expression_for_extracting_vector(dir_path):
    raw_data_path = dir_path + 'raw_data/'
    ge.processing_GTEx_expression_data(raw_data_path)

# Preparing Data
def expression_dimensionality_reduction_vector(dir_path):
    human_gene_expression_in_all_tissues_fp = dir_path + 'raw_data/gene_expression_in_all_tissues_by_total_tpm_values.csv'
    human_gene_expression_in_all_tissues_for_dr_fp = dir_path + 'raw_data/human_protein_expression_265D.csv'
    with open(human_gene_expression_in_all_tissues_for_dr_fp, 'w') as output_file:
        output_file.write('20131' + ' ' + '265' + '\n')
        with codecs.open(human_gene_expression_in_all_tissues_fp, 'rb', 'utf-8') as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if temp[0].strip() == 'gene_symbol':
                    continue
                if len(temp) != 54:
                    continue
                temp_temp = ','.join(temp[1:])
                expression_values = ' '.join(temp_temp.split(','))
                output_file.write(expression_values+'\n')

# Put together what Largevis predicted
def generate_dimension_reduction_vector(dir_path):

    human_protein_symbols = []
    train_gene_expression_in_all_tissues_fp = dir_path + 'raw_data/gene_expression_in_all_tissues_by_total_tpm_values.csv'
    with codecs.open(train_gene_expression_in_all_tissues_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if temp[0] == 'gene_symbol':
                continue
            human_protein_symbols.append(temp[0].strip())
    print('human_protein_symbols: ' + str(len(human_protein_symbols)))

    train_protein_symbols = set()
    train_gene_expression_in_all_tissues_fp = dir_path + 'raw_data/train_gene_expression_in_all_tissues_by_total_tpm_values.csv'
    with codecs.open(train_gene_expression_in_all_tissues_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            train_protein_symbols.add(temp[0].strip())
    print('train_protein_symbols: ' + str(len(train_protein_symbols)))

    human_gene_expression_128D_fp = dir_path + 'LargeVis-master/Examples/protein_expression/human_protein_expression_128D.csv'
    train_protein_gene_expression_128D_vector_fp = dir_path + 'train_protein_gene_expression_128D_vector.csv'
    protein_symbol_index = 0
    with open(train_protein_gene_expression_128D_vector_fp, 'w') as output_file:
        with codecs.open(human_gene_expression_128D_fp, 'rb', 'utf-8') as input_file:
            for line in input_file:
                temp = line.strip().split(' ')
                if len(temp) != 128:
                    print(line.strip())
                    continue
                expression_values = ','.join(temp)
                protein_symbol = human_protein_symbols[protein_symbol_index]
                protein_symbol_index += 1
                if protein_symbol not in train_protein_symbols:
                    continue
                output_file.write(protein_symbol + '$' + expression_values + '\n')
    print('protein_symbol_index: ' + str(protein_symbol_index))

if __name__ == '__main__':

    dir_path = ''
    prepare_gene_expression_for_extracting_vector(dir_path)
    # expression_dimensionality_reduction_vector(dir_path)
    # generate_dimension_reduction_vector(dir_path)