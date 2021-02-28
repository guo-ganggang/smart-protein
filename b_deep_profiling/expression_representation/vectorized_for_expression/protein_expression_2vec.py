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

# Put together what Largevis predicted
def generate_dimension_reduction_vector(dir_path):

    human_protein_symbols = []
    train_gene_expression_in_all_tissues_fp = dir_path + 'gene_expression_in_all_tissues_by_total_tpm_values.csv'
    with codecs.open(train_gene_expression_in_all_tissues_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if temp[0] == 'gene_symbol':
                continue
            human_protein_symbols.append(temp[0].strip())
    print('human_protein_symbols: ' + str(len(human_protein_symbols)))

    train_protein_symbols = set()
    train_gene_expression_in_all_tissues_fp = dir_path + 'train_gene_expression_in_all_tissues_by_total_tpm_values.csv'
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
    generate_dimension_reduction_vector(dir_path)