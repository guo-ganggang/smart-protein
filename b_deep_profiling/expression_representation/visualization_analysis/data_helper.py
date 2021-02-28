#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 17/11/2018 9:13 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : data_helper.py
# @Software: Mining from a specific set of proteins in human sperm

import codecs
import os

import sys

# reload(sys)
# sys.setdefaultencoding('utf8')

def read_protein_atlas_info(dir_path):

    protein_expression_vector_class_fname = dir_path + 'human_protein_expression_vector_class.csv'
    vector_class = {}
    if not os.path.exists(protein_expression_vector_class_fname) or (os.path.getsize(protein_expression_vector_class_fname) == 0):
        human_protein_atlas_pf = dir_path + 'raw_data/the_human_protein_atlas/proteinatlas.tsv'
        protein_symbol_chromosome = {}
        chromosome_codes = []
        protein_symbols = set()
        for i in range(1,23):
            chromosome_codes.append(str(i))
        chromosome_codes.append('X')
        chromosome_codes.append('Y')
        with codecs.open(human_protein_atlas_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('	')
                if len(temp) < 5:
                    continue
                chromosome = temp[4].strip()
                if chromosome not in chromosome_codes:
                    chromosome = temp[3].strip()
                    if chromosome not in chromosome_codes:
                        continue
                if chromosome == 'X':
                    protein_symbol_chromosome[temp[0].strip()] = 23
                elif chromosome == 'Y':
                    protein_symbol_chromosome[temp[0].strip()] = 24
                else:
                    protein_symbol_chromosome[temp[0].strip()] = float(chromosome)
        print('protein_symbol_chromosome: ' + str(len(protein_symbol_chromosome.keys())))
        related_protein_symbols = set(protein_symbol_chromosome.keys())

        # 基因名称与其表达向量
        gene_expression_in_all_tissues_pf = dir_path + 'raw_data/gene_expression_in_all_tissues_by_total_tpm_values.csv'
        human_gene_expression_vector_pf = dir_path + 'LargeVis-master/Examples/protein_expression/human_protein_expression_128D.csv'

        gene_names = []
        with codecs.open(gene_expression_in_all_tissues_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if temp[0].strip() == 'gene_symbol':
                    continue
                if len(temp) < 2:
                    continue
                gene_names.append(temp[0].strip())
        print('gene_names: ' + str(len(gene_names)))

        list_count = 0
        with open(protein_expression_vector_class_fname, 'w') as output_file:
            with codecs.open(human_gene_expression_vector_pf, "rb", "utf-8") as input_file:
                for line in input_file:
                    temp = line.strip().split(' ')
                    if len(temp) != 128:
                        continue
                    gene_name = gene_names[list_count]
                    list_count += 1
                    # print(gene_name)
                    if gene_name not in related_protein_symbols:
                        continue
                    vector_class[','.join(temp[1:])] = protein_symbol_chromosome[gene_name]
                    output_file.write(','.join(temp[1:]) + '$' + str(protein_symbol_chromosome[gene_name]) + '\n')
    else:
        with codecs.open(protein_expression_vector_class_fname, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if len(temp) < 2:
                    continue
                vector_class[temp[0].strip()] = float(temp[1].strip())

    print('vector_class: ' + str(len(vector_class.keys())))
    return vector_class

def temp_statistic_gene_info(dir_path):

    human_protein_atlas_pf = dir_path + 'raw_data/the_human_protein_atlas/proteinatlas.tsv'
    classes_protein = set()
    with codecs.open(human_protein_atlas_pf, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) < 7:
                continue
            protein_class = temp[6].strip()
            if 'proteins' not in protein_class:
                continue
            classes_protein.add(protein_class)

    print('classes_protein: ' + str(len(classes_protein)))

    for conten in classes_protein:
        print(conten)




if __name__ == '__main__':

    dir_path = ''
    # read_protein_atlas_info(dir_path)
    temp_statistic_gene_info(dir_path)

