#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 6/11/2018 5:20 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : protein_semsim_2vec.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import data_helper as dh
import csv

def generate_protein_ppi_vector(dir_path,norm = 'non-norm'):

    print("Loading term vector")
    # vector_values = dh.mapping_gene_id_symbol_by_dbs(dir_path)
    vector_values = dh.prepare_ppi_scores_for_visualize(dir_path)

    with open('train_protein_symbol_ppi_vector.csv', 'w') as f:
        for key in vector_values.keys():
            f.write("%s,%s\n" % (key, vector_values[key]))


if __name__ == '__main__':

    dir_path = ''
    generate_protein_ppi_vector(dir_path)