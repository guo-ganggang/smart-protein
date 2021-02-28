#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 17/11/2018 9:13 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : data_helper.py
# @Software: Mining from a specific set of proteins in human sperm

import codecs
from itertools import islice
import numpy as np
from os import listdir
import re
import os
import nltk
import nltk.data
from nltk.corpus import stopwords
import json
from bs4 import BeautifulSoup
from collections import OrderedDict

import sys

# reload(sys)
# sys.setdefaultencoding('utf8')


# matching ensemble_id / ncbi_id / gene symbol,obtain family
def mapping_gene_id_symbol_by_dbs(dir_path):

    protein_ppi_vector_values_fname = dir_path + 'human_protein_ppi_vector_values_family.csv'
    vector_values = {}
    if not os.path.exists(protein_ppi_vector_values_fname):

        # 匹配 node_id 与 ensemble_id
        node_ensemble_id_mapping_pf = dir_path + 'protein_ensemble_id_for_node_mapping.csv'
        node_id_ppi_vector_pf = dir_path + 'snap-master/examples/node2vec/emb/human_protein_node_vector.csv'

        node_ensemble_id_mapping = {}
        ensemble_id_vector = {}
        with codecs.open(node_ensemble_id_mapping_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split(' ')
                if len(temp) != 2:
                    continue
                node_ensemble_id_mapping[temp[1].strip()] = temp[0].strip()
        print('node_ensemble_id_mapping: ' + str(len(node_ensemble_id_mapping.keys())))

        with codecs.open(node_id_ppi_vector_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split(' ')
                if len(temp) != 129:
                    continue
                ensemble_id_vector[node_ensemble_id_mapping[temp[0].strip()]] = ','.join(temp[1:])
        print('ensemble_id_vector: ' + str(len(ensemble_id_vector.keys())))
        ensemble_needed = set(node_ensemble_id_mapping.values())

        # 匹配 ensemble_id 与 ncbi_id
        ensemble_ncbi_id_mapping = {}
        ncbi_ids_for_ensemble = set()
        ensemble_ids = set()
        ncbi_ensemble_id_mapping_pf = dir_path + 'raw_data/gene2ensembl'
        with codecs.open(ncbi_ensemble_id_mapping_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('	')
                if len(temp) < 7:
                    continue
                ensemble_id = temp[-1].strip().split('.')[0].strip()
                if ensemble_id not in ensemble_needed:
                    continue
                if ensemble_id not in ensemble_ids:
                    ensemble_ncbi_id_mapping[ensemble_id] = set()
                    ensemble_ids.add(ensemble_id)
                ensemble_ncbi_id_mapping[ensemble_id].add(temp[1].strip())
                ncbi_ids_for_ensemble.add(temp[1].strip())
        print('ensemble_ncbi_id_mapping: ' + str(len(ensemble_ncbi_id_mapping.keys())))

        # ncbi_id 与 Family 匹配
        ncbi_id_family_id_mapping_pf = dir_path + 'raw_data/pfamA_ncbi_uniprot.txt'
        ncbi_id_family_ids_mapping = {}
        ncbi_ids_for_family = set()
        family_ids = set()
        with codecs.open(ncbi_id_family_id_mapping_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('	')
                if len(temp) != 3:
                    continue
                if temp[2].strip() not in ncbi_ids_for_ensemble:
                    continue
                if temp[2].strip() not in ncbi_ids_for_family:
                    ncbi_id_family_ids_mapping[temp[2].strip()] = set()
                    ncbi_ids_for_family.add(temp[2].strip())
                ncbi_id_family_ids_mapping[temp[2].strip()].add(temp[0].strip())
                family_ids.add(temp[0].strip())

        print('ncbi_id_family_ids_mapping: ' + str(len(ncbi_id_family_ids_mapping.keys())))
        print('ncbi_ids_for_family: ' + str(len(ncbi_ids_for_family)))
        print('family_ids: ' + str(len(family_ids)))

        # 匹配 family_id 与 protein vector
        protein_vector_family_ids = {}
        for key in ensemble_id_vector.keys():
            if key not in ensemble_ids:
                continue
            temp_family_ids = set()
            for ncbi_id in ensemble_ncbi_id_mapping[key]:
                if ncbi_id in ncbi_ids_for_family:
                    temp_family_ids = temp_family_ids.union(ncbi_id_family_ids_mapping[ncbi_id])
            if len(temp_family_ids) == 0:
                continue
            protein_vector_family_ids[ensemble_id_vector[key]] = temp_family_ids
        print('protein_vector_family_ids: ' + str(len(protein_vector_family_ids.keys())))

        # 匹配 family 相关的值
        # family_types = ['Family', 'Motif', 'Domain', 'Coiled-coil', 'Repeat', 'Disordered']
        family_id_values_pf = dir_path + 'raw_data/browse_pfam_family_details.csv'
        family_id_values = {}
        with codecs.open(family_id_values_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split(' 	')
                if len(temp) < 10:
                    continue
                if temp[1].strip() not in family_ids:
                    continue
                # family_type = [str(family_types.index(temp[2].strip()))]
                # family_id_values[temp[1].strip()] = family_type + temp[3:8]
                family_id_values[temp[1].strip()] = temp[3:8]
        print('family_id_values: ' + str(len(family_id_values.keys())))  # total: 17929，17571

        with open(protein_ppi_vector_values_fname, 'w') as output_file:
            for protein_vector in protein_vector_family_ids.keys():
                term_value = np.zeros(5, dtype=np.float32) #5,6
                for family_id in protein_vector_family_ids[protein_vector]:
                    value = family_id_values[family_id]
                    term_value += np.array(map(float, value), dtype=np.float32)
                count = float(len(protein_vector_family_ids[protein_vector]))
                avg_term_value =  [1.0 / count] + list(term_value / count)
                value_str = ','.join(map(str, avg_term_value))
                vector_values[protein_vector] = avg_term_value
                output_file.write(protein_vector + '$' + value_str + '\n')
    else:
        with codecs.open(protein_ppi_vector_values_fname, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if len(temp) < 2:
                    continue
                vector_values[temp[0].strip()] = map(float, temp[1].strip().split(','))
    return vector_values

# ppi of some scores
def prepare_ppi_scores_for_visualize(dir_path):

    protein_ppi_vector_values_fname = dir_path + 'human_protein_ppi_vector_values_scores.csv'
    vector_values = OrderedDict()
    if not os.path.exists(protein_ppi_vector_values_fname):

        # 匹配 node_id 与 ensemble_id
        node_ensemble_id_mapping_pf = dir_path + 'protein_ensemble_id_for_node_mapping.csv'
        node_id_ppi_vector_pf = dir_path + 'snap-master/examples/node2vec/emb/human_protein_node_vector.csv'

        node_ensemble_id_mapping = {}
        ensemble_id_vector = {}
        with codecs.open(node_ensemble_id_mapping_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split(' ')
                if len(temp) != 2:
                    continue
                node_ensemble_id_mapping[temp[1].strip()] = temp[0].strip()
        print('node_ensemble_id_mapping: ' + str(len(node_ensemble_id_mapping.keys())))

        with codecs.open(node_id_ppi_vector_pf, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split(' ')
                if len(temp) != 129:
                    continue
                ensemble_id_vector[node_ensemble_id_mapping[temp[0].strip()]] = ','.join(temp[1:])
        print('ensemble_id_vector: ' + str(len(ensemble_id_vector.keys())))
        ensemble_needed = set(ensemble_id_vector.keys())

        # 匹配出相关的scores
        subdir_path = '/species_protein_pair_score_update/'
        dir_file_names = [f for f in listdir(subdir_path) if f.endswith('csv')]
        ensemble_id_scores = {}
        ensemble_ids_for_scores = set()
        for file_name in dir_file_names:
            with codecs.open(subdir_path + file_name, "rb", "utf-8") as input_file:
                for line in input_file:
                    if line.strip()[0] == ',':
                        continue
                    temp = line.strip().split(',')
                    protein1 = temp[1].strip().split('.')[1].strip()
                    protein2 = temp[2].strip().split('.')[1].strip()
                    if protein1 in ensemble_needed:
                        if protein1 not in ensemble_ids_for_scores:
                            ensemble_id_scores[protein1] = []
                            ensemble_ids_for_scores.add(protein1)
                        ensemble_id_scores[protein1].append(temp[3:])
                    if protein2 in ensemble_needed:
                        if protein2 not in ensemble_ids_for_scores:
                            ensemble_id_scores[protein2] = []
                            ensemble_ids_for_scores.add(protein2)
                        ensemble_id_scores[protein2].append(temp[3:])
        print('ensemble_id_scores: ' + str(len(ensemble_id_scores.keys())))

        with open(protein_ppi_vector_values_fname, 'w') as output_file:
            for ensemble_id in ensemble_id_scores.keys():
                term_scores = np.zeros(8, dtype=np.float32)
                for scores in ensemble_id_scores[ensemble_id]:
                    term_scores += np.array(map(float, scores), dtype=np.float32)
                sum_scores_list = float(len(ensemble_id_scores[ensemble_id]))
                term_scores_avg = list(term_scores / sum_scores_list)
                value_str = ','.join(map(str, term_scores_avg))
                protein_vector = ensemble_id_vector[ensemble_id]
                vector_values[protein_vector] = term_scores_avg
                output_file.write(protein_vector + '$' + value_str + '\n')
    else:
        with codecs.open(protein_ppi_vector_values_fname, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if len(temp) < 2:
                    continue
                vector_values[temp[0].strip()] = [float(x) for x in temp[1].strip().split(',')]

    print('vector_values: ' + str(len(vector_values.keys())))
    return vector_values


if __name__ == '__main__':

    dir_path = ''
    # mapping_gene_id_symbol_by_dbs(dir_path)
    prepare_ppi_scores_for_visualize(dir_path)

