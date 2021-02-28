#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 28/12/2017 下午 8:43
# @Author  : GUO Ganggang
# @email   : ganggangguo@csu.edu.cn
# @Site    : 
# @File    : protein_literature_2vec.py
# @Software: PyCharm

###############
###Intermediate process code, for user reference only
###############


import warnings
warnings.filterwarnings(action='ignore', category=UserWarning, module='gensim')

import nltk
import nltk.data
from nltk.corpus import stopwords
import os
import codecs
import logging
import os.path
import sys
import multiprocessing
from gensim.models import Word2Vec
from gensim.models.word2vec import Text8Corpus
import numpy as np
from nltk import data
import data_helper as dh

# Generates a vector representation of the protein Mesh
def get_entity_vectors(dir_path):

    text_corpus_vector_fp = dir_path + 'w2v_model/human_protein_literature_corpus_sg1_s128_w5_m1.vector'

    protein_symbol_mesh_words = {}
    mesh_words_all = set()
    protein_entry_pmid_title_abstract_fp = dir_path + 'protein_official_symbol_pmids_mesh_words_for_train_modify_mesh_terms_clean.csv'
    if not os.path.isfile(protein_entry_pmid_title_abstract_fp):
        dh.clean_mesh_words(dir_path)
    with codecs.open(protein_entry_pmid_title_abstract_fp, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            words = temp[1].strip().split(' ')
            for word in words:
                mesh_words_all.add(word.strip())
            protein_symbol_mesh_words[temp[0].strip()] = words
    print('protein_symbol_mesh_words: ' + str(len(protein_symbol_mesh_words.keys()))) # v1 6866,v2
    print('mesh_words_all: ' + str(len(mesh_words_all)))

    entity_vectors = {}
    vector_words = set()
    with codecs.open(text_corpus_vector_fp, "rb", "utf-8") as infile:
        for line in infile:
            line_parts = line.rstrip().split()
            # skip first line with metadata in word2vec text file format
            if len(line_parts) > 2:
                entity, vector_values = line_parts[0], line_parts[1:]
                if entity in mesh_words_all:
                    vector_words.add(entity)
                    entity_vectors[entity] = np.array(map(float, vector_values), dtype=np.float32)
    print('vector_words: ' + str(len(vector_words)))

    protein_entity_vector_fp = dir_path + 'train_protein_sympol_mesh_vector.csv'
    if not os.path.isfile(protein_entity_vector_fp):
        with open(protein_entity_vector_fp, 'a') as output_file:
            for key in protein_symbol_mesh_words.keys():
                temp_vector = np.zeros((128))
                for word in protein_symbol_mesh_words[key]:
                    if word in vector_words:
                        temp_vector += entity_vectors[word]
                values = ','.join([str(x) for x in (temp_vector/float(len(protein_symbol_mesh_words[key])))])
                output_file.write(key + '$' + values + '\n')


if __name__ == '__main__':

    dir_path = ''
    get_entity_vectors(dir_path)