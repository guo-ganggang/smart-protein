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
import data_helper as dh

# Obtain a vector representation of the protein GO
def get_entity_vectors(dir_path):

    if not os.path.exists(dir_path + 'w2v_model/'):
        os.makedirs(dir_path + 'w2v_model/')

    text_corpus_vector_fp = dir_path + 'w2v_model/protein_go_term_corpus_sg1_s128_w5_m1.vector'
    if not os.path.isfile(text_corpus_vector_fp):
        dh.train_w2v_model(dir_path)

    go_id_term_words = {}
    train_go_id_term_fname = dir_path + 'text_corpus/train_go_id_term_clean.csv'
    if not os.path.isfile(train_go_id_term_fname):
        dh.clean_go_term_corpus(dir_path)
    with codecs.open(train_go_id_term_fname, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            # Names and aliases are stored separately in the list
            go_id_term_words[temp[0].strip()] = temp[1].strip().split(' ')
    print('go_id_term_words: ' + str(len(go_id_term_words.keys())))

    word_vectors = {}
    with codecs.open(text_corpus_vector_fp, "rb", "utf-8") as infile:
        for line in infile:
            line_parts = line.rstrip().split()
            # skip first line with metadata in word2vec text file format
            if len(line_parts) > 2:
                entity, vector_values = line_parts[0], line_parts[1:]
                word_vectors[entity] = np.array(map(float, vector_values), dtype=np.float32)
    print('word_vectors: ' + str(len(word_vectors.keys())))
    total_corpus_words = set(word_vectors.keys())

    total_go_id_vec = {}
    for key in go_id_term_words.keys():
        term_vec = np.zeros(128, dtype=np.float32)
        flag = 0
        for word in go_id_term_words[key]:
            if word in total_corpus_words:
                word_vector = word_vectors[word]
                term_vec += word_vector
                flag = 1
        if flag == 0:
            continue
        total_go_id_vec[key] = term_vec
    print('total_go_id_vec: ' + str(len(total_go_id_vec.keys())))
    total_go_ids = set(total_go_id_vec.keys())

    write_out_fp = dir_path + 'train_protein_symbol_go_ids_vector.csv'
    train_protein_symbol_go_ids_fp = dir_path + 'text_corpus/train_protein_symbol_go_ids.csv'
    with codecs.open(write_out_fp, 'w') as output_file:
        with codecs.open(train_protein_symbol_go_ids_fp, "rb", "utf-8") as infile:
            for line in infile:
                temp = line.strip().split('$')
                temp_ids = temp[1].strip().split(';')
                term_vec = np.zeros(128, dtype=np.float32)
                temp_count = 0
                for go_id in temp_ids:
                    if go_id.strip() == '':
                        continue
                    if go_id.strip() in total_go_ids:
                        go_id_vector = total_go_id_vec[go_id.strip()]
                        term_vec += go_id_vector
                        temp_count += 1
                if temp_count == 0:
                    continue
                avg_term_vec = term_vec / float(temp_count)
                values = ','.join(map(str, avg_term_vec))
                output_file.write(temp[0].strip() + '$' + values + '\n')


if __name__ == '__main__':

    dir_path = ''
    get_entity_vectors(dir_path)