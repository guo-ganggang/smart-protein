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

# import sys
#
# reload(sys)
# sys.setdefaultencoding('utf8')


def get_entity_vectors(dir_path):

    text_corpus_vector_fp = dir_path + 'w2v_model/protein_go_term_corpus_sg1_s128_w5_m1.vector'

    go_id_term_words = {}
    train_go_id_term_fname = dir_path + 'text_corpus/train_go_id_term_clean.csv'
    with codecs.open(train_go_id_term_fname, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            # 名称与别名分别存入列表中
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
    return total_go_id_vec

def go_function_vec_classification(dir_path):

    vector_class = {}
    write_out_fp = dir_path + 'human_protein_go_vector_function_class.csv'
    if not os.path.exists(write_out_fp):
        total_go_id_vec = get_entity_vectors(dir_path)
        go_ids = set(total_go_id_vec.keys())
        go_classifications = ['Component', 'Process', 'Function']
        with open(write_out_fp, 'w') as output_file:
            ncbi_go_term_fname = dir_path + 'raw_data/gene2go'
            with codecs.open(ncbi_go_term_fname, "rb", "utf-8") as input_file:
                for line in input_file:
                    if line.strip()[0] == '#':
                        continue
                    temp = line.strip().split('	')
                    if len(temp) < 8:
                        continue
                    if (temp[2].strip() == '-') or (temp[-1].strip() == '-'):
                        continue
                    if temp[2].strip() in go_ids:
                        vector_str = ','.join([str(x) for x in total_go_id_vec[temp[2].strip()]])
                        vector_class[vector_str] = go_classifications.index(temp[-1].strip())
                        output_file.write(vector_str + '$' + str(vector_class[vector_str]) + '\n')

    else:
        with codecs.open(write_out_fp, 'rb', 'utf-8') as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if len(temp) < 2:
                    continue
                vector_class[temp[0].strip()] = temp[1].strip()
    print('vector_class: ' + str(len(vector_class.keys())))

    return vector_class

if __name__ == '__main__':

    dir_path = ''
    go_function_vec_classification(dir_path)
