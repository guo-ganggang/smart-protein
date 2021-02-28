#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 17/11/2018 9:13 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : data_helper.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import codecs
from itertools import islice
import numpy as np
from os import listdir
import re
import os
import nltk
import nltk.data
from nltk.corpus import stopwords
from collections import Counter
from collections import OrderedDict

# import sys
# reload(sys)
# sys.setdefaultencoding('utf8')


# Classify the MeSH terms
def compute_mesh_terms_classification_score(dir_path):

    raw_data_mesh_term_fp = dir_path + 'MeSH_terms/parse_MeSH_terms_info.csv'
    pubmed_mesh_term_fp = dir_path + 'MeSH_terms/pubmed_id_mesh_words.csv'

    # 读入标准文件
    mesh_term_tree = {}
    tree_classifications = OrderedDict()
    biggest_class = set()
    with codecs.open(raw_data_mesh_term_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) != 4:
                continue
            if temp[-1].strip() == '':
                continue
            term_trees = temp[-1].strip().split(';')
            while '' in term_trees:
                term_trees.remove('')
            if len(term_trees) == 0:
                continue
            for tree in term_trees:
                classes = tree.split('.')
                if len(classes) < 2:
                    one_class = classes[0].strip()
                    if one_class not in biggest_class:
                        tree_classifications[one_class] = []
                        biggest_class.add(one_class)
                else:
                    one_class = classes[0].strip()
                    two_class = classes[1].strip()
                    if one_class not in biggest_class:
                        tree_classifications[one_class] = []
                        biggest_class.add(one_class)
                    tree_classifications[one_class].append(two_class)
            mesh_term_tree[temp[0].strip()] = term_trees

    target_tree_classes = OrderedDict()
    for key in tree_classifications.keys():
        if len(tree_classifications[key]) > 500:
            print(key,str(len(tree_classifications[key])))
            target_tree_classes[key] = tree_classifications[key]
    print('mesh_term_tree %d' % len(mesh_term_tree.keys()))
    print('target_tree_classes %d' % len(target_tree_classes))

    # 读入语料文献关键词
    pubmed_mesh_terms = set()
    tree_mesh_terms = set(mesh_term_tree.keys())
    with codecs.open(pubmed_mesh_term_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if temp[1] == '':
                continue
            temp_terms = set(temp[1].strip().split(';;'))
            for tt in temp_terms:
                if tt.strip() == '':
                    continue
                term = re.sub('\*','',tt.strip().split('/')[0].strip())
                if term in tree_mesh_terms:
                    pubmed_mesh_terms.add(term)
    print('pubmed_mesh_terms %d' % len(pubmed_mesh_terms)) #23287

    # 计算每个MeSH term 的分值

    pubmed_mesh_term_score = {}
    target_tree_classes_list = list(target_tree_classes.keys())
    for term in pubmed_mesh_terms:
        tree_classes = mesh_term_tree[term]
        temp_ids = []
        # print(tree_classes)
        for tree_class in tree_classes:
            level = tree_class.split('.')[0].strip()
            temp_ids.append(level)
        if len(temp_ids) == 0:
            continue
        class_index = max(temp_ids, key=temp_ids.count)
        if class_index not in target_tree_classes_list:
            continue
        score = target_tree_classes_list.index(class_index)+1
        # pubmed_mesh_term_score[term]=float(score) / len(target_tree_classes_list)
        pubmed_mesh_term_score[term] = float(score)
    print('pubmed_mesh_term_score %d' % len(pubmed_mesh_term_score.keys()))

    return pubmed_mesh_term_score

# Clean up the MeSH terms
def clean_mesh_words(dir_path):

    clean_pubmed_mesh_term_score = {}

    # Load the punkt tokenizer
    download_nltk_corpus_path = '/nltk_data/'
    if not os.path.exists(download_nltk_corpus_path):
        os.makedirs(download_nltk_corpus_path)
        nltk.download() # 下载到上面那个路径下
    english_stopwords = stopwords.words('english') # 自动检索上面那个路径
    english_punctuations = [',', '.', ':', ';', '?', '(', ')', '[', ']', '&', '!', '*', '@', '#', '$', '%',
                            '>','<','\'s','\'','-/-','``',"''",'--','-']
    pubmed_mesh_term_score = compute_mesh_terms_classification_score(dir_path)
    for key in pubmed_mesh_term_score.keys():

        corpus = key.strip()
        texts_tokenized = [nltk.word_tokenize(sentence) for sentence in nltk.sent_tokenize(corpus)]
        texts_filtered_stopwords = [[word for word in document if not word in english_stopwords] for document in \
                                    texts_tokenized]
        texts_filtered = [[word for word in document if not word in english_punctuations] for document in \
                          texts_filtered_stopwords]
        # texts_stemmed = [[st.stem(word) for word in docment] for docment in texts_filtered]
        for text in texts_filtered:
            if len(text) == 0:
                continue
            words = []
            for word in text:
                if word.strip() == '':
                    continue
                if word.strip().istitle():
                    words.append(word.strip().lower())
                else:
                    words.append(word.strip())
            sent = ' '.join(words)
            if ''.join(words) == '':
                continue
            clean_pubmed_mesh_term_score[sent] = pubmed_mesh_term_score[key]
    print('clean_pubmed_mesh_term_score %d' % len(clean_pubmed_mesh_term_score.keys()))

    return clean_pubmed_mesh_term_score

# Read in the training word vector
def get_word_vectors(corpus_vector):

    term_vectors = {}
    with open(corpus_vector) as infile:
        for line in infile:
            line_parts = line.rstrip().split()
            # skip first line with metadata in word2vec text file format
            if len(line_parts) > 2:
                term, vector_values = line_parts[0], line_parts[1:]
                term_vectors[term] = np.array(map(float, vector_values), dtype=np.float32)
    return term_vectors

# Prepare the Mesh term vector
def get_term_vectors(dir_path):
    corpus_vector = "/w2v_model/human_protein_literature_corpus_sg1_s128_w5_m1.vector"
    vector_score = {}

    write_out_fp = dir_path + 'human_protein_mesh_vector_values_concept.csv'

    if not os.path.exists(write_out_fp):
        with open(write_out_fp, 'w') as output_file:
            clean_pubmed_mesh_term_score = clean_mesh_words(dir_path)
            term_vectors = get_word_vectors(corpus_vector)
            model_words = set(term_vectors.keys())

            for key in clean_pubmed_mesh_term_score.keys():
                words = key.split(' ')
                init_vector = np.zeros((128))
                count_words = 0
                for word in words:
                    if word in model_words:
                        init_vector += term_vectors[word]
                        count_words += 1
                if count_words == 0:
                    continue
                vector = map(str, init_vector / float(count_words))
                vector_score[','.join(vector)] = clean_pubmed_mesh_term_score[key]
                output_file.write(','.join(vector) + '$' + str(clean_pubmed_mesh_term_score[key]) + '\n')
    else:
        with codecs.open(write_out_fp, 'rb', 'utf-8') as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if len(temp) < 2:
                    continue
                vector_score[temp[0].strip()] = temp[1].strip()
    print('vector_score %d' % len(vector_score.keys()))

    return vector_score




if __name__ == '__main__':

    dir_path = ''
    get_term_vectors(dir_path)
