#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 4/11/2018 3:59 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : seq_bio2vec_visualize.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import tsne_3gram as t3
import ngrams_properties as pro
import pickle
import os
from gensim.models import word2vec
import numpy as np
import data_helper

def protein_space(dir_path,norm='no-norm'):

    print("Loading 3gram vector")
    model_bin = dir_path+"bio2vec_model/protein_sequence_corpus_3letter_sg1_s128_w25_m1.bin"
    model = word2vec.Word2Vec.load(model_bin)
    # filter_model = data_helper.get_ngram_vectors(model_3gram)
    print("... Ok\n")

    print("Making tsne")
    tsne = t3.BioTsne()
    labels = model.wv.vocab.keys()
    # print(type(labels))
    # labels = list(filter_model.keys())

    file_path_pkl = dir_path + "/visualize/ngram_2D_vector.pkl"
    tsne.make_tsne(file_path_pkl,model)
    f = open(file_path_pkl, "rb")
    vectors = pickle.load(f)
    print(vectors.shape)

    property_list = pro.make_property_list(labels)
    print(len(property_list))
    final_embedding = tsne.link_with_vector(vectors, property_list)

    if norm != 'norm':
        X_tsne = final_embedding
    else:
        x_min, x_max = final_embedding.min(0), final_embedding.max(0)
        X_tsne = (final_embedding - x_min) / (x_max - x_min)  # 归一化
    print("... OK\n")

    color_sets = ['RdYlBu','CMRmap','RdYlGn','afmhot','autumn','spring','PiYG','YlGn','YlGnBu','Set1','Set2','Set3']
    print("Visualization:Normalized distributions of biochemical and biophysical properties in protein-space")

    color_set = color_sets[-3] #3,-3
    fig_path = dir_path+'/visualize/Distributions of amino acid biochemical properties in protein-sequence-space.png'
    tsne.visualization(fig_path,X_tsne,color_sets=color_set)


if __name__ == '__main__':

    dir_path = ''
    protein_space(dir_path)