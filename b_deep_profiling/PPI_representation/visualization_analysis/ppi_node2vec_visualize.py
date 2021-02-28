#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 4/11/2018 3:59 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : seq_bio2vec_visualize.py
# @Software: Mining from a specific set of proteins in human sperm

from sklearn.manifold import TSNE
import pickle
import os
import numpy as np
import data_helper as dh
import tsne_visual as tf
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn import preprocessing

def generate_cmap(colors):
    """custom colors"""
    values = range(len(colors))
    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)

def protein_space(dir_path,norm = 'non-norm'):

    print("Loading term vector")
    # vector_values = dh.mapping_gene_id_symbol_by_dbs(dir_path)
    vector_values = dh.prepare_ppi_scores_for_visualize(dir_path)

    print("Making tsne")
    temp_vectors = []
    property_list = np.array(list(vector_values.values())) # raw

    # min_max_scaler = preprocessing.MinMaxScaler()
    # property_list = min_max_scaler.fit_transform(np.array(list(vector_values.values()))) # min_max

    print(property_list.shape)

    for vector in vector_values.keys():
        temp_vectors.append([float(x) for x in vector.split(',')])

    tsne = tf.BioTsne()

    file_path_pkl = dir_path + 'visualize/scores_term_2D_vector.pkl'
    # file_path_pkl = dir_path + 'visualize/family_term_2D_vector.pkl'
    if not os.path.exists(file_path_pkl):
        tsne.make_tsne(file_path_pkl,temp_vectors)
    print("X_tsne... Ok\n")

    # load X_tsne
    f = open(file_path_pkl, "rb")
    vectors = pickle.load(f)
    print(vectors.shape)

    final_embedding = tsne.link_with_vector(vectors, property_list)

    if norm != 'norm':
        X_tsne = final_embedding
    else:
        x_min, x_max = final_embedding.min(0), final_embedding.max(0)
        X_tsne = (final_embedding - x_min) / (x_max - x_min)  # 归一化
    print("... OK\n")

    # color_set = generate_cmap(['#99CC66', '#FF6600', '#663366']) ##FFFF00,#CCCCFF,#339999,#CCFF00
    color_sets = ['YlOrRd','RdYlBu','CMRmap','RdYlGn','PRGn_r','autumn_r','spring','PiYG','YlGn','YlGnBu','Set1','Set2']
    color_set = color_sets[-1]  # -2,-1
    fig_path = dir_path + '/Distribution of network edges weight proterties in protein-PPI-space.png'
    # fig_path = dir_path + '/scores_ppi_node2vec_tsne_' + 'custom' + '_' + norm + '.png'
    # tsne.visualization_for_family(fig_path, X_tsne, color_sets=color_set)
    tsne.visualization_for_scores(fig_path, X_tsne, color_sets=color_set)


if __name__ == '__main__':

    dir_path = ''
    protein_space(dir_path)