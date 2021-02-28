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

from sklearn.manifold import TSNE
import pickle
import os
import numpy as np
import data_helper as dh
import matplotlib.pyplot as plt
# from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import warnings
from matplotlib.colors import LinearSegmentedColormap

def generate_cmap(colors):
    """custom colors"""
    values = range(len(colors))
    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


def protein_space(dir_path):

    print("Loading term vector")
    vector_score = dh.get_term_vectors(dir_path)

    print("Making tsne")
    temp_vectors = []
    property_list = []
    for vector in vector_score.keys():
        temp_vectors.append([float(x) for x in vector.split(',')])
        property_list.append(vector_score[vector])

    file_path_pkl = dir_path + 'visualize/mesh_term_3D_vector.pkl'
    if not os.path.exists(file_path_pkl):
        X = np.array(temp_vectors)
        tsne = TSNE(n_components=3)
        X_tsne = tsne.fit_transform(X)

        # save X_tsne
        f = open(file_path_pkl, "wb")
        pickle.dump(X_tsne, f)
        f.close()
    print("X_tsne... Ok\n")

    # load X_tsne
    f = open(file_path_pkl, "rb")
    vectors = pickle.load(f)
    print(vectors.shape)

    # final_embedding = np.append(vectors, property_list, axis=1)
    # print(final_embedding.shape)

    # color_sets = ['YlOrRd','RdYlBu','CMRmap','RdYlGn','PRGn_r','autumn_r','spring','PiYG','YlGn','YlGnBu','Set1','Set2']
    # color_set = color_sets[-3]
    color_set = generate_cmap(['#339999', '#CCCCFF', '#663366'])  ##FFFF00,#CCCCFF,#339999,#CCFF00,#663366
    cm = plt.cm.get_cmap(color_set)
    marker_size = 1
    # plt.style.use('bmh')

    print('property_list: '+str(len(property_list)))

    fig_path = dir_path + 'visualize/Distribution of MeSH terms in literature-space.png' #color_set,'custom'

    # # 2D 图
    # fig = plt.figure(figsize=(13, 7))
    # sc = plt.scatter(vectors[:,0], vectors[:,1],marker_size, c=property_list, cmap=cm)
    # fig.colorbar(sc)
    # plt.savefig(fig_path)
    # plt.show()

    # 3D 图
    fig = plt.figure(figsize=(13,10))
    # ax = fig.add_subplot(111, projection='3d')
    ax = fig.gca(projection='3d')
    # ax = Axes3D(fig)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()

    sc = ax.scatter(vectors[:,0], vectors[:,1], vectors[:,2], alpha=0.4,s=marker_size,c=property_list, cmap=cm)
    # ax.set_title("Distribution of MeSH concept trees in literature-space")
    ax.set_xlabel('X',fontsize=10).set_fontname('Times New Roman')
    ax.set_ylabel('Y',fontsize=10).set_fontname('Times New Roman')
    ax.set_zlabel('Z',fontsize=10).set_fontname('Times New Roman')
    ax.view_init(elev=10, azim=-35)
    fig.colorbar(sc,shrink=0.5, aspect=15)
    plt.savefig(fig_path) # ,bbox_inches='tight'
    plt.show()


if __name__ == '__main__':

    dir_path = ''
    protein_space(dir_path)