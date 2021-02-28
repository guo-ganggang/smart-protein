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
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
import warnings


def generate_cmap(colors):
    """自定义颜色"""
    values = range(len(colors))
    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)

def protein_space(dir_path,algorithm = 'tsne'):

    print("Loading term vector")
    vector_class = dh.go_function_vec_classification(dir_path)

    print("Making tsne")
    temp_vectors = []
    property_list = []
    for vector in vector_class.keys():
        temp_vectors.append([float(x) for x in vector.split(',')])
        property_list.append(vector_class[vector])
    print('property_list: ' + str(len(property_list)))

    file_path_pkl = dir_path + 'visualize/go_term_3D_vector.pkl'
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

    # color_list = map(set_color, property_list)
    cm = generate_cmap(['#FF6666', '#2E8B57', '#0066CC'])
    marker_size = 1
    # color_sets = ['YlOrRd','RdYlBu','CMRmap','RdYlGn','PRGn_r','autumn_r','spring','PiYG','YlGn','YlGnBu','Set1','Set2']
    # cm = plt.cm.get_cmap(color_sets[4])
    print("Visualization:Normalized distributions of biochemical and biophysical properties in protein-space")
    fig_path = dir_path + 'visualize/Distribution of GO terms in literature-space.pdf'
    fig = plt.figure(figsize=(13, 10))

    # 2D 可视化
    # sc = plt.scatter(vectors[:,0], vectors[:,1],marker_size, c=property_list, cmap=cm)
    # fig.colorbar(sc)
    # plt.savefig(fig_path)
    # plt.show()

    # 3D 可视化
    ax = Axes3D(fig)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()

    sc = ax.scatter(vectors[:, 0], vectors[:, 1], vectors[:, 2], alpha=0.4, s=marker_size, c=property_list, cmap=cm)
    # ax.set_title("Distribution of GO terms in literature-space")
    ax.set_xlabel('X',fontsize=10).set_fontname('Times New Roman')
    ax.set_ylabel('Y',fontsize=10).set_fontname('Times New Roman')
    ax.set_zlabel('Z',fontsize=10).set_fontname('Times New Roman')
    ax.view_init(elev=10, azim=-35)
    fig.colorbar(sc, shrink=0.5, aspect=15)
    plt.savefig(fig_path)  # ,bbox_inches='tight'
    plt.show()




if __name__ == '__main__':

    dir_path = ''
    protein_space(dir_path)