#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 4/11/2018 3:59 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : tsne_visual.py
# @Software: Mining from a specific set of proteins in human sperm

import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pickle
import os

class BioTsne:
    def __init__(self):
        print('TSNE is running..')

    # making tsne
    def make_tsne(self, file_path , model):

        if not os.path.isfile(file_path):
            # make tsne
            X = model[model.wv.vocab]
            tsne = TSNE(n_components=2)
            X_tsne = tsne.fit_transform(X)
            print(X_tsne.shape)

            # save X_tsne
            f = open(file_path,"wb")
            pickle.dump(X_tsne , f)

            f.close()

    def link_with_vector(self, vectors, property_list):
        # print np.append(vectors , property_list,axis=1)
        append_vec_pro_np = np.append(vectors , property_list,axis=1)
        print(append_vec_pro_np.shape)
        return append_vec_pro_np

    def visualization(self, file_path,X_tsne,color_sets='YlGnBu'):
        # load final_embedding data
        plt.style.use('bmh')
        fig , axarr = plt.subplots(3,3, figsize=(15,15))
        #set marker size
        marker_size=1

        cm = plt.cm.get_cmap(color_sets)
        # set scatter
        g1 = axarr[0,0].scatter(X_tsne[:,0], X_tsne[:, 1], marker_size ,X_tsne[:,2],cmap=cm)
        axarr[0,0].set_title("Hydrophobicity",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g1, ax=axarr[0,0])

        g2 = axarr[0,1].scatter(X_tsne[:,0], X_tsne[:, 1], marker_size ,X_tsne[:,3],cmap=cm)
        axarr[0,1].set_title("Hydrophicility",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g2, ax=axarr[0,1])

        g3 = axarr[0,2].scatter(X_tsne[:,0], X_tsne[:, 1], marker_size ,X_tsne[:,4],cmap=cm)
        axarr[0,2].set_title("Hydrogen bond",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g3, ax=axarr[0,2])

        g4 = axarr[1,0].scatter(X_tsne[:,0], X_tsne[:, 1], marker_size ,X_tsne[:,5],cmap=cm)
        axarr[1,0].set_title("Volumes of side chains",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g4, ax=axarr[1,0])

        g5 = axarr[1,1].scatter(X_tsne[:,0], X_tsne[:, 1], marker_size ,X_tsne[:,6],cmap=cm)
        axarr[1,1].set_title("Polarity",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g5, ax=axarr[1,1])

        g6 = axarr[1,2].scatter(X_tsne[:,0], X_tsne[:, 1], marker_size ,X_tsne[:,7],cmap=cm)
        axarr[1,2].set_title("Polarizability",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g6, ax=axarr[1,2])

        g7 = axarr[2, 0].scatter(X_tsne[:, 0], X_tsne[:, 1], marker_size, X_tsne[:, 8],cmap=cm)
        axarr[2, 0].set_title("SASA",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g7, ax=axarr[2,0])

        g8 = axarr[2, 1].scatter(X_tsne[:, 0], X_tsne[:, 1], marker_size, X_tsne[:, 9],cmap=cm)
        axarr[2, 1].set_title("NCI",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g8, ax=axarr[2,1])

        g9 = axarr[2, 2].scatter(X_tsne[:, 0], X_tsne[:, 1], marker_size, X_tsne[:, 10],cmap=cm)
        axarr[2, 2].set_title("MASS",fontsize=14).set_fontname('Times New Roman')
        fig.colorbar(g9, ax=axarr[2,2])

        plt.savefig(file_path)

        plt.show()
