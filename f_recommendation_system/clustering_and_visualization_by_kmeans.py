#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 14/4/2019 11:09 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : clustering_and_visualization_by_kmeans.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import codecs
from itertools import islice
import time


# Read in the compressed gene characteristic data set and get the selected gene name set
def prepare_dataset_for_predicting(dir_path, protein_symbols_selected):

    # 读入全部蛋白质质及其相应特征向量
    dataset_for_clustering = []
    raw_data_fp = dir_path + 'protein_symbol_2D_features_with_dimensionality_reduction.csv'
    with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if len(temp) < 3:
                continue
            if temp[0].strip() not in protein_symbols_selected:
                continue
            feature_values = [float(x) for x in temp[1:]]
            dataset_for_clustering.append(feature_values)
    dataset_for_clustering = np.array(dataset_for_clustering)
    print('dataset_for_clustering: ' + str(dataset_for_clustering.shape))

    return dataset_for_clustering

# Clustering and visualization
def dimensionality_reduction_clustering_visualize(dir_path,reduced_data,n_clusters,gene_names):

    # k-means++ 聚类
    kmeans = KMeans(init='k-means++', n_clusters=n_clusters)
    kmeans.fit(reduced_data)

    # Step size of the mesh. Decrease to increase the quality of the VQ.
    h = .02  # point in the mesh [x_min, x_max]x[y_min, y_max].

    # Plot the decision boundary. For that, we will assign a color to each
    x_min, x_max = reduced_data[:, 0].min() - 1, reduced_data[:, 0].max() + 1
    y_min, y_max = reduced_data[:, 1].min() - 1, reduced_data[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))# 将两个一维数组变为二维矩阵

    # Obtain labels for each point in mesh. Use last trained model.
    Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    plt.figure(1)
    plt.clf()
    plt.imshow(Z, interpolation='nearest',
               extent=(xx.min(), xx.max(), yy.min(), yy.max()),
               cmap=plt.cm.Paired,
               aspect='auto', origin='lower')

    plt.plot(reduced_data[:, 0], reduced_data[:, 1], 'g*', markersize=10)
    for i in range(len(gene_names)):
        print(reduced_data[:, 0][i], reduced_data[:, 1][i])
        plt.text(reduced_data[:, 0][i], reduced_data[:, 1][i], gene_names[i], fontsize=10, color="k", style="italic", weight="light", verticalalignment='center',
              horizontalalignment='right')
    # Plot the centroids as a white X
    centroids = kmeans.cluster_centers_
    plt.scatter(centroids[:, 0], centroids[:, 1],
                marker='x', s=100, linewidths=3,
                color='w', zorder=10)
    plt.title('k-means++ clustering and centroids are marked by white cross')
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xticks(())
    plt.yticks(())

    # 保存可视化结果，并且show
    now_time = str(int(time.time()))
    save_fig_path = dir_path + now_time + '_' + str(n_clusters) + 'clusters_visualiaztion_by_kmeans.png'
    plt.savefig(save_fig_path)
    plt.show()


if __name__ == '__main__':

    dir_path = ''
    n_clusters = 5
    protein_symbols_selected = ['RNF14','NDP','GRIN1','RAB40B','ATP2A3','SDK1','ERICH2','FRG1','KHDRBS2','SMOC2',
                                'PPIB', 'CD177', 'GOLGA8J', 'DNAJC5', 'CHD8', 'RAB3GAP2', 'TTC8', 'CSNK2A2', 'FAM49B',
                                'PKD1','FABP9', 'CILP', 'VPS13A', 'TMEM164', 'ACTG1', 'DAZ2', 'AKAP7', 'SOX9', 'PCLO',
                                'DNAJC2','CDK17','LRRC4C','SYNE2','UBQLN1','GABARAP','ZNRF4','PSMD1','ROPN1L','MCTS1']
    dataset_for_clustering = prepare_dataset_for_predicting(dir_path, protein_symbols_selected)
    dimensionality_reduction_clustering_visualize(dir_path, dataset_for_clustering, n_clusters,protein_symbols_selected)

