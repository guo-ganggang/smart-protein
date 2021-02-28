#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 6/11/2018 5:20 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : protein_semsim_2vec.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import networkx as nx
from node2vec import Node2Vec
import matplotlib.pyplot as plt
from os import listdir
import re
import codecs
from itertools import islice
import platform

### read edge list to networkx ###
# the format of each line: (src dst whole_duration/total_duration total_duration)
# stringId_A	stringId_B	preferredName_A	preferredName_B	ncbiTaxonId	score	nscore
# 	fscore	pscore	ascore	escore	dscore	tscore

def create_graph(dir_path,path_splitting):

    file_path = dir_path + 'genes_scores.csv'
    protein_name_index = {}
    count_index = 0
    G = nx.Graph() # networkx有四种图 Graph 、DiGraph、MultiGraph、MultiDiGraph，分别为无多重边无向图、无多重边有向图、有多重边无向图、有多重边有向图。
    with codecs.open(file_path, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) != 13:
                continue
            if (temp[5].strip() == '') or (float(temp[5].strip()) >= 0.8):
                continue
            for preferred_name in temp[2:4]:
                if preferred_name.strip() not in protein_name_index.keys():
                    protein_name_index[preferred_name.strip()] = count_index
                    count_index += 1
            n1 = int(protein_name_index[temp[2].strip()])
            n2 = int(protein_name_index[temp[3].strip()])
            weight = float(temp[5].strip())
            G.add_weighted_edges_from([(n1, n2, weight)]) #G.add_edges_from([(n1, n2)])
    number_nodes = len(protein_name_index.keys())
    print(number_nodes)

    pos = nx.spring_layout(G)
    # 按权重划分为重权值得边和轻权值的边
    # elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0.5]
    # esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= 0.5]
    nx.draw(G,pos,edge_color='k',style='solid',node_color = 'y',node_shape='o',node_size =5)
    # nx.draw_networkx_edges(G, pos, edgelist=elarge, edge_color='k', style='solid', node_color='y', node_shape='o',
    #                        node_size=5)
    # nx.draw_networkx_edges(G, pos, edgelist=esmall, alpha=0.5, edge_color='b', style='dashed', node_shape='o',
    #                        node_size=5)
    plt.savefig(dir_path + 'vector_model'+path_splitting+'gene_ppi_score.pdf')  # 输出方式1: 将图像存为一个png格式的图片文件
    # plt.show()  # 输出方式2: 在窗口中显示这幅图像

    return G

def read_train_save(dir_path,path_splitting='/'):
    sysstr = platform.system() # 判断操作系统
    if sysstr != "Windows":
        path_splitting = '/'
    # FILES
    EMBEDDING_FILENAME = dir_path + 'vector_model'+path_splitting+'ppi_score_embeddings.emb'
    EMBEDDING_MODEL_FILENAME = dir_path + 'vector_model'+path_splitting+'ppi_score_embeddings.model'

    # Create a graph
    # graph = nx.fast_gnp_random_graph(n=10, p=0.5)
    graph = create_graph(dir_path,path_splitting)

    # Precompute probabilities and generate walks
    node2vec = Node2Vec(graph, dimensions=100, walk_length=44, num_walks=300, workers=4)

    # Embed
    # Any keywords acceptable by gensim.Word2Vec can be passed, `diemnsions` and `workers` are automatically passed
    #  (from the Node2Vec constructor)

    model = node2vec.fit(window=10, min_count=1,batch_words=4)

    # Look for most similar nodes
    model.wv.most_similar('2')  # Output node names are always strings

    # Save embeddings for later use
    model.wv.save_word2vec_format(EMBEDDING_FILENAME)

    # Save model for later use
    model.save(EMBEDDING_MODEL_FILENAME)

if __name__ == '__main__':


    dir_path = ''
    read_train_save(dir_path)
    # create_graph(dir_path)