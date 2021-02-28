#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 6/2/2021 11:35 AM
# @Author  : ggguo
# @Site    : SYM PROJECT
# @File    : causal_genes_prediction_for_new_disease.py
# @Software: SYM application case
###############
###Intermediate process code, for user reference only
###############

from __future__ import print_function
from keras.layers import Dense, Conv2D,Dropout, Activation, Flatten,BatchNormalization
from keras.layers import Convolution2D, MaxPooling2D
from keras.utils import np_utils
from keras import backend as K
import matplotlib.pyplot as plt
import argparse
import numpy as np
from keras.models import Sequential
from keras.optimizers import Adam
from keras.layers.advanced_activations import PReLU
from keras.optimizers import SGD, Adadelta, Adagrad
import codecs
from os import listdir
from itertools import islice
from collections import OrderedDict
from sklearn.model_selection import train_test_split
import keras
import time
from keras.models import model_from_json
import scipy.misc
import os
import re
import json
from keras.callbacks import TensorBoard
from keras.preprocessing.image import ImageDataGenerator
from sklearn.model_selection import train_test_split, StratifiedKFold
import numpy
from keras.callbacks import Callback
from sklearn.metrics import f1_score, precision_score, recall_score
from sklearn import preprocessing
import keras_metrics
from keras.models import load_model

seed = 1337
np.random.seed(1337)

batch_size = 64
nb_classes = 2
nb_epoch = 60 # 12,15,20,30

# The dimension of the input image, in this case MNIST image, is 5*128
img_rows, img_cols = 5, 128
# size of convolution kernels used in the convolution layer
nb_filters = (32,64) #(32,64)
# The scope of the pooling layer operation
pool_sizes = (2,2)
# The size of the convolution kernel
kernel_sizes = (1,5) # (1,5)
INIT_LR = 1e-3
strides = (2,2)
chanDim = 1

map_characters = {'negative': 0, 'positive': 1}

# Standardization of gene names
def human_protein_official_symbol_names_standardization(dir_path, gene_symbols):

    official_protein_symbol_fp = dir_path + 'human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    human_protein_official_symbol_names = {}
    with codecs.open(official_protein_symbol_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            human_protein_official_symbol_names[temp[0]] = set(temp[1:])
    print('human_protein_official_symbol_names: ' + str(len(human_protein_official_symbol_names.keys())))

    gene_symbols_standardization = set()
    for gene_name in gene_symbols:
        if gene_name in human_protein_official_symbol_names.keys():
            gene_symbols_standardization.add(gene_name)
            continue
        for key in human_protein_official_symbol_names.keys():
            if gene_name in human_protein_official_symbol_names[key]:
                gene_symbols_standardization.add(gene_name)
                break
    print('gene_symbols_standardization: ' + str(len(gene_symbols_standardization)))

    return gene_symbols_standardization

# Read in the target data set and prepare the training data
def prepare_dataset_for_predicting(dir_path):

    # 读入全部蛋白质质及其相应特征向量
    embedding_data_dir = dir_path + 'data_embedding/'
    file_names = [f for f in listdir(embedding_data_dir) if f.endswith('csv')]
    feature_type_protein_symbol_embedding = OrderedDict()
    for file_name in file_names:
        key_name = file_name.strip().split('_')[3].strip()
        feature_type_protein_symbol_embedding[key_name] = {}
        raw_data_fp = embedding_data_dir + file_name
        with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 1, None):
                temp = line.strip().split('$')
                if len(temp) < 2:
                    continue
                feature_type_protein_symbol_embedding[key_name][temp[0].strip()] = temp[1].strip().split(',')
    print('feature_type_protein_symbol_embedding: ' + str(len(feature_type_protein_symbol_embedding.keys())))

    # 匹配出五种特征都存在的蛋白质集合
    final_target_dataset_including_protein_symbols = set()
    feature_types = list(feature_type_protein_symbol_embedding.keys())
    for i in range(len(feature_types)):
        key = feature_types[i]
        temp_protein_symbols = set()
        for protein_symbol in feature_type_protein_symbol_embedding[key]:
            temp_protein_symbols.add(protein_symbol)
        if i == 0:
            final_target_dataset_including_protein_symbols = temp_protein_symbols
        else:
            final_target_dataset_including_protein_symbols = final_target_dataset_including_protein_symbols.intersection(temp_protein_symbols)
    print('final_target_dataset_including_protein_symbols: ' + str(len(final_target_dataset_including_protein_symbols)))

    # 读入已经打上标签的蛋白质集合
    label_data_dir = dir_path + 'data_label/'
    file_names = [f for f in listdir(label_data_dir) if f.startswith('gene_names')]
    classification_protein_symbol_label = {}
    for file_name in file_names:
        key_name = file_name.strip().split('_')[3].strip()
        classification_protein_symbol_label[key_name] = []
        raw_data_fp = label_data_dir + file_name
        with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip()
                if temp == '':
                    continue
                if temp not in final_target_dataset_including_protein_symbols:
                    continue
                classification_protein_symbol_label[key_name].append(temp)
    protein_symbols_labeled = set(classification_protein_symbol_label['positive']).union(classification_protein_symbol_label['negative'])
    print('classification_protein_symbol_label: ' + str(len(classification_protein_symbol_label.keys())))
    print('positive_classification_protein_symbol_label: ' + str(len(classification_protein_symbol_label['positive'])))
    print('negative_classification_protein_symbol_label: ' + str(len(classification_protein_symbol_label['negative'])))
    print('protein_symbols_labeled: ' + str(len(protein_symbols_labeled)))

    # 读入女性不孕的疾病基因域
    seed_universe_gene_symbols_dir = dir_path + 'seed_universe_gene_symbols_female_infertility.csv'
    seed_universe_gene_symbols = set()
    with codecs.open(seed_universe_gene_symbols_dir, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            if temp == '':
                continue
            seed_universe_gene_symbols.add(temp)
    print('seed_universe_gene_symbols: ' + str(len(seed_universe_gene_symbols)))

    # 准备需要预测的蛋白质集合及其特征
    dataset_x_for_predicting = []
    protein_symbols_for_predicting = []
    for protein_symbol in final_target_dataset_including_protein_symbols:
        if protein_symbol in protein_symbols_labeled:
            continue
        if protein_symbol not in seed_universe_gene_symbols:
            continue
        temp_features_embedding_array = []
        for feature in feature_type_protein_symbol_embedding.keys():
            feature_embedding = [float(x) for x in feature_type_protein_symbol_embedding[feature][protein_symbol]]
            temp_features_embedding_array.append(feature_embedding)
        dataset_x_for_predicting.append(temp_features_embedding_array)
        protein_symbols_for_predicting.append(protein_symbol)
    print('protein_symbols_for_predicting: ' + str(len(protein_symbols_for_predicting)))
    print('dataset_x_for_predicting: ' + str(len(dataset_x_for_predicting)))
    return dataset_x_for_predicting,protein_symbols_for_predicting

# model prediction
def single_label_classification_model_predicting(dir_path):

    # 采用四折交叉训练的结果
    dataset_x_for_predicting,protein_symbols_for_predicting = prepare_dataset_for_predicting(dir_path)
    X_test = np.array(dataset_x_for_predicting)
    print('X_test shape: ' + str(X_test.shape))
    output_fp = dir_path + 'recommendation_results/labeled_causal_genes_by_model_predicting.csv'
    positive_causal_genes_set_from_prediction = []
    if not os.path.exists(output_fp):
        model_file_name = {0: ['1612105399_weights_model.h5', '1612105399_serialize_model.json'],
                           1: ['1612105830_weights_model.h5', '1612105830_serialize_model.json'],
                           2: ['1612106269_weights_model.h5', '1612106269_serialize_model.json'],
                           3: ['1612106713_weights_model.h5', '1612106713_serialize_model.json'],
                           4: ['1612107160_weights_model.h5', '1612107160_serialize_model.json'],
                           5: ['1612107608_weights_model.h5', '1612107608_serialize_model.json'],
                           6: ['1612108046_weights_model.h5', '1612108046_serialize_model.json'],
                           7: ['1612108489_weights_model.h5', '1612108489_serialize_model.json'],
                           8: ['1612108933_weights_model.h5', '1612108933_serialize_model.json']}

        results_output = []
        for fold_num in range(len(model_file_name)):
            if K.image_data_format() == 'channels_first':
                X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
            else:
                X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)
            X_test = X_test.astype('float32')
            model_path = dir_path + '9fold_cross/'
            json_file = open(model_path + model_file_name[fold_num][1], 'r')
            loaded_model_json = json_file.read()
            json_file.close()
            loaded_model = model_from_json(loaded_model_json)
            loaded_model.load_weights(model_path + model_file_name[fold_num][0])
            opt = Adam(lr=INIT_LR, decay=INIT_LR / nb_epoch)
            loaded_model.compile(loss='binary_crossentropy',
                                 optimizer=opt, metrics=['accuracy', keras_metrics.precision(), keras_metrics.recall()])
            predicting_output = loaded_model.predict_classes(X_test, verbose=1)
            results_output.append(predicting_output)

        # 最总四个模型结果，按照“少数服从多数的原则”统计，采用保守预测，即将持平的结果预测为0
        with open(output_fp, 'w') as output_file:
            for i in range(len(protein_symbols_for_predicting)):
                protein_symbol = protein_symbols_for_predicting[i]
                score = 0
                for j in range(len(results_output)):
                    score += int(results_output[j][i])
                if score >= 6:
                    positive_causal_genes_set_from_prediction.append(protein_symbol)
                output_file.write(str(protein_symbol) + ',' + str(score) + '\n')
    else:
        with codecs.open(output_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                if len(temp) != 2:
                    continue
                if int(temp[1].strip()) >= 6:
                    positive_causal_genes_set_from_prediction.append(temp[0].strip())
    print('positive_causal_genes_set_from_prediction: ' + str(len(positive_causal_genes_set_from_prediction)))

    return positive_causal_genes_set_from_prediction

# Prepare etiological gene recommendation classifications
def prepare_causal_genes_recommendation_results(dir_path):
    recommendation_results_dp = dir_path + 'recommendation_results/'

    # Five stars: a positive label based on a disease and gene association database
    file_for_label = dir_path + 'data_label/positive_samples_from_disease_db.csv'
    positive_samples_from_disease = set()
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            positive_samples_from_disease.add(temp[0].strip())

    five_stars_protein_symbols = human_protein_official_symbol_names_standardization(dir_path, positive_samples_from_disease)
    print('five_stars_protein_symbols: '+str(len(five_stars_protein_symbols)))

    # Four stars: a positive label based on a mouse phenotype associated with a gene
    file_for_label = dir_path + 'data_label/positive_samples_from_phenotype_db.csv'
    positive_samples_from_phenotype = set()
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if temp[0].strip() in five_stars_protein_symbols:
                continue
            positive_samples_from_phenotype.add(temp[0].strip())

    four_stars_protein_symbols = human_protein_official_symbol_names_standardization(dir_path, positive_samples_from_phenotype)
    print('four_stars_protein_symbols: ' + str(len(four_stars_protein_symbols)))

    # 三颗星： k 个模型预测全为阳性标签
    # 二颗星： k-1个模型预测为阳性标签
    # 半颗星： 不确定推荐的阳性标签
    # 零颗星： 不推荐的病因基因
    file_for_label = recommendation_results_dp + 'labeled_causal_genes_by_model_predicting.csv'
    three_stars_protein_symbols = set()
    two_stars_protein_symbols = set()
    half_star_protein_symbols = set()
    negtive_protein_symbols = set()
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if temp[0].strip() in five_stars_protein_symbols.union(four_stars_protein_symbols):
                continue
            if temp[1].strip() in ['6','7']:
                two_stars_protein_symbols.add(temp[0].strip())
            elif temp[1].strip() in ['8','9']:
                three_stars_protein_symbols.add(temp[0].strip())
            elif temp[1].strip() == '5':
                half_star_protein_symbols.add(temp[0].strip())
            else:
                negtive_protein_symbols.add(temp[0].strip())
    print('three_stars_protein_symbols: ' + str(len(three_stars_protein_symbols)))
    print('two_stars_protein_symbols: ' + str(len(two_stars_protein_symbols)))
    print('negtive_protein_symbols: ' + str(len(negtive_protein_symbols)))

    file_for_label = dir_path + 'data_label/gene_names_for_negative_example_standardization.csv'
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            negtive_protein_symbols.add(temp)
    print('negtive_protein_symbols: ' + str(len(negtive_protein_symbols)))

    # 不确定: 属于潜在病因基因种子集合，但不属于上述四类
    file_for_label = dir_path + 'seed_universe_gene_symbols_female_infertility.csv'
    prepare_gene_symbols_for_train_model = set()
    two_three_four_five_protein_symbols = two_stars_protein_symbols.union(three_stars_protein_symbols.union(four_stars_protein_symbols.union(five_stars_protein_symbols)))
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            prepare_gene_symbols_for_train_model.add(temp)
            if temp in two_three_four_five_protein_symbols:
                continue
            if temp in negtive_protein_symbols:
                continue
            half_star_protein_symbols.add(temp)
    print('half_star_protein_symbols: ' + str(len(half_star_protein_symbols)))
    print('prepare_gene_symbols_for_train_model: ' + str(len(prepare_gene_symbols_for_train_model)))

    total_protein_symbols = half_star_protein_symbols.union(negtive_protein_symbols.union(two_stars_protein_symbols.
             union(three_stars_protein_symbols.union(four_stars_protein_symbols.union(five_stars_protein_symbols)))))
    print('total_protein_symbols: ' + str(len(total_protein_symbols)))

    output_fp = recommendation_results_dp + 'causal_genes_recommendation_results.csv'
    with open(output_fp, 'a') as output_file:
        for protein_symbol in half_star_protein_symbols:
            output_file.write(str(protein_symbol) + ',' + '0.5 star' + '\n')
        for protein_symbol in negtive_protein_symbols:
            output_file.write(str(protein_symbol) + ',' + '0 star' + '\n')
        for protein_symbol in two_stars_protein_symbols:
            output_file.write(str(protein_symbol) + ',' + '2 stars' + '\n')
        for protein_symbol in three_stars_protein_symbols:
            output_file.write(str(protein_symbol) + ',' + '3 stars' + '\n')
        for protein_symbol in four_stars_protein_symbols:
            output_file.write(str(protein_symbol) + ',' + '4 stars' + '\n')
        for protein_symbol in five_stars_protein_symbols:
            output_file.write(str(protein_symbol) + ',' + '5 stars' + '\n')


if __name__ == '__main__':

    dir_path = ''
    single_label_classification_model_predicting(dir_path)
    # prepare_causal_genes_recommendation_results(dir_path)


