#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 30/3/2019 5:17 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : models_prediction_for_recommendation.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

from __future__ import print_function
from keras import backend as K
import numpy as np
from keras.optimizers import Adam
import codecs
import os
from os import listdir
from itertools import islice
from collections import OrderedDict
from keras.models import model_from_json
import keras_metrics
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
K.set_image_dim_ordering('th')

seed = 1337
np.random.seed(1337)


batch_size = 64
nb_classes = 2
nb_epoch = 20 # 12,15,20,30

# The dimension of the input image, in this case MNIST image, is 5*128
img_rows, img_cols = 5, 128
# Size of convolution kernels used in the convolution layer
nb_filters = (32,64) #(32,64)
# The scope of the pooling layer operation
pool_sizes = (2,2)
# The size of the convolution kernel
kernel_sizes = (1,5) # (1,5)
INIT_LR = 1e-3
strides = (2,2)
chanDim = 1
EPOCHS = 15

map_characters = {0: 'negative', 1: 'positive'}

# Read in the target data set and prepare the training data
def prepare_dataset_for_predicting(dir_path):

    # 读入全部蛋白质质及其相应特征向量
    embedding_data_dir = dir_path + '/'
    file_names = [f for f in listdir(embedding_data_dir) if f.endswith('csv')]
    feature_type_protein_symbol_embedding = OrderedDict()
    for file_name in file_names:
        key_name = file_name.strip().split('_')[2].strip()
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
    label_data_dir = dir_path + 'causal_genes_prediction/'
    file_names = [f for f in listdir(label_data_dir) if f.startswith('total_train')]
    classification_protein_symbol_label = {}
    for file_name in file_names:
        key_name = file_name.strip().split('_')[4].strip()
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


    # 准备需要预测的蛋白质集合及其特征
    dataset_x_for_predicting = []
    protein_symbols_for_predicting = []
    for protein_symbol in final_target_dataset_including_protein_symbols:
        if protein_symbol in protein_symbols_labeled:
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
    output_fp = dir_path + 'causal_genes_prediction/labeled_causal_genes_by_model_predicting.csv'
    positive_causal_genes_set_from_prediction = []
    if not os.path.exists(output_fp):
        model_file_name = {0: ['1550416011_weights_model.h5', '1550416011_serialize_model.json'],
                           1: ['1550416674_weights_model.h5',
                               '1550416674_serialize_model.json'],
                           2: ['1550417332_weights_model.h5', '1550417332_serialize_model.json'],
                           3: ['1550417988_weights_model.h5', '1550417988_serialize_model.json']}

        results_output = []
        for fold_num in range(len(model_file_name)):
            if K.image_data_format() == 'channels_first':
                X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
            else:
                X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)
            X_test = X_test.astype('float32')
            model_path = dir_path + 'causal_genes_prediction/4fold_cross/'
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
                if score >= 3:
                    positive_causal_genes_set_from_prediction.append(protein_symbol)
                output_file.write(str(protein_symbol) + ',' + str(score) + '\n')
    else:
        with codecs.open(output_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                if len(temp) != 2:
                    continue
                if int(temp[1].strip()) >= 3:
                    positive_causal_genes_set_from_prediction.append(temp[0].strip())
    print('positive_causal_genes_set_from_prediction: ' + str(len(positive_causal_genes_set_from_prediction)))

    return positive_causal_genes_set_from_prediction

# Read in the target data set and prepare the training data
def prepare_dataset_for_train():
    dir_path = ''
    embedding_data_dir = dir_path + 'data_embedding/'
    label_data_dir = dir_path + 'data_label/'

    # 读入多标签
    file_for_label = label_data_dir + 'protein_multi_labels.csv'
    protein_symbol_multi_labels = {}
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) != 2:
                continue
            protein_symbol_multi_labels[temp[0].strip()] = [int(x) for x in temp[1].split(',')]
    print('protein_symbol_multi_labels: ' + str(len(protein_symbol_multi_labels.keys())))
    protein_symbol_multi_labeled = set(protein_symbol_multi_labels.keys())

    # 读入样本特征，过滤掉多余的
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
                if temp[0].strip() not in protein_symbol_multi_labeled:
                    continue
                feature_type_protein_symbol_embedding[key_name][temp[0].strip()] = temp[1].strip().split(',')
    print('feature_type_protein_symbol_embedding: ' + str(len(feature_type_protein_symbol_embedding.keys())))

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
            final_target_dataset_including_protein_symbols = final_target_dataset_including_protein_symbols.intersection(
                temp_protein_symbols)
    print('final_target_dataset_including_protein_symbols: ' + str(len(final_target_dataset_including_protein_symbols)))

    train_test_dataset_x = []
    train_test_dataset_y = []
    for protein_symbol in final_target_dataset_including_protein_symbols:
        temp_features_embedding_array = []
        for feature in feature_type_protein_symbol_embedding.keys():
            feature_embedding = [float(x) for x in feature_type_protein_symbol_embedding[feature][protein_symbol]]
            temp_features_embedding_array.append(feature_embedding)

        train_test_dataset_x.append(temp_features_embedding_array)
        train_test_dataset_y.append(protein_symbol_multi_labels[protein_symbol])
    X = np.array(train_test_dataset_x)
    Y = np.array(train_test_dataset_y)

    # 归一化
    for i, brain_slice in enumerate(X):
        brain_slice = (brain_slice - np.mean(brain_slice)) / np.std(brain_slice)
        # brain_slice = (brain_slice - np.min(brain_slice)) / (np.max(brain_slice)-np.min(brain_slice))
        # 下面的if...else很关键，如果没有这个叠加操作，你会发现for循环结束后imgs里面的数据还是未归一化的数据
        if i == 0:
            X = np.reshape(brain_slice, [1, brain_slice.shape[0], brain_slice.shape[1]])
        else:
            X = np.concatenate((X, np.reshape(brain_slice, [1, brain_slice.shape[0], brain_slice.shape[1]])),
                               axis=0)
    print('X shape:', X.shape)

    return X, Y


# multi label model prediction
def multi_label_classification_model_predicting(dir_path):

    positive_causal_genes_set_from_prediction = single_label_classification_model_predicting(dir_path)
    dataset_x_for_predicting, protein_symbols_for_predicting = prepare_dataset_for_predicting(dir_path)
    X_for_predicting = []
    for i in range(len(protein_symbols_for_predicting)):
        protein_symbol = protein_symbols_for_predicting[i]
        if protein_symbol not in positive_causal_genes_set_from_prediction:
            continue
        X_for_predicting.append(dataset_x_for_predicting[i])
    X_predicting = np.array(X_for_predicting)

    # 归一化
    for i, brain_slice in enumerate(X_predicting):
        brain_slice = (brain_slice - np.mean(brain_slice)) / np.std(brain_slice)
        # brain_slice = (brain_slice - np.min(brain_slice)) / (np.max(brain_slice)-np.min(brain_slice))
        # 下面的if...else很关键，如果没有这个叠加操作，你会发现for循环结束后imgs里面的数据还是未归一化的数据
        if i == 0:
            X_predicting = np.reshape(brain_slice, [1, brain_slice.shape[0], brain_slice.shape[1]])
        else:
            X_predicting = np.concatenate((X_predicting, np.reshape(brain_slice, [1, brain_slice.shape[0], brain_slice.shape[1]])),
                               axis=0)
    print('X_predicting shape:', X_predicting.shape)

    model_file_name = {0: ['1551097132_weights_model.h5', '1551097132_serialize_model.json'],
                       1: ['1551097269_weights_model.h5','1551097269_serialize_model.json'],
                       2: ['1551097407_weights_model.h5', '1551097407_serialize_model.json'],
                       3: ['1551097544_weights_model.h5', '1551097544_serialize_model.json'],
                       4: ['1551097682_weights_model.h5', '1551097682_serialize_model.json']}

    # 由于没有保留每个模型预测的概率阈值，因此，拿训练好的模型预测有标签的数据，获得概率阈值
    nb_classes = 8
    X, Y = prepare_dataset_for_train()
    kf = KFold(n_splits=len(model_file_name))
    fold_k = 0
    output_results = []
    for train_index, test_index in kf.split(X):

        _, X_test = X[train_index], X[test_index]
        _, y_test = Y[train_index], Y[test_index]

        if K.image_data_format() == 'channels_first':
            X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
        else:
            X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)

        X_test = X_test.astype('float64')
        model_path = dir_path + 'pathological_processes_prediction/5-fold/'
        json_file = open(model_path + model_file_name[fold_k][1], 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights(model_path + model_file_name[fold_k][0])
        adam = Adam(lr=INIT_LR, decay=INIT_LR / EPOCHS)
        loaded_model.compile(loss='binary_crossentropy', optimizer=adam, metrics=['accuracy'])
        y_probs = np.array(loaded_model.predict_proba(X_test))

        threshold = np.arange(0.001, 0.9, 0.001) #0.1, 0.9, 0.1
        acc = []
        accuracies = []
        best_threshold = np.zeros(y_probs.shape[1])
        for i in range(y_probs.shape[1]):
            y_prob = np.array(y_probs[:, i])
            for j in threshold:
                y_pred = [1 if prob >= j else 0 for prob in y_prob]
                # 马修斯相关系数是在使用机器学习作为二进制（2类）的质量的度量的分类
                acc.append(matthews_corrcoef(y_test[:, i], y_pred))
            acc = np.array(acc)
            index = np.where(acc == acc.max())
            accuracies.append(acc.max())
            best_threshold[i] = threshold[index[0][0]]
            acc = []
        print("best thresholds", best_threshold)

        # 开始预测
        X_predicting = X_predicting.astype('float64')
        if K.image_data_format() == 'channels_first':
            X_predicting = X_predicting.reshape(X_predicting.shape[0], 1, img_rows, img_cols)
        else:
            X_predicting = X_predicting.reshape(X_predicting.shape[0], img_rows, img_cols, 1)
        y_probs_predicting = np.array(loaded_model.predict_proba(X_predicting))
        print(y_probs_predicting[0])
        y_pred_predicting = np.array([[1 if y_probs_predicting[i, j] >= best_threshold[j] else 0 for j in range(nb_classes)] for i in
                           range(len(y_probs_predicting))])
        print("y_pred_predicting", y_pred_predicting)
        output_results.append(y_pred_predicting)
        fold_k += 1

    protein_symbol_pathological_processes = {}
    for p in range(len(positive_causal_genes_set_from_prediction)):
        sum_protein_vectors = np.zeros((nb_classes))
        for k in range(len(output_results)):
            sum_protein_vectors += output_results[k][p]
        protein_symbol_pathological_processes[positive_causal_genes_set_from_prediction[p]] = sum_protein_vectors

    output_fp = dir_path + 'pathological_processes_prediction/labeled_pathological_processes_by_model_predicting.csv'
    with open(output_fp, 'w') as output_file:
        for protein_symbol in protein_symbol_pathological_processes.keys():
            classes_string = ','.join([str(int(cl)) for cl in protein_symbol_pathological_processes[protein_symbol]])
            output_file.write(str(protein_symbol) + ',' + classes_string + '\n')


# Prepare etiological gene recommendation classification
def prepare_causal_genes_recommendation_results(dir_path):
    recommendation_results_dp = dir_path + 'causal_genes_prediction/'

    # 五颗星：基于疾病与基因关联数据库的阳性标签
    file_for_label = recommendation_results_dp + 'total_gene_symbols_labeled_by_disease_name_and_keywords.csv'
    five_stars_protein_symbols = set()
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            five_stars_protein_symbols.add(temp)
    print('five_stars_protein_symbols: '+str(len(five_stars_protein_symbols)))

    # 四颗星： 基于小鼠表型与基因关联的阳性标签
    file_for_label = recommendation_results_dp + 'total_train_dataset_with_positive_samples_by_disease_phenotype.csv'
    four_stars_protein_symbols = set()
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            if temp in five_stars_protein_symbols:
                continue
            four_stars_protein_symbols.add(temp)
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
            if temp[1].strip() == '3':
                two_stars_protein_symbols.add(temp[0].strip())
            elif temp[1].strip() == '4':
                three_stars_protein_symbols.add(temp[0].strip())
            elif temp[1].strip() == '2':
                half_star_protein_symbols.add(temp[0].strip())
            else:
                negtive_protein_symbols.add(temp[0].strip())
    print('three_stars_protein_symbols: ' + str(len(three_stars_protein_symbols)))
    print('two_stars_protein_symbols: ' + str(len(two_stars_protein_symbols)))
    print('negtive_protein_symbols: ' + str(len(negtive_protein_symbols)))

    file_for_label = recommendation_results_dp + 'total_train_dataset_with_negative_samples_by_phenotype.csv'
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            negtive_protein_symbols.add(temp)
    print('negtive_protein_symbols: ' + str(len(negtive_protein_symbols)))

    # 不确定: 属于潜在病因基因种子集合，但不属于上述四类
    file_for_label = recommendation_results_dp + 'prepare_gene_symbols_for_train_model_0115.csv'
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

    output_fp = recommendation_results_dp + 'recommendation_results/causal_genes_recommendation_results.csv'
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


# Prepare recommended classification of pathologic processes
# in which etiological genes are involved
def prepare_pathological_processes_recommendation_results(dir_path):

    recommendation_results_dp = dir_path + 'models_predicting_results/pathological_processes_prediction/'

    pathological_processes = [u'male infertility',
                              u'spermatogenesis abnormalities',
                              u'fertilization and early embryonic development',
                              u'abnormal testicular development and/or related diseases',
                              u'pathological types and/or structural abnormalities of sperm',
                              u'underlying syndromes affecting endocrine and/or urogenital systems',
                              u'malignant tumors of the urogenital system',
                              u'abnormal testicular development and/or related diseases']

    # 三颗星：专家标注
    file_for_label = recommendation_results_dp + 'protein_multi_labels.csv'
    three_stars_protein_symbol_pathological_processes = {}
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            pathological_processes_indexes = temp[1].strip().split(',')
            pathological_processes_names = []
            for i in range(len(pathological_processes_indexes)):
                if pathological_processes_indexes[i] == '1':
                    pathological_processes_names.append(pathological_processes[i].strip())
            three_stars_protein_symbol_pathological_processes[temp[0].strip()]=pathological_processes_names
    print('three_stars_protein_symbol_pathological_processes: '+str(len(three_stars_protein_symbol_pathological_processes)))

    # 两颗星：算法预测标注
    file_for_label = recommendation_results_dp + 'labeled_pathological_processes_by_model_predicting.csv'
    two_stars_protein_symbol_pathological_processes = {}
    with codecs.open(file_for_label, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            pathological_processes_indexes = temp[1:]
            pathological_processes_names = []
            for i in range(len(pathological_processes_indexes)):
                if int(pathological_processes_indexes[i]) >= 3:
                    pathological_processes_names.append(pathological_processes[i].strip())
            if len(pathological_processes_names) == 0:
                for i in range(len(pathological_processes_indexes)):
                    if int(pathological_processes_indexes[i]) >= 2:
                        pathological_processes_names.append(pathological_processes[i].strip())
            if len(pathological_processes_names) == 0:
                for i in range(len(pathological_processes_indexes)):
                    if int(pathological_processes_indexes[i]) >= 1:
                        pathological_processes_names.append(pathological_processes[i].strip())
            two_stars_protein_symbol_pathological_processes[temp[0].strip()] = pathological_processes_names
    print('two_stars_protein_symbol_pathological_processes: ' + str(len(two_stars_protein_symbol_pathological_processes)))

    output_fp = recommendation_results_dp + 'recommendation_results/pathological_processes_recommendation_results.csv'
    with open(output_fp, 'a', encoding='utf-8') as output_file:
        for protein_symbol in three_stars_protein_symbol_pathological_processes:
            pathological_processes_set = ','.join(three_stars_protein_symbol_pathological_processes[protein_symbol])
            output_file.write(str(protein_symbol) + '$' + pathological_processes_set + '$' + '3 star' + '\n')
        for protein_symbol in two_stars_protein_symbol_pathological_processes:
            pathological_processes_set = ','.join(two_stars_protein_symbol_pathological_processes[protein_symbol])
            output_file.write(str(protein_symbol) + '$' + pathological_processes_set + '$' + '2 star' + '\n')


if __name__ == '__main__':

    dir_path = ''
    # single_label_classification_model_predicting(dir_path)
    # multi_label_classification_model_predicting(dir_path)
    # prepare_causal_genes_recommendation_results(dir_path)
    prepare_pathological_processes_recommendation_results(dir_path)

