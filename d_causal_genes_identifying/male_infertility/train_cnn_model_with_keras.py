#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 4/1/2019 9:35 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : train_cnn_model_by_keras_v1.py
# @Software: Mining from a specific set of proteins in human sperm

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
from sklearn.cross_validation import train_test_split
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
nb_epoch = 20 # 12,15,20,30

# The dimension of the input image, in this case MNIST image, is 5*128
img_rows, img_cols = 5, 128
# Number of convolution kernels used in the convolution layer
nb_filters = (32,64) #(32,64)
# The scope of the pooling layer operation
pool_sizes = (2,2)
# The size of the convolution kernel
kernel_sizes = (1,5) # (1,5)
INIT_LR = 1e-3
strides = (2,2)
chanDim = 1

map_characters = {0: 'negative', 1: 'positive'}

# Read in the target data set and prepare the training data
def prepare_dataset_for_train(dir_path):

    embedding_data_dir = dir_path + 'data_embedding/'
    label_data_dir = dir_path + 'data_label/'
    file_names = [f for f in listdir(embedding_data_dir) if f.endswith('csv')]
    feature_type_protein_symbol_embedding = OrderedDict()
    for file_name in file_names:
        if 'Copy' in file_name:
            continue
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

    file_names = [f for f in listdir(label_data_dir) if f.endswith('csv')]
    classification_protein_symbol_label = {}
    for file_name in file_names:
        if 'Copy' in file_name:
            continue
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
    print('classification_protein_symbol_label: ' + str(len(classification_protein_symbol_label.keys())))
    print('positive_classification_protein_symbol_label: ' + str(len(classification_protein_symbol_label['positive'])))
    print('negative_classification_protein_symbol_label: ' + str(len(classification_protein_symbol_label['negative'])))

    # # 截取阴性例子，使其与阳性例子对等
    # for key in classification_protein_symbol_label.keys():
    #     if key == 'negative':
    #         positive_count = len(classification_protein_symbol_label['positive'])
    #         classification_protein_symbol_label[key] = classification_protein_symbol_label[key][:positive_count]

    train_protein_symbol_label = {}
    for key in classification_protein_symbol_label.keys():
        if key == 'positive':
            for protein_symbol in classification_protein_symbol_label[key]:
                train_protein_symbol_label[protein_symbol] = '1'
        elif key == 'negative':
            for protein_symbol in classification_protein_symbol_label[key]:
                train_protein_symbol_label[protein_symbol] = '0'
        else:
            print('please input the correct key word!')
    print('train_protein_symbol_label: ' + str(
                len(train_protein_symbol_label.keys())))

    train_test_dataset_x = []
    train_test_dataset_y = []
    for protein_symbol in train_protein_symbol_label.keys():
        temp_features_embedding_array = []
        for feature in feature_type_protein_symbol_embedding.keys():
            feature_embedding = [float(x) for x in feature_type_protein_symbol_embedding[feature][protein_symbol]]
            temp_features_embedding_array.append(feature_embedding)

        train_test_dataset_x.append(temp_features_embedding_array)
        train_test_dataset_y.append(train_protein_symbol_label[protein_symbol])

    return train_test_dataset_x,train_test_dataset_y

# building model
def build_cnn_two_conv_model(input_shape):

    # 建立序贯模型
    model = Sequential()

    # 卷积层，对二维输入进行滑动窗卷积
    # 当使用该层为第一层时，应提供input_shape参数，在tf模式中，通道维位于第三个位置
    # border_mode：边界模式，为"valid","same"或"full"，即图像外的边缘点是补0
    # 还是补成相同像素，或者是补1
    model.add(Convolution2D(nb_filters[0], kernel_sizes,
                            padding='same',
                            dim_ordering='tf', #th
                            input_shape=input_shape))
    model.add(BatchNormalization(axis=chanDim))
    model.add(Activation('relu')) # tanh,relu

    # 池化层，选用Maxpooling，给定pool_size，dropout比例为0.25
    model.add(MaxPooling2D(pool_size=pool_sizes, strides=strides, padding='same'))
    model.add(Dropout(0.75))

    # 卷积层，激活函数是ReLu
    model.add(Convolution2D(nb_filters[1], kernel_sizes,padding='same'))
    model.add(BatchNormalization(axis=chanDim))
    model.add(Activation('relu')) # tanh,relu

    # 池化层，选用Maxpooling，给定pool_size，dropout比例为0.25
    model.add(MaxPooling2D(pool_size=pool_sizes,strides=strides, padding='same'))
    model.add(Dropout(0.75))

    # Flatten层，把多维输入进行一维化，常用在卷积层到全连接层的过渡
    model.add(Flatten())

    # 包含128个神经元的全连接层，激活函数为ReLu，dropout比例为0.5
    model.add(Dense(64)) # 64, 128,512,1024
    model.add(BatchNormalization())
    model.add(Activation('relu')) #tanh,relu
    model.add(Dropout(0.5))

    # 包含2个神经元的输出层
    # 习惯性的认为，SigmoidCrossEntropyLoss 用于二类问题；
    # SoftmaxCrossEntropyLoss 用于多类问题. 但，在二分类情况时，
    # SoftmaxCrossEntropyLoss 与 SigmoidCrossEntropyLoss 作用等价.
    # 但是，Softmax 会比 Sigmoid 浪费 2 倍的权值空间
    model.add(Dense(nb_classes))
    model.add(Activation('sigmoid')) # sigmoid,softmax

    # 输出模型的参数信息
    model.summary()
    # 配置模型的学习过程
    opt = Adam(lr=INIT_LR, decay=INIT_LR / nb_epoch)
    model.compile(loss='binary_crossentropy', #squared_hinge,categorical_crossentropy
                  optimizer=opt,  # RMSprop,Adagrad,adadelta,Adam
                  metrics=['accuracy',keras_metrics.precision(), keras_metrics.recall()])
    return model

# Training output
def training_process(dir_path,X_train,X_test,Y_train,Y_test,kfold,fold_num):

    if K.image_data_format() == 'channels_first':
        X_train = X_train.reshape(X_train.shape[0], 1, img_rows, img_cols)
        X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
        input_shape = (1, img_rows, img_cols)
    else:
        X_train = X_train.reshape(X_train.shape[0], img_rows, img_cols, 1)
        X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)
        input_shape = (img_rows, img_cols, 1)

    # 将X_train, X_test的数据格式转为float32
    X_train = X_train.astype('float32')
    X_test = X_test.astype('float32')

    # 打印出相关信息
    print('X_train shape:', X_train.shape)
    print(X_train.shape[0], 'train samples')
    print(X_test.shape[0], 'test samples')

    # 将类别向量(从0到nb_classes的整数向量)映射为二值类别矩阵，相当于将向量用one-hot重新编码
    Y_train = np_utils.to_categorical(Y_train, nb_classes)
    Y_test = np_utils.to_categorical(Y_test, nb_classes)

    model = build_cnn_two_conv_model(input_shape)

    # 在训练期间运行并输出可用于张量板的文件

    save_dir = dir_path + 'results/' + str(kfold) + 'fold_cross/'
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    completed_time = str(int(time.time()))
    print('log_folder_name: ' + completed_time)
    tbCallBack = TensorBoard(log_dir=save_dir + completed_time,  # log 目录
                             histogram_freq=0,  # 按照何等频率（epoch）来计算直方图，0为不计算 nb_epoch
                             # batch_size=batch_size,     # 用多大量的数据计算直方图
                             write_graph=True,  # 是否存储网络结构图
                             write_grads=True,  # 是否可视化梯度直方图
                             write_images=True,  # 是否可视化参数
                             embeddings_freq=0,
                             embeddings_layer_names=True,
                             embeddings_metadata=True)
    # metrics = Metrics()
    # 训练模型
    # 数据经过随机打乱shuffle=True。verbose=1，训练过程中输出的信息，0、1、2三种方式都可以，无关紧要。
    # show_accuracy=True，训练时每一个epoch都输出accuracy。
    H = model.fit(X_train, Y_train, batch_size=batch_size, epochs=nb_epoch, shuffle=True,
                  verbose=2, validation_data=(X_test, Y_test), callbacks=[tbCallBack])
    train_out_log_path = save_dir + str(fold_num) + '_model_fit_history'
    with open(train_out_log_path, 'w') as f:
        f.write(str(H.history))

    # plot the training loss and accuracy
    plt.style.use('ggplot')
    plt.figure()
    N = nb_epoch
    plt.plot(np.arange(0, N), H.history["loss"], label="train_loss")
    plt.plot(np.arange(0, N), H.history["val_loss"], label="val_loss")
    plt.plot(np.arange(0, N), H.history["acc"], label="train_acc")
    plt.plot(np.arange(0, N), H.history["val_acc"], label="val_acc")
    plt.title("Training Loss and Accuracy")
    plt.xlabel("Epoch #")
    plt.ylabel("Loss/Accuracy")
    plt.legend(loc="upper left")
    image_path = save_dir + completed_time + '_loss_accu_train_var.png'
    plt.savefig(image_path)
    plt.show()

    # 按batch计算在某些输入数据上模型的误差
    scores = model.evaluate(X_test, Y_test, verbose=1)
    print(scores)

    # 输出训练好的模型在测试集上的表现
    print('Test loss score:', scores[0])
    print('Test accuracy:', scores[1])

    # serialize weights to HDF5
    model_path = save_dir + completed_time
    json_model_path = model_path + '_serialize_model.json'
    weights_model_path = model_path + '_weights_model.h5'

    # serialize model to JSON
    model_json = model.to_json()
    with open(json_model_path, "w") as json_file:
        json_file.write(model_json)
    # serialize weights to HDF5
    model.save_weights(weights_model_path)
    print("Saved model to disk")

    return scores

# training model
def train_cnn_model(dir_path,kfold =7):

    # K-fold cross 训练模型
    train_test_dataset_x, train_test_dataset_y = prepare_dataset_for_train(dir_path)
    X = np.array(train_test_dataset_x)
    Y = np.array(train_test_dataset_y)

    # # 归一化
    # for i, brain_slice in enumerate(X):
    #     brain_slice = (brain_slice - np.mean(brain_slice)) / np.std(brain_slice)
    #     # 下面的if...else很关键，如果没有这个叠加操作，你会发现for循环结束后X里面的数据还是未归一化的数据
    #     if i == 0:
    #         X = np.reshape(brain_slice, [1, brain_slice.shape[0], brain_slice.shape[1]])
    #     else:
    #         X = np.concatenate((X, np.reshape(brain_slice, [1, brain_slice.shape[0], brain_slice.shape[1]])),
    #                            axis=0)
    # print('X shape:', X.shape)

    # 分层采样，确保训练集，测试集中各类别样本的比例与原始数据集中相同。
    kfold_cross = StratifiedKFold(n_splits=kfold, shuffle=True, random_state=seed).split(X, Y)
    statistic_loss_list = []
    statistic_accuracy_list = []
    statistic_precision_list = []
    statistic_recall_list = []
    # statistic_f1_list = []
    for fold_num, (train_idx, test_idx) in enumerate(kfold_cross):

        print('\nFold ', fold_num)
        X_train = X[train_idx]
        Y_train = Y[train_idx]
        X_test = X[test_idx]
        Y_test = Y[test_idx]

        scores = training_process(dir_path,X_train, X_test, Y_train, Y_test, kfold,fold_num)
        statistic_loss_list.append(scores[0])
        statistic_accuracy_list.append(scores[1])
        statistic_precision_list.append(scores[2])
        statistic_recall_list.append(scores[3])
        # statistic_f1_list.append(scores[4])

    print("%.2f%% (+/- %.2f%%)" % (np.array(statistic_loss_list).mean(), np.array(statistic_loss_list).std()))
    print("%.2f%% (+/- %.2f%%)" % (np.array(statistic_accuracy_list).mean(), np.array(statistic_accuracy_list).std()))

    scores_output_fp = dir_path + 'results/average_test_data_scores_by_kfold.csv'
    with open(scores_output_fp, 'a') as output_file:
        output_file.write(str(kfold) + ',' + 'loss' + ',' + str(np.array(statistic_loss_list).mean()) + ',' + str(np.array(statistic_loss_list).std()) + '\n')
        output_file.write(str(kfold) + ',' + 'accuracy' + ',' + str(np.array(statistic_accuracy_list).mean()) + ',' + str(
            np.array(statistic_accuracy_list).std()) + '\n')
        output_file.write(
            str(kfold) + ',' + 'precision' + ',' + str(np.array(statistic_precision_list).mean()) + ',' + str(
                np.array(statistic_precision_list).std()) + '\n')
        output_file.write(
            str(kfold) + ',' + 'recall' + ',' + str(np.array(statistic_recall_list).mean()) + ',' + str(
                np.array(statistic_recall_list).std()) + '\n')
        # output_file.write(
        #     str(kfold) + ',' + 'f1' + ',' + str(np.array(statistic_f1_list).mean()) + ',' + str(
        #         np.array(statistic_f1_list).std()) + '\n')

# Set image annotations
def annotate_style(rad):
    connectionstyle = "arc3,rad=%s" %rad
    arrowprops = dict(arrowstyle="->",
                      color="c",
                      shrinkA=5, shrinkB=5,
                      patchA=None,
                      patchB=None,
                      connectionstyle=connectionstyle,
                      )
    return arrowprops

# Visually compare the results of different folds
def visualize_k_fold_results(dir_path,kfolds):

    k_fold_results_fp = dir_path + 'results/average_test_data_scores_by_kfold.csv'
    loss_k_folds_mean = []
    accuracy_k_folds_mean = []
    loss_k_folds_std = []
    accuracy_k_folds_std = []
    precision_k_folds_mean = []
    recall_k_folds_mean = []
    precision_k_folds_std = []
    recall_k_folds_std = []
    total_evaluation_indicator_k_folds_mean = []
    total_evaluation_indicator_k_folds_std = []
    with codecs.open(k_fold_results_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if 'loss' in temp:
                loss_k_folds_mean.append(float(temp[2]))
                loss_k_folds_std.append(float(temp[-1]))
            elif 'accuracy' in temp:
                accuracy_k_folds_mean.append(float(temp[2]))
                accuracy_k_folds_std.append(float(temp[-1]))
            elif 'precision' in temp:
                precision_k_folds_mean.append(float(temp[2]))
                precision_k_folds_std.append(float(temp[-1]))
            elif 'recall' in temp:
                recall_k_folds_mean.append(float(temp[2]))
                recall_k_folds_std.append(float(temp[-1]))
            else:
                print('Error.................')

    for k in range(len(kfolds)):
        sum_mean = -loss_k_folds_mean[k]+accuracy_k_folds_mean[k]+precision_k_folds_mean[k]+recall_k_folds_mean[k]
        total_evaluation_indicator_k_folds_mean.append(sum_mean/2.0)
        sum_std = loss_k_folds_std[k] + accuracy_k_folds_std[k] + precision_k_folds_std[k] + recall_k_folds_std[k]
        total_evaluation_indicator_k_folds_std.append(sum_std /3.0)

    # plot the training loss and accuracy
    font = {'family': 'Times New Roman',
            'weight': 'normal',
            'size': 10,
            }

    plt.style.use('ggplot')
    plt.figure()
    N = len(kfolds) + 3
    A, = plt.plot(np.arange(3, N), loss_k_folds_mean, label="K-fold loss functions", linewidth=1)
    B, = plt.plot(np.arange(3, N), accuracy_k_folds_mean, label="K-fold accuracies", linewidth=1)
    C, = plt.plot(np.arange(3, N), precision_k_folds_mean, label="K-fold precisions", linewidth=1)
    D, = plt.plot(np.arange(3, N), recall_k_folds_mean, label="K-fold recalls", linewidth=1)
    E, = plt.plot(np.arange(3, N), total_evaluation_indicator_k_folds_mean, label="Comprehensive index",linestyle=':', linewidth=2)
    # plt.title("Training Loss and Accuracy").set_fontname('Times New Roman')
    plt.xlabel("K-Fold", fontsize=10).set_fontname('Times New Roman')
    bbox = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    arrowprops = annotate_style(0.3)
    plt.annotate('(%s, %.4f)' % ('10-fold', total_evaluation_indicator_k_folds_mean[-1]),
                        xy=(kfolds[-1], total_evaluation_indicator_k_folds_mean[-1]), xycoords='data',  # figure pixels
                        xytext=(0.85 * kfolds[-1], 0.85 * total_evaluation_indicator_k_folds_mean[1]), textcoords='data',  # offset points
                        bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    plt.annotate('(%s, %.4f)' % ('4-fold', total_evaluation_indicator_k_folds_mean[1]),
                        xy=(kfolds[1], total_evaluation_indicator_k_folds_mean[1]), xycoords='data',  # figure pixels
                        xytext=(0.85 * kfolds[1], 0.85 * total_evaluation_indicator_k_folds_mean[1]), textcoords='data',  # offset points
                        bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    # plt.annotate('(%s, %.4f)' % ('8-fold', accuracy_k_folds_mean[-3]),
    #                     xy=(kfolds[-3], accuracy_k_folds_mean[-3]), xycoords='data',  # figure pixels
    #                     xytext=(0.6 * kfolds[-1], 0.96 * accuracy_k_folds_mean[2]), textcoords='data',  # offset points
    #                     bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    #
    # arrowprops = annotate_style(-0.3)
    # plt.annotate('(%s, %.4f)' % ('10-fold', loss_k_folds_mean[-1]),
    #                    xy=(kfolds[-1], loss_k_folds_mean[-1]),xycoords='data', # figure pixels
    #                    xytext=(0.85*kfolds[-1], 1.046 * loss_k_folds_mean[2]),textcoords='data', # offset points
    #                    bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    # plt.annotate('(%s, %.4f)' % ('5-fold', loss_k_folds_mean[2]),
    #                     xy=(kfolds[2], loss_k_folds_mean[2]), xycoords='data',  # figure pixels
    #                     xytext=(0.65 * kfolds[2], 1.046 * loss_k_folds_mean[2]), textcoords='data',  # offset points
    #                     bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    # plt.annotate('(%s, %.4f)' % ('8-fold', loss_k_folds_mean[-3]),
    #                     xy=(kfolds[-3], loss_k_folds_mean[-3]), xycoords='data',  # figure pixels
    #                     xytext=(0.6 * kfolds[-1], 1.046 * loss_k_folds_mean[2]), textcoords='data',
    #                     bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    plt.ylabel("Mean", fontsize=10).set_fontname('Times New Roman')
    # plt.legend(loc="center")
    plt.legend(handles=[A, B, C, D,E], prop=font,loc=2, bbox_to_anchor=(1.01,1.0),borderaxespad = 0.0)
    image_path = dir_path + 'results/k-fold_loss_accu_mean.pdf'
    # fig.subplots_adjust(right=0.98)
    plt.savefig(image_path,bbox_inches='tight')
    plt.show()

    # plot the training loss and accuracy
    plt.style.use('ggplot')
    plt.figure()
    N = len(kfolds) + 3

    A, = plt.plot(np.arange(3, N), loss_k_folds_std, label="K-fold loss functions", linewidth=1)
    B, = plt.plot(np.arange(3, N), accuracy_k_folds_std, label="K-fold accuracies", linewidth=1)
    C, = plt.plot(np.arange(3, N), precision_k_folds_std, label="K-fold precisions", linewidth=1)
    D, = plt.plot(np.arange(3, N), recall_k_folds_std, label="K-fold recalls", linewidth=1)
    E, = plt.plot(np.arange(3, N), total_evaluation_indicator_k_folds_std, label="Comprehensive index",linestyle=':', linewidth=2)
    # plt.title("Training Loss and Accuracy").set_fontname('Times New Roman')
    plt.xlabel("K-Fold", fontsize=10).set_fontname('Times New Roman')
    bbox = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    arrowprops = annotate_style(0.3)
    plt.annotate('(%s, %.4f)' % ('4-fold', total_evaluation_indicator_k_folds_std[1]),
                        xy=(kfolds[1], total_evaluation_indicator_k_folds_std[1]), xycoords='data',  # figure pixels
                        xytext=(0.85 * kfolds[1], 0.7 * total_evaluation_indicator_k_folds_std[1]), textcoords='data',  # offset points
                        bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    # plt.annotate('(%s, %.2f, %.2f)' % ('8-fold', loss_k_folds_std[-3], accuracy_k_folds_std[-3]),
    #                     xy=(kfolds[-3], accuracy_k_folds_std[-3]), xycoords='data',  # figure pixels
    #                     xytext=(0.67 * kfolds[-3], 0.65 * accuracy_k_folds_std[2]), textcoords='data',  # offset points
    #                     bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    plt.annotate('(%s, %.4f)' % ('10-fold', total_evaluation_indicator_k_folds_std[-1]),
                        xy=(kfolds[-1], total_evaluation_indicator_k_folds_std[-1]), xycoords='data',  # figure pixels
                        xytext=(0.85 * kfolds[-1], 0.7 * total_evaluation_indicator_k_folds_std[1]), textcoords='data',  # offset points
                        bbox=bbox, arrowprops=arrowprops,size=10).set_fontname('Times New Roman')
    plt.ylabel("Standard deviation", fontsize=10).set_fontname('Times New Roman')
    # plt.legend(loc="upper left")
    plt.legend(handles=[A, B, C, D, E], prop=font,loc=2, bbox_to_anchor=(1.01,1.0),borderaxespad = 0.0)
    image_path = dir_path + 'results/k-fold_loss_accu_std.pdf'
    # fig.subplots_adjust(right=0.98)
    plt.savefig(image_path,bbox_inches='tight')
    plt.show()

# Load the training composition history file and visualize the training results
def visualize_history_model_results(dir_path):

    history_training_fp = dir_path + 'results/4fold_cross/'
    file_names = [f for f in listdir(history_training_fp) if f.endswith('_model_fit_history')]
    history_training_acc = {}
    history_training_loss = {}
    history_training_val_acc = {}
    history_training_val_loss = {}
    for file_name in file_names:
        json_file = open(history_training_fp+file_name, 'r')
        json_file_str = json_file.read()
        loaded_model_json = re.sub('\'', '\"', json_file_str)
        history_dict = json.loads(loaded_model_json)
        json_file.close()
        key = file_name.split('_')[0].strip()
        print(history_dict.keys())
        history_training_acc[key] = history_dict['acc']
        history_training_val_acc[key] = history_dict['val_acc']
        history_training_loss[key] = history_dict['loss']
        history_training_val_loss[key] = history_dict['val_loss']

    # 带有多个轴域刻度的 plot
    plt.close('all')
    plt.style.use('ggplot')
    plt.figure(1)
    font = {'family': 'Times New Roman',
            'weight': 'normal',
            'size': 10,
            }
    y_ticks = np.arange(0.55, 0.80, 0.05)

    plt.subplot(231)
    epochs = range(1, len(history_training_acc['0']) + 1)
    plt.plot(epochs, history_training_acc['0'])
    plt.plot(epochs, history_training_loss['0'])
    plt.plot(epochs, history_training_val_acc['0'])
    plt.plot(epochs, history_training_val_loss['0'])
    # plt.xticks([])
    # plt.legend(loc="upper center", prop=font)
    plt.xlabel("Epoch", fontsize=10).set_fontname('Times New Roman')
    plt.ylabel("Accuracy/Loss function", fontsize=10).set_fontname('Times New Roman')
    plt.yticks(y_ticks, fontsize=10)
    plt.title('Fold %s' %'0',fontsize=10).set_fontname('Times New Roman')
    plt.grid(True)

    plt.subplot(232)
    epochs = range(1, len(history_training_acc['1']) + 1)
    plt.plot(epochs, history_training_acc['1'], label="Train accuracy")
    plt.plot(epochs, history_training_loss['1'], label="Train loss function")
    plt.plot(epochs, history_training_val_acc['1'], label="Validation accuracy")
    plt.plot(epochs, history_training_val_loss['1'], label="Validation loss function")
    plt.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0., prop=font)
    # plt.xticks([])
    plt.xlabel("Epoch", fontsize=10).set_fontname('Times New Roman')
    plt.ylabel("Accuracy/Loss function", fontsize=10).set_fontname('Times New Roman')
    plt.yticks(y_ticks, fontsize=10)
    plt.title('Fold %s' % '1',fontsize=10).set_fontname('Times New Roman')
    plt.grid(True)

    plt.subplot(234)
    epochs = range(1, len(history_training_acc['2']) + 1)
    plt.plot(epochs, history_training_acc['2'])
    plt.plot(epochs, history_training_loss['2'])
    plt.plot(epochs, history_training_val_acc['2'])
    plt.plot(epochs, history_training_val_loss['2'])
    # plt.xticks([])
    plt.xlabel("Epoch", fontsize=10).set_fontname('Times New Roman')
    plt.ylabel("Accuracy/Loss function", fontsize=10).set_fontname('Times New Roman')
    plt.yticks(y_ticks, fontsize=10)
    plt.title('Fold %s' % '2',fontsize=10).set_fontname('Times New Roman')
    plt.grid(True)

    plt.subplot(235)
    epochs = range(1, len(history_training_acc['3']) + 1)
    plt.plot(epochs, history_training_acc['3'], label="Train accuracy")
    plt.plot(epochs, history_training_loss['3'], label="Train loss function")
    plt.plot(epochs, history_training_val_acc['3'], label="Validation accuracy")
    plt.plot(epochs, history_training_val_loss['3'], label="Validation loss function")
    plt.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0., prop=font)
    # plt.xticks([])
    plt.xlabel("Epoch",fontsize=10).set_fontname('Times New Roman')
    plt.ylabel("Accuracy/Loss function", fontsize=10).set_fontname('Times New Roman')
    plt.yticks(y_ticks,fontsize=10)
    plt.title('Fold %s' % '3',fontsize=10).set_fontname('Times New Roman')
    plt.grid(True)

    plt.tight_layout(pad=0.5)
    image_path = dir_path + 'results/4-fold_train_val_loss_acc.pdf'
    plt.savefig(image_path)

    plt.show()

# Visualize the results of different models
def visualize_models_compared_results(dir_path):

    # DF-CNN-SLC
    model_results = OrderedDict()
    dl_model_results = []
    indicators = ['accuracy','precision','recall','f1 score']
    for i in range(len(indicators)):
        dl_model_results.append([])
    k_fold_results_fp = dir_path + 'results/test_data_scores_by_kfold.csv'
    with codecs.open(k_fold_results_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            for e in range(len(indicators)):
                dl_model_results[e].append(float(temp[e+1]))
    model_results['DP-CNN-SLBC'] = dl_model_results

    # Read in to compare the accuracy and recall rates of the models

    compare_model_results_fp = dir_path + 'compare_models_results/4-fold/'
    file_names = [f for f in listdir(compare_model_results_fp) if f.endswith('_total_results_2019-02-19.txt')]
    print('The amount of files: %d' %len(file_names))
    for file_name in file_names:
        ml_model_results = []
        model_results_fp = compare_model_results_fp + file_name
        algorithm = file_name.split('_')[0].strip()
        with codecs.open(model_results_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('\t')
                ml_model_results.append([float(x) for x in temp[1:]])
        model_results[algorithm] = ml_model_results
    print('model_results: ' + str(len(model_results.keys())))

    # Plot Precision-Recall curve for each model
    font = {'family': 'Times New Roman',
            'weight': 'normal',
            'size': 8,
            }
    plt.style.use('ggplot')
    plt.figure()
    f_scores = np.linspace(0.2, 0.8, num=4)
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
        plt.annotate('F1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02),size=10).set_fontname('Times New Roman')
    markers = ['d','+','o','*','x','1']
    m = 0
    marker_size = [50,40,30,40,30,40]
    for algorithm in model_results.keys():
        # print(model_results[algorithm][2], model_results[algorithm][1])
        plt.scatter(model_results[algorithm][2], model_results[algorithm][1],s=marker_size[m],marker=markers[m],label=algorithm)
        m += 1
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall',fontsize=10).set_fontname('Times New Roman')
    plt.ylabel('Precision',fontsize=10).set_fontname('Times New Roman')
    # plt.title('Extension of Precision-Recall curve to models')
    plt.legend(prop=font,loc="upper left")

    plt.tight_layout(pad=0.5)
    image_path = dir_path + 'results/4-fold_comparison_of_prediction_results_of_different_models.pdf'
    plt.savefig(image_path)
    plt.show()

    # --------------------------------------------------------

    # model_mean_results = {}
    # for al in model_results.keys():
    #     mean_indicators = []
    #     for values in model_results[al]:
    #         narray = numpy.array(values)
    #         mean_indicators.append(narray.mean())
    #     model_mean_results[al] = mean_indicators
    #
    #
    # N = len(indicators)
    # angles = np.linspace(0, 2 * np.pi, N, endpoint=False)  # 设置雷达图的角度，用于平分切开一个圆面
    # angles = np.concatenate((angles, [angles[0]]))  # 为了使雷达图一圈封闭起来
    # plt.close()
    # plt.style.use('ggplot')
    # fig = plt.figure(figsize=(7, 7))  # 设置画布大小
    # ax = fig.add_subplot(111, polar=True)  # 这里一定要设置为极坐标格式
    # sam = ['x-', '-.', '--', '-', '--', '-.']  # 样式
    # s = 0
    # lab = []  # 图例标签名
    # for key in model_mean_results.keys():
    #     values = model_mean_results[key]
    #     # 设置各指标名称
    #     feature = [x.capitalize() for x in indicators]
    #     # 为了使雷达图一圈封闭起来，需要下面的步骤
    #     values = np.concatenate((values, [values[0]]))
    #     ax.plot(angles, values, sam[s], linewidth=1)  # 绘制折线图
    #     ax.fill(angles, values, alpha=0.1)  # 填充颜色
    #     ax.set_thetagrids(angles * 180 / np.pi, feature,fontsize=15,family='Times New Roman')  # 添加每个特征的标签
    #     ax.set_ylim(0.3, 0.75)  # 设置雷达图的范围
    #     ax.grid(linestyle='-', linewidth=0.5)  # 添加网格线
    #     lab.append(key)
    #     s += 1
    # plt.legend(lab, prop=font)
    # image_path = dir_path + 'results/4-fold_comparison_of_prediction_results_of_different_models_total_indicators.pdf'
    # plt.savefig(image_path)  # 保存图片到本地
    # plt.show()  # 显示图形

# Load the model and predict the test set data
def load_model_for_prediction(dir_path,kfold):
    # K-fold cross 训练模型
    train_test_dataset_x, train_test_dataset_y = prepare_dataset_for_train(dir_path)
    X = np.array(train_test_dataset_x)
    Y = np.array(train_test_dataset_y)

    # 分层采样，确保训练集，测试集中各类别样本的比例与原始数据集中相同。
    kfold_cross = StratifiedKFold(n_splits=kfold, shuffle=True, random_state=seed).split(X, Y)

    model_file_name = {0:['1550416011_weights_model.h5','1550416011_serialize_model.json'],1:['1550416674_weights_model.h5',
                  '1550416674_serialize_model.json'],2:['1550417332_weights_model.h5','1550417332_serialize_model.json'],
                  3:['1550417988_weights_model.h5','1550417988_serialize_model.json']}

    for fold_num, (train_idx, test_idx) in enumerate(kfold_cross):

        print('\nFold ', fold_num)
        X_test = X[test_idx]
        Y_test = Y[test_idx]

        if K.image_data_format() == 'channels_first':
            X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
        else:
            X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)

        X_test = X_test.astype('float32')
        print(X_test.shape[0], 'test samples')

        # 将类别向量(从0到nb_classes的整数向量)映射为二值类别矩阵，相当于将向量用one-hot重新编码
        Y_test = np_utils.to_categorical(Y_test, nb_classes)

        model_path = dir_path + 'results/4fold_cross/'
        json_file = open(model_path+model_file_name[fold_num][1], 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights(model_path+model_file_name[fold_num][0])
        opt = Adam(lr=INIT_LR, decay=INIT_LR / nb_epoch)
        loaded_model.compile(loss='binary_crossentropy',
                  optimizer=opt, metrics=['accuracy',keras_metrics.precision(), keras_metrics.recall()])
        scores = loaded_model.evaluate(X_test, Y_test, verbose=1)
        print(scores)

        f1_score = 2* float(scores[2])*float(scores[3]) / (float(scores[2])+float(scores[3]))
        temp_evaluate_scores = scores[1:] + [f1_score]

        scores_output_fp = dir_path + 'results/test_data_scores_by_kfold.csv'
        with open(scores_output_fp, 'a') as output_file:
            temp_values = [str(x) for x in temp_evaluate_scores]
            output_file.write(str(fold_num) + ',' + ','.join(temp_values) + '\n')


if __name__ == '__main__':

    dir_path = ''
    kfolds = [3,4,5,6,7,8,9,10]
    # for kfold in kfolds:
    #     train_cnn_model(dir_path,kfold)
    visualize_k_fold_results(dir_path,kfolds)
    # visualize_history_model_results(dir_path)
    # visualize_models_compared_results(dir_path)
    # kfold = 4
    # train_cnn_model(dir_path, kfold)
    # load_model_for_prediction(dir_path, kfold)

