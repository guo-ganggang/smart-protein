#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 19/1/2019 9:22 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : multi-label_classification_with_dl_by_keras.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import numpy as np
from keras.layers import Dense, Conv2D,Dropout, Activation, Flatten,BatchNormalization
from keras.layers import Convolution2D, MaxPooling2D
from keras.optimizers import SGD
from keras.optimizers import Adam
from keras.callbacks import ModelCheckpoint
from sklearn.metrics import matthews_corrcoef
import codecs
from os import listdir
from itertools import islice
from collections import OrderedDict
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import hamming_loss
# from sklearn.metrics import classification_report
from sklearn.metrics import jaccard_similarity_score
# from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from scipy import interp
from sklearn.cross_validation import train_test_split
from sklearn.model_selection import KFold
from keras.models import Sequential
import matplotlib.pyplot as plt
from itertools import cycle
from keras import backend as K
import keras_metrics
import os
from sklearn import preprocessing
import numpy
from keras.callbacks import Callback
from sklearn.metrics import f1_score, precision_score, recall_score
from keras.callbacks import TensorBoard
import time
from bokeh.io import save
from math import log, sqrt
import pandas as pd
from six.moves import cStringIO as StringIO
from bokeh.plotting import figure, show, output_file
from bokeh.io import export_png

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
K.set_image_dim_ordering('th')

img_rows, img_cols = 5, 128
kernel_sizes = (1,8)
nb_filters = (16,32)
chanDim = 1
pool_sizes = (2,2)
strides = (2,2)
nb_classes = 8
INIT_LR = 1e-3
EPOCHS = 15

# Read in the target data set and prepare the training data
def prepare_dataset_for_train(dir_path):

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

# building model
def build_model_architecture(input_shape):

    model = Sequential()
    model.add(Convolution2D(nb_filters[0], kernel_size=kernel_sizes,padding='same',input_shape=input_shape))
    # 通常在线性向非线性转变时使用，能够保证权重的尺度不变，因为BatchNormalization在激活函数前对输入进行了标准化
    # "channels_first"（通道在前）的 Conv2D 层，这时特征轴是编号为 1 的轴
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(Convolution2D(nb_filters[1], kernel_sizes))
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.005))

    # model.add(Convolution2D(nb_filters[1],kernel_sizes, padding='same'))
    # model.add(BatchNormalization(axis=chanDim))
    # model.add(Activation('relu'))
    # model.add(Convolution2D(nb_filters[1], kernel_sizes))
    # model.add(BatchNormalization(axis=chanDim))
    # model.add(Activation('relu'))
    # model.add(MaxPooling2D(pool_size=(2, 2)))
    # model.add(Dropout(0.15))

    model.add(Flatten())
    model.add(Dense(128))
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(Dropout(0.5))
    model.add(Dense(nb_classes))
    model.add(Activation('sigmoid'))
    # model.summary()

    # # Calculate precision for the k label.
    # precisions = []
    # recalls = []
    # for i in range(nb_classes):
    #     recalls.append(keras_metrics.recall(label=i))
    #     precisions.append(keras_metrics.recall(label=i))

    # let's train the model using SGD + momentum (how original).

    # sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    # model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])

    adam = Adam(lr=INIT_LR, decay=INIT_LR / EPOCHS)
    model.compile(loss='binary_crossentropy', optimizer=adam, metrics=['accuracy'])
    return model

# Visual Training Results
def visualize_training_performance(history,kfold,fold_k):
    font = {'family': 'Times New Roman',
            'weight': 'normal',
            'size': 10,
            }
    acc = history.history['acc']
    val_acc = history.history['val_acc']
    loss = history.history['loss']
    val_loss = history.history['val_loss']

    epochs = range(len(acc))
    plt.figure()
    plt.plot(epochs, acc, 'c-', label='Training acc')
    plt.plot(epochs, val_acc, 'm', label='Validation acc')
    # plt.title('Training and validation accuracy')
    plt.xlabel('Epoch', fontsize=10).set_fontname('Times New Roman')
    plt.ylabel('Accuracy', fontsize=10).set_fontname('Times New Roman')
    plt.legend(loc="upper right", prop=font)
    image_path = dir_path + 'results/' + str(kfold) + '-fold/' + 'accuracy_performance_in_training_' + str(fold_k) + '.png'
    plt.savefig(image_path, bbox_inches='tight')

    plt.figure()

    plt.plot(epochs, loss, 'c-', label='Training loss')
    plt.plot(epochs, val_loss, 'm', label='Validation loss')
    # plt.title('Training and validation loss')
    plt.xlabel('Epoch', fontsize=10).set_fontname('Times New Roman')
    plt.ylabel('Loss', fontsize=10).set_fontname('Times New Roman')
    plt.legend(loc="upper right", prop=font)
    image_path = dir_path + 'results/' + str(kfold) + '-fold/'+'loss_performance_in_training_'+str(fold_k)+'.png'
    plt.savefig(image_path, bbox_inches='tight')
    # plt.show()

# training model
def train_model(dir_path,kfold):

    model_checkpoint_path = dir_path + 'results/' + str(kfold) + '-fold/'
    if not os.path.isdir(model_checkpoint_path):
        os.makedirs(model_checkpoint_path)

    X, Y = prepare_dataset_for_train(dir_path)
    kf = KFold(n_splits=kfold)
    fold_k = 0
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = Y[train_index], Y[test_index]

        if K.image_data_format() == 'channels_first':
            X_train = X_train.reshape(X_train.shape[0], 1, img_rows, img_cols)
            X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
            input_shape = (1, img_rows, img_cols)
        else:
            X_train = X_train.reshape(X_train.shape[0], img_rows, img_cols, 1)
            X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)
            input_shape = (img_rows, img_cols, 1)

        # 将X_train, X_test的数据格式转为float32
        X_train = X_train.astype('float64')
        X_test = X_test.astype('float64')

        # 打印出相关信息
        print('X_train shape:', X_train.shape)
        print(X_train.shape[0], 'train samples')
        print(X_test.shape[0], 'test samples')

        print(X_train.shape)
        print(X_test.shape)
        print(y_train.shape)
        print(y_test.shape)

        model = build_model_architecture(input_shape)
        completed_time = str(int(time.time()))
        print('log_folder_name: ' + completed_time)
        tbCallBack = TensorBoard(log_dir=model_checkpoint_path + completed_time,  # log 目录
                                 histogram_freq=0,  # 按照何等频率（epoch）来计算直方图，0为不计算 nb_epoch
                                 # batch_size=batch_size,     # 用多大量的数据计算直方图
                                 write_graph=True,  # 是否存储网络结构图
                                 write_grads=True,  # 是否可视化梯度直方图
                                 write_images=True,  # 是否可视化参数
                                 embeddings_freq=0,
                                 embeddings_layer_names=True,
                                 embeddings_metadata=True)
        # check = ModelCheckpoint(model_checkpoint_path+"weights.{epoch:02d}-{val_acc:.5f}.hdf5", monitor='val_acc', verbose=1, save_best_only=True, save_weights_only=True, mode='auto')
        history = model.fit(X_train, y_train, batch_size=32, epochs=EPOCHS,callbacks=[tbCallBack],verbose=2,validation_data=(X_test,y_test))
        train_out_log_path = model_checkpoint_path + str(fold_k) + '_model_fit_history'
        with open(train_out_log_path, 'w') as f:
            f.write(str(history.history))
            # serialize weights to HDF5

        visualize_training_performance(history,kfold,fold_k)
        y_probs = np.array(model.predict_proba(X_test))
        # serialize model to JSON
        model_path = model_checkpoint_path + completed_time
        json_model_path = model_path + '_serialize_model.json'
        weights_model_path = model_path + '_weights_model.h5'
        model_json = model.to_json()
        with open(json_model_path, "w") as json_file:
            json_file.write(model_json)
        # serialize weights to HDF5
        model.save_weights(weights_model_path)
        print("Saved model to disk")

        average_precision_micro = average_precision_score(y_test, y_probs, average="micro")
        # average_precision_macro = average_precision_score(y_test, y_probs, average="macro")
        auc_micro = roc_auc_score(y_test, y_probs, average='micro')
        # auc_macro = roc_auc_score(y_test, y_probs, average='macro')

        print('average_precision_micro: %.4f' % average_precision_micro)
        # print('average_precision_macro: %.4f' % average_precision_macro)
        print('roc_auc_score_micro: %.4f' % auc_micro)
        # print('roc_auc_score_macro: %.4f' % auc_macro)

        threshold = np.arange(0.1, 0.9, 0.1)
        acc = []
        accuracies = []
        best_threshold = np.zeros(y_probs.shape[1])
        for i in range(y_probs.shape[1]):
            y_prob = np.array(y_probs[:, i])
            for j in threshold:
                y_pred = [1 if prob >= j else 0 for prob in y_prob]
                # 马修斯相关系数是在使用机器学习作为二进制（2类）的质量的度量的分类
                # print(y_test[:,i])
                # print(y_pred)
                acc.append(matthews_corrcoef(y_test[:, i], y_pred))
            # print(acc)
            acc = np.array(acc)
            index = np.where(acc == acc.max())
            accuracies.append(acc.max())
            best_threshold[i] = threshold[index[0][0]]
            acc = []

        # print("best thresholds", best_threshold)
        y_pred = np.array([[1 if y_probs[i, j] >= best_threshold[j] else 0 for j in range(y_test.shape[1])] for i in
                           range(len(y_test))])
        # print(y_pred)

        accuracy_for_n_classes = []
        hamming_loss_for_n_classes = []
        for i in range(nb_classes):
            # print('accuracy_score of class %d ' %i + str(accuracy_score(y_test[:, i],y_pred[:, i])))
            accuracy_for_n_classes.append(str(accuracy_score(y_test[:, i], y_pred[:, i])))
            hamming_loss_for_n_classes.append(str(hamming_loss(y_test[:, i], y_pred[:, i])))

        results_output_fp = dir_path + 'results/n_classes_accuracy_results_by_kfold.csv'
        with open(results_output_fp, 'a') as output_file:
            output_file.write(str(kfold) + ',' + 'accuracy' + ',' + ','.join(accuracy_for_n_classes) + '\n')
            output_file.write(str(kfold) + ',' + 'hamming_loss' + ',' + ','.join(hamming_loss_for_n_classes) + '\n')

        hamming_loss_score = hamming_loss(y_test, y_pred)
        print('hamming_loss_score: %.4f' % hamming_loss_score)
        jaccard_similarity = jaccard_similarity_score(y_test, y_pred)
        print('jaccard_similarity: %.4f' % jaccard_similarity)

        f1_score_micro = f1_score(y_test, y_pred, average='micro')
        print('f1_score_micro: %.4f' % f1_score_micro)

        accuracy = accuracy_score(y_test, y_pred)
        print('accuracy_score: %.4f' % accuracy)

        evaluate_scores = []
        evaluate_scores.append(hamming_loss_score)
        evaluate_scores.append(jaccard_similarity)
        evaluate_scores.append(auc_micro)
        evaluate_scores.append(average_precision_micro)
        evaluate_scores.append(f1_score_micro)
        evaluate_scores.append(accuracy)
        scores_output_fp = dir_path + 'results/k-fold_model_evaluation_indicators_results.csv'
        with open(scores_output_fp, 'a') as output_file:
            temp_values = [str(x) for x in evaluate_scores]
            output_file.write(str(kfold) + ',' + ','.join(temp_values) + '\n')
        auc_figure_fp = model_checkpoint_path + completed_time + '_auc.png'
        visualize_ROC_performance(y_probs, y_test, auc_figure_fp)
        pr_figure_fp = model_checkpoint_path + completed_time + '_pr.png'
        visualize_PR_performance(y_probs, y_test, pr_figure_fp)
        fold_k += 1

# Set image annotations
def annotate_style(rad):
    connectionstyle = "arc3,rad=%s" %rad
    arrowprops = dict(arrowstyle="->",
                      color="thistle",
                      linestyle=':',
                      shrinkA=5, shrinkB=5,
                      patchA=None,
                      patchB=None,
                      connectionstyle=connectionstyle,
                      )
    return arrowprops

# Discover the parameters and further optimize the prediction
def visualize_ACC_performance(dir_path,kfolds):

    # 读入准确率与汉明损失
    model_results_fp = dir_path + 'results/'
    file_names = [f for f in listdir(model_results_fp) if f.endswith('by_kfold.csv')]
    print('The amount of files: %d' % len(file_names))

    for file_name in file_names:
        model_results_acc = OrderedDict()
        model_results_acc_mean = np.zeros((nb_classes))
        model_results_loss = OrderedDict()
        model_results_loss_mean = np.zeros((nb_classes))
        file_path = model_results_fp + file_name
        opt = file_name.split('_')[0].strip()
        # if opt != 'adam':
        #     continue
        k_folds = set()

        model_total_results_acc = OrderedDict()
        model_total_results_loss = OrderedDict()
        with codecs.open(file_path, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                if temp[0] not in k_folds:
                    model_results_acc[temp[0]] = np.zeros((nb_classes))
                    model_results_loss[temp[0]] = np.zeros((nb_classes))
                    model_total_results_acc[temp[0]] = []
                    model_total_results_loss[temp[0]] = []
                    k_folds.add(temp[0])
                result_values = np.array([float(x) for x in temp[2:]])
                if temp[1] == 'accuracy':
                    model_results_acc[temp[0]] += result_values
                    model_total_results_acc[temp[0]].append([float(x) for x in temp[2:]])
                else:
                    model_results_loss[temp[0]] += result_values
                    model_total_results_loss[temp[0]].append([float(x) for x in temp[2:]])


        for key in model_results_acc.keys():
            model_results_acc_mean += (model_results_acc[key] / int(key))

        for key in model_results_loss.keys():
            model_results_loss_mean += (model_results_loss[key] / int(key))

        # 计算全局标准差
        model_results_acc_std = []
        for key in model_total_results_acc.keys():
            temp_nparray = np.array(model_total_results_acc[key])
            # print(key,temp_nparray.shape)
            model_results_acc_std.append(np.std(temp_nparray))
        model_results_loss_std = []
        for key in model_total_results_loss.keys():
            temp_nparray = np.array(model_total_results_loss[key])
            # print(key,temp_nparray.shape)
            model_results_loss_std.append(np.std(temp_nparray))

        # plot the training loss and accuracy
        font = {'family': 'Times New Roman',
                'weight': 'normal',
                'size': 10,
                }

        plt.close()
        plt.style.use('ggplot')
        plt.figure()
        ax = plt.subplot(111)
        N = len(kfolds) + 3
        for i in range(nb_classes):
            if i == 7:
                ax.plot(np.arange(3, N), model_results_acc[str(i + 3)] / (i + 3), color='olive', label="class %d" % (i+1), linewidth=0.7)
                continue
            ax.plot(np.arange(3, N), model_results_acc[str(i+3)] / (i+3), label="Class %d" % (i+1), linewidth=0.7)
        ax.plot(np.arange(3, N), model_results_acc_mean / nb_classes, 'dodgerblue',linestyle='--',label="mean", linewidth=1.7)
        ax.plot(np.arange(3, N), model_results_acc_std, 'tomato', linestyle='-.',
                label="Standard deviation", linewidth=1.7)
        # plt.title("Training Loss and Accuracy").set_fontname('Times New Roman')
        plt.xlabel("K-Fold", fontsize=10).set_fontname('Times New Roman')
        plt.ylabel("Accuracy", fontsize=10).set_fontname('Times New Roman')

        bbox = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
        arrowprops = annotate_style(-0.3)
        plt.annotate('(%s, %.4f)' % ('MAX: 5-fold', list(model_results_acc_mean / nb_classes)[2]),
                     xy=(kfolds[2], list(model_results_acc_mean / nb_classes)[2]), xycoords='data',  # figure pixels
                     xytext=(0.7 * kfolds[2], 0.5 * list(model_results_acc_mean / nb_classes)[2]), textcoords='data',
                     # offset points
                     bbox=bbox, arrowprops=arrowprops, size=10).set_fontname('Times New Roman')
        arrowprops = annotate_style(0.3)
        plt.annotate('(%s, %.4f)' % ('MAX: 5-fold', model_results_acc_std[2]),
                     xy=(kfolds[2], model_results_acc_std[2]), xycoords='data',  # figure pixels
                     xytext=(0.7 * kfolds[2], 0.4 * list(model_results_acc_mean / nb_classes)[2]), textcoords='data',
                     # offset points
                     bbox=bbox, arrowprops=arrowprops, size=10).set_fontname('Times New Roman')

        # arrowprops = annotate_style(-0.3)
        # plt.annotate('(%s, %.4f)' % ('8-fold', list(model_results_acc_mean / nb_classes)[-3]),
        #              xy=(kfolds[-3], list(model_results_acc_mean / nb_classes)[-3]), xycoords='data',  # figure pixels
        #              xytext=(0.75 * kfolds[-3], 0.35 * list(model_results_acc_mean / nb_classes)[0]),
        #              textcoords='data',
        #              # offset points
        #              bbox=bbox, arrowprops=arrowprops, size=10).set_fontname('Times New Roman')
        # arrowprops = annotate_style(0.3)
        # plt.annotate('(%s, %.4f)' % ('8-fold', model_results_acc_std[-3]),
        #              xy=(kfolds[-3], model_results_acc_std[-3]), xycoords='data',  # figure pixels
        #              xytext=(0.75 * kfolds[-3], 0.28 * list(model_results_acc_mean / nb_classes)[0]),
        #              textcoords='data',
        #              # offset points
        #              bbox=bbox, arrowprops=arrowprops, size=10).set_fontname('Times New Roman')


        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        # ax.set_position([box.x0, box.y0 + box.height * 0.1,
        #                  box.width, box.height * 0.9])
        ax.legend(prop=font, loc='center left', bbox_to_anchor=(1., 0.7))
        # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
        #           fancybox=True, shadow=True, ncol=5)
        image_path = dir_path + 'results/'+opt+'_k-fold_accuracy.pdf'
        plt.savefig(image_path,bbox_inches='tight')
        plt.show()

        plt.close()
        plt.style.use('ggplot')
        plt.figure()
        ax = plt.subplot(111)
        N = len(kfolds) + 3
        for i in range(nb_classes):
            if i == 7:
                ax.plot(np.arange(3, N), model_results_loss[str(i + 3)] / (i + 3), color='olive', label="class %d" % (i+1),
                         linewidth=0.7)
                continue
            ax.plot(np.arange(3, N), model_results_loss[str(i + 3)] / (i + 3), label="Class %d" % (i+1), linewidth=0.7)
        ax.plot(np.arange(3, N), model_results_loss_mean / nb_classes, color='dodgerblue', linestyle='--',
                 label="mean", linewidth=1.7)
        ax.plot(np.arange(3, N), model_results_loss_std, 'tomato', linestyle='-.',
                label="Standard deviation", linewidth=1.7)
        # plt.title("Training Loss and Accuracy").set_fontname('Times New Roman')
        plt.xlabel("K-Fold", fontsize=10).set_fontname('Times New Roman')
        plt.ylabel("Hamming loss", fontsize=10).set_fontname('Times New Roman')

        bbox = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
        arrowprops = annotate_style(0.3)
        plt.annotate('(%s, %.4f)' % ('MIN: 5-fold', list(model_results_loss_mean / nb_classes)[2]),
                     xy=(kfolds[2], list(model_results_loss_mean / nb_classes)[2]), xycoords='data',  # figure pixels
                     xytext=(1.3*kfolds[2], 1.7 * list(model_results_loss_mean / nb_classes)[2]), textcoords='data',
                     # offset points
                     bbox=bbox, arrowprops=arrowprops, size=10).set_fontname('Times New Roman')
        arrowprops = annotate_style(-0.3)
        plt.annotate('(%s, %.4f)' % ('MIN: 5-fold', model_results_loss_std[2]),
                     xy=(kfolds[2], model_results_loss_std[2]), xycoords='data',  # figure pixels
                     xytext=(1.3 * kfolds[2], 1.0 * list(model_results_loss_mean / nb_classes)[2]), textcoords='data',
                     # offset points
                     bbox=bbox, arrowprops=arrowprops, size=10).set_fontname('Times New Roman')

        # arrowprops = annotate_style(-0.3)
        # plt.annotate('(%s, %.4f)' % ('8-fold', list(model_results_loss_mean / nb_classes)[-3]),
        #              xy=(kfolds[-3], list(model_results_loss_mean / nb_classes)[-3]), xycoords='data',  # figure pixels
        #              xytext=(0.8 * kfolds[-3], 0.8 * list(model_results_loss_mean / nb_classes)[0]), textcoords='data',
        #              # offset points
        #              bbox=bbox, arrowprops=arrowprops, size=10).set_fontname('Times New Roman')
        # arrowprops = annotate_style(0.3)
        # plt.annotate('(%s, %.4f)' % ('8-fold', model_results_loss_std[-3]),
        #              xy=(kfolds[-3], model_results_loss_std[-3]), xycoords='data',  # figure pixels
        #              xytext=(0.8 * kfolds[-3], 0.7 * list(model_results_loss_mean / nb_classes)[0]), textcoords='data',
        #              # offset points
        #              bbox=bbox, arrowprops=arrowprops, size=10).set_fontname('Times New Roman')

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        # ax.set_position([box.x0, box.y0 + box.height * 0.1,
        #                  box.width, box.height * 0.9])
        ax.legend(prop=font, loc='center left', bbox_to_anchor=(1., 0.7))
        # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
        #           fancybox=True, shadow=True, ncol=5)
        image_path = dir_path + 'results/' + opt + '_k-fold_hanmming_loss.pdf'
        plt.savefig(image_path,bbox_inches='tight')
        plt.show()

# Visual Model Evaluation
def visualize_PR_performance(y_probs,y_test,image_path):

    precision = dict()
    recall = dict()
    average_precision = dict()

    # Compute micro-average ROC curve and ROC area
    precision_micro, recall_micro, _ = precision_recall_curve(y_test.ravel(), y_probs.ravel())  # 把多维的数组降为1维
    average_precision_micro = average_precision_score(y_test, y_probs, average="micro")

    precision['micro'] = precision_micro
    recall['micro'] = recall_micro
    average_precision['micro'] = average_precision_micro

    # Compute Precision-Recall and plot curve
    for i in range(nb_classes):
        precision[i], recall[i], _ = precision_recall_curve(y_test[:, i],
                                                            y_probs[:, i])
        # print(precision[i])
        average_precision[i] = average_precision_score(y_test[:, i], y_probs[:, i])
        # print(average_precision[i])

    # Plot Precision-Recall curve for each class
    font = {'family': 'Times New Roman',
            'weight': 'normal',
            'size': 9,
            }
    lw = 0.5
    colors = cycle(["lightpink", "g", "powderblue", 'b', 'darkorange', 'cornflowerblue', "peachpuff", "fuchsia"])
    plt.clf()
    plt.plot(recall['micro'], precision['micro'],
             label='micro-average Precision-recall curve (area = {0:0.2f})'
                   ''.format(average_precision['micro']),color='deeppink', linestyle=':', linewidth=1.0)
    for i, color in zip(range(nb_classes), colors):
        plt.plot(recall[i], precision[i], color=color, lw=lw,
                 label='Precision-recall curve of class {0} (area = {1:0.2f})'
                       ''.format(i, average_precision[i]))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall', fontsize=10).set_fontname('Times New Roman')
    plt.ylabel('Precision', fontsize=10).set_fontname('Times New Roman')
    # plt.title('Extension of Precision-Recall curve to multi-class')
    plt.legend(loc="upper right", prop=font)
    plt.savefig(image_path)  # 保存图片到本地
    # plt.show()

# Visual Model Evaluation
def visualize_ROC_performance(y_probs_list,y_test_list,kfold):


    plt.close('all')
    plt.style.use('ggplot')
    plt.figure(1)
    font = {'family': 'Times New Roman',
            'weight': 'normal',
            'size': 6,
            }
    for kf in range(kfold):
        y_test = y_test_list[kf]
        y_probs = y_probs_list[kf]
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(nb_classes):
            fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_probs[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_probs.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

        # First aggregate all false positive rates
        all_fpr = np.unique(np.concatenate([fpr[i] for i in range(nb_classes)]))
        # Then interpolate all ROC curves at this points
        mean_tpr = np.zeros_like(all_fpr)
        for i in range(nb_classes):
            mean_tpr += interp(all_fpr, fpr[i], tpr[i])

        # Finally average it and compute AUC
        mean_tpr /= nb_classes

        fpr['macro'] = all_fpr
        tpr['macro'] = mean_tpr
        roc_auc['macro'] = auc(fpr["macro"], tpr["macro"])

        # Plot all ROC curves
        plt.subplot(231+kf)
        lw = 0.5
        plt.plot(fpr["micro"], tpr["micro"],
                 label='Micro-average ROC curve (AUC = {0:0.2f})'
                       ''.format(roc_auc["micro"]),
                 color='deeppink', linestyle=':', linewidth=1.0)

        plt.plot(fpr["macro"], tpr["macro"],
                 label='Macro-average ROC curve (AUC = {0:0.2f})'
                       ''.format(roc_auc["macro"]),
                 color='navy', linestyle=':', linewidth=1.0)

        colors = cycle(["lightpink", "g", "powderblue", 'b', 'darkorange', 'cornflowerblue',"peachpuff", "fuchsia"])
        for i, color in zip(range(nb_classes), colors):
            plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                     label='ROC curve of class {0} (AUC = {1:0.2f})'
                           ''.format(i+1, roc_auc[i]))
        plt.plot([0, 1], [0, 1], 'k--', lw=lw)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False positive rate', fontsize=9).set_fontname('Times New Roman')
        plt.ylabel('True positive rate', fontsize=9).set_fontname('Times New Roman')
        plt.title('Fold %d' % kf, fontsize=8).set_fontname('Times New Roman')
        # plt.legend(loc="lower right", prop=font)
        plt.grid(True)

    plt.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0., prop=font)
    plt.tight_layout(pad=0.5)
    image_path = dir_path + 'results/' + str(kfold) + '-fold_ROC_curves.png'
    plt.savefig(image_path)  # 保存图片到本地
    image_path = dir_path + 'results/' + str(kfold) + '-fold_ROC_curves.pdf'
    plt.savefig(image_path)  # 保存图片到本地
    # plt.show()

# Load the model and predict the test set data
def load_model_for_prediction(dir_path,kfold):

    # K-fold cross 训练模型
    X, Y = prepare_dataset_for_train(dir_path)
    kf = KFold(n_splits=kfold)

    best_models_hdf5s = ['1551097132_weights_model.h5','1551097269_weights_model.h5','1551097407_weights_model.h5',
                         '1551097544_weights_model.h5','1551097682_weights_model.h5']

    # for best_models_hdf5 in best_models_hdf5s:
    fold_k = 0
    y_test_list = []
    y_probs_list = []
    for train_index, test_index in kf.split(X):
        _, X_test = X[train_index], X[test_index]
        _, y_test = Y[train_index], Y[test_index]

        if K.image_data_format() == 'channels_first':
            X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
            input_shape = (1, img_rows, img_cols)
        else:
            X_test = X_test.reshape(X_test.shape[0], img_rows, img_cols, 1)
            input_shape = (img_rows, img_cols, 1)

        # 将X_train, X_test的数据格式转为float32
        X_test = X_test.astype('float64')

        y_test_list.append(y_test)

        # 打印出相关信息
        print(X_test.shape)
        print(y_test.shape)
        model = build_model_architecture(input_shape)
        weights_path = dir_path + 'results/'+str(kfold)+'-fold/' + best_models_hdf5s[fold_k]
        model.load_weights(weights_path)
        y_probs = np.array(model.predict_proba(X_test))
        y_probs_list.append(y_probs)

        average_precision_micro = average_precision_score(y_test, y_probs, average="micro")
        average_precision_macro = average_precision_score(y_test, y_probs, average="macro")
        auc_micro = roc_auc_score(y_test, y_probs, average='micro')
        auc_macro = roc_auc_score(y_test, y_probs, average='macro')

        print('average_precision_micro: %.4f' %average_precision_micro)
        print('average_precision_macro: %.4f' % average_precision_macro)
        print('roc_auc_score_micro: %.4f' % auc_micro)
        print('roc_auc_score_macro: %.4f' % auc_macro)

        threshold = np.arange(0.1, 0.9, 0.1)
        acc = []
        accuracies = []
        best_threshold = np.zeros(y_probs.shape[1])
        for i in range(y_probs.shape[1]):
            y_prob = np.array(y_probs[:, i])
            for j in threshold:
                y_pred = [1 if prob >= j else 0 for prob in y_prob]
                # 马修斯相关系数是在使用机器学习作为二进制（2类）的质量的度量的分类
                # print(y_test[:,i])
                # print(y_pred)
                acc.append(matthews_corrcoef(y_test[:, i], y_pred))
            # print(acc)
            acc = np.array(acc)
            index = np.where(acc == acc.max())
            accuracies.append(acc.max())
            best_threshold[i] = threshold[index[0][0]]
            acc = []

        # print("best thresholds", best_threshold)
        y_pred = np.array([[1 if y_probs[i, j] >= best_threshold[j] else 0 for j in range(y_test.shape[1])] for i in
                           range(len(y_test))])

        print(y_test.shape, y_pred.shape)
        hamming_loss_score = hamming_loss(y_test, y_pred)
        print('hamming_loss_score: %.4f' % hamming_loss_score)
        jaccard_similarity = jaccard_similarity_score(y_test, y_pred)
        print('jaccard_similarity: %.4f' % jaccard_similarity)

        f1_score_micro = f1_score(y_test, y_pred, average='micro')
        print('f1_score_micro: %.4f' % f1_score_micro)
        f1_score_macro = f1_score(y_test, y_pred, average='macro')
        print('f1_score_macro: %.4f' % f1_score_macro)

        evaluate_scores= []
        evaluate_scores.append(hamming_loss_score)
        evaluate_scores.append(jaccard_similarity)
        evaluate_scores.append(auc_micro)
        evaluate_scores.append(auc_macro)
        evaluate_scores.append(average_precision_micro)
        evaluate_scores.append(average_precision_macro)
        evaluate_scores.append(f1_score_micro)
        evaluate_scores.append(f1_score_macro)
        scores_output_fp = dir_path + 'results/k-fold_model_evaluation_indicators_results_reload.csv'
        with open(scores_output_fp, 'a') as output_file:
            temp_values = [str(x) for x in evaluate_scores]
            output_file.write(str(kfold) + ',' + ','.join(temp_values) + '\n')
        fold_k += 1
    visualize_ROC_performance(y_probs_list, y_test_list, kfold)

# Visualize the results of different models
def visualize_models_compared_results(dir_path,kfold):

    # DF-CNN-SLC
    model_mean_results = OrderedDict()
    dl_model_results = []
    indicators = ['Hamming loss','Jaccard similarity','Micro AUC','Macro AUC','Micro Avg. precision',
                  'Macro Avg. precision','Micro F1 score','Macro F1 Score']
    dl_model_results_fp = dir_path + 'results/k-fold_model_evaluation_indicators_results_reload.csv'
    with codecs.open(dl_model_results_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if temp[0] != str(kfold):
                continue
            dl_model_results.append([float(x) for x in temp[1:]])
    sum_values = np.zeros((len(indicators)))
    if len(dl_model_results) == kfold:
        for model_result in dl_model_results:
            sum_values += np.array(model_result)
    model_mean_results['DP-CNN-MLMC'] = list(sum_values / kfold)

    # 读入对比模型的精确率和召回率
    compare_models_results_fp = dir_path + 'compare_models_results/'+str(kfold)+'-fold_compared_models_evaluation_indicators_results.csv'
    compare_models_results = []
    with codecs.open(compare_models_results_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            compare_models_results.append([float(x) for x in temp[1:]])
            if len(compare_models_results) == kfold:
                sum_values = np.zeros((len(indicators)))
                for model_result in compare_models_results:
                    sum_values += np.array(model_result)
                model_mean_results[temp[0]] = list(sum_values/kfold)
                compare_models_results = []
    print('model_mean_results: ' + str(len(model_mean_results.keys())))

    # Plot Precision-Recall curve for each model
    # font = {'family': 'Times New Roman',
    #         'weight': 'normal',
    #         'size': 10,
    #         }

    # # model_mean_results = {}
    # # for al in model_results.keys():
    # #     mean_indicators = []
    # #     for values in model_results[al]:
    # #         narray = numpy.array(values)
    # #         mean_indicators.append(narray.mean())
    # #     model_mean_results[al] = mean_indicators
    # for key in model_mean_results.keys():
    #     print(key,model_mean_results[key])
    #
    # N = len(indicators)
    # angles = np.linspace(0, 2 * np.pi, N, endpoint=False)  # 设置雷达图的角度，用于平分切开一个圆面
    # angles = np.concatenate((angles, [angles[0]]))  # 为了使雷达图一圈封闭起来
    # plt.close()
    # plt.style.use('ggplot')
    # fig = plt.figure(figsize=(7, 7))  # 设置画布大小
    # ax = fig.add_subplot(111, polar=True)  # 这里一定要设置为极坐标格式
    # sam = ['x-', '-.', '--', '-', ':', '-.']  # 样式
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
    #     ax.set_thetagrids(angles * 180 / np.pi, feature,fontsize=10,family='Times New Roman')  # 添加每个特征的标签
    #     ax.set_ylim(0.17, 0.71)  # 设置雷达图的范围
    #     ax.grid(linestyle='-', linewidth=0.5)  # 添加网格线
    #     lab.append(key)
    #     s += 1
    # plt.legend(lab, prop=font)
    # image_path = dir_path + 'results/'+str(kfold)+'-fold_comparison_of_prediction_results_of_different_models_total_indicators.png'
    # plt.savefig(image_path)  # 保存图片到本地
    # plt.show()  # 显示图形

    models = list(model_mean_results.keys())

    first_row = ['evaluation_indexes'] + models + ['bigger_smaller_better_sign']
    sign_values = ['smaller is better','bigger is better','bigger is better','bigger is better',
                    'bigger is better','bigger is better','bigger is better','bigger is better']

    string_IO = ','.join(first_row) + '\n'
    for t in range(len(indicators)):
        row_t = indicators[t] + ',' + str(model_mean_results[models[0]][t]) + ',' + str(model_mean_results[models[1]][t]) + ',' + \
                str(model_mean_results[models[2]][t]) + ',' + str(model_mean_results[models[3]][t]) + ',' + \
                str(model_mean_results[models[4]][t]) + ',' + str(model_mean_results[models[5]][t]) + ',' + \
                sign_values[t] + '\n'
        string_IO += row_t

    models_color = OrderedDict([
        (models[0], "#FF6666"),
        (models[1], "#99CC33"),
        (models[2], "#33CC99"),
        (models[3], "#FFFF00"),
        (models[4], "#66CCFF"),
        (models[5], "#0099FF"),
    ])

    sign_color = {
        "bigger is better": "#FF99CC",
        "smaller is better": "#CCCCCC",
    }

    df = pd.read_csv(StringIO(string_IO),
                     skiprows=0,
                     skipinitialspace=True,
                     engine='python')
    df_output = df.stack()
    df_output.to_csv(dir_path + 'results/' + str(
        kfold) + '-fold_comparison_of_prediction_results_of_different_models_total_indicators.csv')

    width = 800
    height = 800
    inner_radius = 90
    outer_radius = 300 - 10

    minr = sqrt(log(.1 * 1E4))
    maxr = sqrt(log(.8 * 1E4))
    print(minr, maxr)
    a = (outer_radius - inner_radius) / (maxr - minr)
    b = inner_radius - a * minr
    print(a, b)

    def rad(mic):
        return a * np.sqrt(np.log(mic * 1E4)) + b

    big_angle = 2.0 * np.pi / (len(df) + 1)
    small_angle = big_angle / (len(df) - 1.5)
    print(big_angle,small_angle)

    p = figure(plot_width=width, plot_height=height, title="",
               x_axis_type=None, y_axis_type=None,
               x_range=(-420, 420), y_range=(-420, 420),
               min_border=0, outline_line_color="#000000",
               background_fill_color="white") #f0e1d2

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None

    # annular wedges
    angles = np.pi / 2 - big_angle / 2 - df.index.to_series() * big_angle #计算角度
    colors = [sign_color[bigger_smaller_better_sign] for bigger_smaller_better_sign in df.bigger_smaller_better_sign]# 配置颜色
    p.annular_wedge(
        0, 0, inner_radius, outer_radius, -big_angle + angles, angles, color=colors,
    )

    # small wedges
    p.annular_wedge(0, 0, inner_radius, rad(df['DP-CNN-MLMC']),
                    -big_angle + angles + 5.5 * small_angle, -big_angle + angles + 6 * small_angle,
                    color=models_color[models[0]])
    p.annular_wedge(0, 0, inner_radius, rad(df.MK),
                    -big_angle + angles + 4.5 * small_angle, -big_angle + angles + 5 * small_angle,
                    color=models_color[models[1]])
    p.annular_wedge(0, 0, inner_radius, rad(df.NB),
                    -big_angle + angles + 3.5 * small_angle, -big_angle + angles + 4 * small_angle,
                    color=models_color[models[2]])
    p.annular_wedge(0, 0, inner_radius, rad(df.LR),
                    -big_angle + angles + 2.5 * small_angle, -big_angle + angles + 3 * small_angle,
                    color=models_color[models[3]])
    p.annular_wedge(0, 0, inner_radius, rad(df.SVM),
                    -big_angle + angles + 1.5 * small_angle, -big_angle + angles + 2 * small_angle,
                    color=models_color[models[4]])
    p.annular_wedge(0, 0, inner_radius, rad(df.RF),
                    -big_angle + angles + 0.5 * small_angle, -big_angle + angles + 1 * small_angle,
                    color=models_color[models[5]])

    # 绘制大圆和标签
    # circular axes and lables
    # labels = np.power(10.0, np.arange(-3, 4))
    labels = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    print(labels)
    radii = a * np.sqrt(np.log(labels * 1E4)) + b
    # print(np.sqrt(np.log(labels * 1E4)))
    # radii = np.array([90., 116.17467948, 142.17467948, 168.17467948, 194.17467948, 230.17467948,
    #                   256.17467948, 290.])
    print(radii)
    p.circle(0, 0, radius=radii, fill_color=None, line_color="#CCFF66")
    p.text(0, radii[:], [str(r) for r in labels[:]],text_font='Times New Roman',
           text_font_size="10pt", text_align="center", text_baseline="middle")

    # radial axes
    p.annular_wedge(0, 0, inner_radius - 10, outer_radius + 10,
                    -big_angle + angles, -big_angle + angles, color="black")

    # bacteria labels
    xr = radii[-1] * np.cos(np.array(-big_angle / 2 + angles))
    yr = radii[-1] * np.sin(np.array(-big_angle / 2 + angles))
    label_angle = np.array(-big_angle / 2 + angles)
    label_angle[label_angle < -np.pi / 2] += np.pi  # easier to read labels on the left side
    p.text(xr, yr, df.evaluation_indexes, angle=label_angle,text_font='Times New Roman',
           text_font_size="10pt", text_align="center", text_baseline="middle")

    # OK, these hand drawn legends are pretty clunky, will be improved in future release
    p.circle([-40, -40], [-370, -390], color=list(sign_color.values()), radius=5)
    p.text([-30, -30], [-370, -390], text=["The " + gr for gr in sign_color.keys()],text_font='Times New Roman',
           text_font_size="10pt", text_align="left", text_baseline="middle")

    p.rect([-40, -40, -40, -40, -40, -40], [36,18, 0, -18,-36,-54], width=20, height=9,
           color=list(models_color.values()))
    p.text([-15, -15, -15, -15, -15, -15], [36,18, 0, -18,-36,-54], text=list(models_color),text_font='Times New Roman',
           text_font_size="8pt", text_align="left", text_baseline="middle")
    show(p)
    image_path_html = dir_path + 'results/' + str(
        kfold) + '-fold_comparison_of_prediction_results_of_different_models_total_indicators.html'
    output_file(image_path_html)
    save(obj=p,filename=image_path_html)  # 保存图片到本地

    # image_path = dir_path + 'results/' + str(
    #     kfold) + '-fold_comparison_of_prediction_results_of_different_models_total_indicators.svg'
    # from bokeh.io import export_svgs
    # p.output_backend = "svg"
    # export_svgs(p, image_path)


if __name__ == '__main__':

    dir_path = ''
    # kfolds = [3,4,5,6,7,8,9,10]
    # for kfold in kfolds:
    #     train_model(dir_path, kfold)
    # visualize_ACC_performance(dir_path, kfolds)
    # visualize_PR_performance(dir_path)
    kfold = 5
    # load_model_for_prediction(dir_path, kfold)
    visualize_models_compared_results(dir_path,kfold)

