#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 6/1/2019 2:50 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : train_compare_models_by_graphlab.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############


from __future__ import division
import turicreate as gl
from turicreate import SFrame
import codecs
import sys
import os
import pandas as pd
import os
import codecs
from os import listdir
from itertools import islice
from collections import OrderedDict
import datetime
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from functools import reduce
from sklearn.model_selection import KFold

# pip install --upgrade --no-cache-dir https://get.graphlab.com/GraphLab-Create/2.1/rjxyggg@csu.edu.cn/
# 9580-1D65-E031-BC57-527E-E66B-1BE6-F70D/GraphLab-Create-License.tar.gz

def prepare_train_data(dir_path):

    feature_input = dir_path + 'data_embedding/'
    file_names = [f for f in listdir(feature_input) if f.endswith('csv')]

    # train_data = gl.SFrame()
    feature_names = []
    raw_feature_data_out_path = dir_path + 'compare_models_results/train_data_protein_symbol_features.csv'

    if not os.path.exists(raw_feature_data_out_path):
        # 读入并合并feature
        train_data_features_labels = pd.DataFrame()
        for i in range(len(file_names)):
            file_name = file_names[i]
            print(file_name)
            feature_name = file_name.strip().split('_')[3].strip()
            feature_names.append(feature_name)
            train_data_path = feature_input + file_name
            train_data_feature = pd.read_csv(train_data_path, sep="$", header=0)

            # 分开
            train_data_feature_values = pd.DataFrame()
            flag = 0
            train_data_feature_split = train_data_feature['feature_values'].str.split(',').str
            for feature_value in train_data_feature_split:
                temp_feature_name = feature_name + '_' + str(flag)
                train_data_feature_values[temp_feature_name] = feature_value
                flag += 1

            # 降维
            train_data_feature_values_np = train_data_feature_values.values
            print(train_data_feature_values_np.shape)
            raw_features_pca = PCA(n_components=50).fit_transform(train_data_feature_values_np.astype('float'))
            raw_features_pca_tsne = TSNE(n_components=3, init='pca').fit_transform(raw_features_pca)  # .as_matrix()
            print(raw_features_pca_tsne.shape)

            train_data_feature_values_tsne = pd.DataFrame()
            train_data_feature_values_tsne['gene_name'] = train_data_feature['gene_name']

            new_flag = 0
            for feature_value in pd.DataFrame(raw_features_pca_tsne):
                temp_feature_name = feature_name + '_' + str(new_flag)
                train_data_feature_values_tsne[temp_feature_name] = pd.DataFrame(raw_features_pca_tsne)[feature_value]
                new_flag += 1
                # print(train_data_feature_values_tsne.head())
            print(train_data_feature_values_tsne.keys())

            if len(train_data_features_labels.keys()) == 0:
                train_data_features_labels = train_data_feature_values_tsne
            else:
                train_data_features_labels = pd.merge(train_data_features_labels, train_data_feature_values_tsne,
                                                      on=["gene_name", "gene_name"])
        # print(train_data_features.head)

        # 读入并合并label
        feature_input = dir_path + 'data_label/'
        file_names = [f for f in listdir(feature_input) if f.endswith('csv')]
        classification_protein_symbol_label = {}
        for file_name in file_names:
            print(file_name)
            key_name = file_name.strip().split('_')[3].strip()
            raw_data_fp = feature_input + file_name
            with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
                for line in islice(input_file.readlines(), 0, None):
                    temp = line.strip()
                    if temp == '':
                        continue
                    if key_name == 'positive':
                        classification_protein_symbol_label[temp] = 1
                    elif key_name == 'negative':
                        classification_protein_symbol_label[temp] = 0
                    else:
                        print('Error..................')
        train_data_labels = pd.DataFrame(list(classification_protein_symbol_label.items()), columns=['gene_name', 'label'])
        # print(train_data_labels.head)
        train_data_features_labels = pd.merge(train_data_features_labels, train_data_labels,
                                              on=["gene_name", "gene_name"])
        # train_data_features_labels.to_csv(raw_feature_data_out_path)

        # print(train_data_features.head)
        print(train_data_features_labels.shape)
        print('--------------------------------------')
        train_data = gl.SFrame(data=train_data_features_labels)
        print(train_data.shape)
        train_data.save(raw_feature_data_out_path,format='csv')

    else:
        train_data = gl.SFrame.read_csv(raw_feature_data_out_path)
        print(train_data.shape)
        for i in range(len(file_names)):
            file_name = file_names[i]
            feature_name = file_name.strip().split('_')[3].strip()
            if feature_name in feature_names:
                continue
            feature_names.append(feature_name)

    # 准备feature 名称
    print(feature_names)
    features_combination = OrderedDict()
    temp_feature_combination = []
    for i in range(len(feature_names)):
        for j in range(3):
            temp_feature_combination.append(feature_names[i]+'_'+str(j))
        features_combination[i] = temp_feature_combination
        print('temp_feature_combination: ' + str(len(temp_feature_combination)))

    return train_data, features_combination


class LeadsGenerationModelValidation:
    def __init__(self,train_data_all,feature,target,validation_outfile,key,k_folds,algorithm):
        self.train_data_all = train_data_all
        self.validation_outfile =  validation_outfile
        self.feature = feature
        self.target = target
        self.key = key
        self.k_folds = k_folds
        self.algorithm = algorithm

    def validate(self):
        acc_list = []
        pre_list = []
        rec_list = []
        f1_list = []

        with codecs.open(self.validation_outfile,"a+","utf-8") as output_handler:
            output_handler.write("Features: " + " ".join(self.feature) + "\n")
            iteration = 0

            data_zero = pd.DataFrame(self.train_data_all[self.train_data_all[self.target] == 0])
            data_one = pd.DataFrame(self.train_data_all[self.train_data_all[self.target] == 1])

            kf = KFold(n_splits=self.k_folds, shuffle = True, random_state = 2)
            folds_zero = kf.split(data_zero)
            folds_one = kf.split(data_one)

            for i in range(self.k_folds):
                iteration += 1
                fz = next(folds_zero, None)
                fo = next(folds_one, None)
                (train_data_0, test_data_0) = data_zero.iloc[fz[0]],data_zero.iloc[fz[1]]
                (train_data_1, test_data_1) = data_one.iloc[fo[0]],data_one.iloc[fo[1]]

                train_data_1 = pd.concat([train_data_1] * 3, ignore_index=True)
                test_data_1 = pd.concat([test_data_1] * 3, ignore_index=True)

                test_data = SFrame(data=test_data_0.append(test_data_1))
                train_data = SFrame(data=train_data_0.append(train_data_1))

                test_data = test_data.dropna()
                train_data = train_data.dropna()

                print("length of train_data:    %s and length of test_data: %s" % (len(train_data),len(test_data)))

                # Feature rescaling is an important pre-processing step that ensures that all features are on the same scale.
                # An l2-norm rescaling is performed to make sure that all features are of the same norm.
                if 'LR' == self.algorithm:
                    model = gl.logistic_classifier.create(train_data, target = self.target,features = self.feature, \
                                                              class_weights='auto')
                elif 'DT' == self.algorithm:
                    model = gl.decision_tree_classifier.create(train_data, target = self.target,features = self.feature, \
                                                              class_weights='auto')
                elif 'SVM' == self.algorithm:
                    model = gl.svm_classifier.create(train_data, target=self.target, features=self.feature, \
                                                              class_weights='auto')
                elif 'RF' == self.algorithm:
                    model = gl.random_forest_classifier.create(train_data, target = self.target,features = self.feature, \
                                                              class_weights='auto')
                elif 'BT' == self.algorithm:
                    model = gl.boosted_trees_classifier.create(train_data, target = self.target,features = self.feature, \
                                                              class_weights='auto')
                else:
                    print('Please input the correct algorithm!')

                classifier = model.classify(test_data)
                print(classifier.column_names())
                if 'SVM' == self.algorithm:
                    test_data["classify_class"] =  classifier["class"]
                else:
                    test_data["classify_probability"] = classifier["probability"]
                    test_data["classify_class"] = classifier["class"]

                results = model.evaluate(test_data)
                acc_list.append(results['accuracy'])
                pre_list.append(results['precision'])
                rec_list.append(results['recall'])
                f1_list.append(results['f1_score'])

                output_handler.write("Accuracy         :%s" % results['accuracy'] + "\n")
                output_handler.write("precision         :%s" % results['precision'] + "\n")
                output_handler.write("recall         :%s" % results['recall'] + "\n")
                output_handler.write("f1_score         :%s" % results['f1_score'] + "\n")
                output_handler.write("Confusion Matrix :%s" % results['confusion_matrix'] + "\n")

            total_acc = reduce(lambda x, y: x + y, acc_list)
            total_pre = reduce(lambda x, y: x + y, pre_list)
            total_rec = reduce(lambda x, y: x + y, rec_list)
            total_f1 = reduce(lambda x, y: x + y, f1_list)

            generate_date = str(datetime.date.today())
            save_dir = dir_path + 'compare_models_results/' + str(self.k_folds) + '-fold/'

            oFile = save_dir + self.algorithm + '_Evaluate_'+str(k_folds)+'folds_' + generate_date +'.txt'
            with codecs.open(oFile, "a", "utf-8") as oFile_handler:
                oFile_handler.write(str(key) + '\t' + str(float(total_acc) / iteration) \
                                    + '\t' + str(float(total_pre) / iteration) + '\t' + str(float(total_rec) / iteration) \
                                    + '\t' + str(float(total_f1) / iteration) \
                                    + '\n')
            if key == 4:
                oFile = save_dir + self.algorithm + '_Evaluate_' + str(k_folds) + 'folds_total_results_' + generate_date + '.txt'
                with codecs.open(oFile, "a", "utf-8") as oFile_handler:
                    oFile_handler.write('accuracy' + '\t' + '\t'.join([str(x) for x in acc_list])+ '\n')
                    oFile_handler.write('precision' + '\t' + '\t'.join([str(x) for x in pre_list]) + '\n')
                    oFile_handler.write('recall' + '\t' + '\t'.join([str(x) for x in rec_list]) + '\n')
                    oFile_handler.write('f1_score' + '\t' + '\t'.join([str(x) for x in f1_list]) + '\n')


if __name__ == "__main__":


    dir_path = ''
    train_data_feature_label, features_combination = prepare_train_data(dir_path)
    print(train_data_feature_label.shape)
    generate_date = str(datetime.date.today())
    algorithms = ['LR','DT','SVM','RF','BT']
    for algorithm in algorithms:
        for key, feature in features_combination.items():
            print(feature)
            k_folds = 9
            save_dir = dir_path + 'compare_models_results/' + str(k_folds) + '-fold/' + 'intermediate_results/'
            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)
            output_file = save_dir + algorithm + '_validatation_'+str(k_folds)+'folds_'+generate_date+'_%s.txt' % key
            model = LeadsGenerationModelValidation(train_data_feature_label, feature, 'label', output_file, key,k_folds,algorithm)
            model.validate()