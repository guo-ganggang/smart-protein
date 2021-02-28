#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 12/1/2021 9:55 PM
# @Author  : ggguo
# @Site    : SYM PROJECT
# @File    : train_compare_models_by_xgboost.py
# @Software: SYM application case
"""
.py:
"""

from __future__ import division

import pandas as pd
import numpy as np
import warnings
from sklearn.preprocessing import scale
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from xgboost.sklearn import XGBClassifier
from sklearn.preprocessing import normalize
import lightgbm as lgb
from collections import OrderedDict
import codecs

import sys

reload(sys)
sys.setdefaultencoding('utf8')


def build_model(model_name):
    model = None
    if model_name == 'xgb':
        model = XGBClassifier(random_state=2020, n_estimators=1000, learning_rate=0.02)
    elif model_name == 'svm':
        model = SVC(random_state=2020, tol=1e-6, probability=True)
    elif model_name == 'dt':
        model = DecisionTreeClassifier(random_state=2020)
    elif model_name == 'lr':
        model = LogisticRegression(random_state=2020, C=1, penalty='l1', tol=1e-6)
    elif model_name == 'rf':
        model = RandomForestClassifier(n_estimators=60, random_state=2020)
    elif model_name == 'gbdt':
        model = GradientBoostingClassifier(random_state=2020)
    elif model_name == 'lgbm':
        model = lgb.LGBMClassifier()
    else:
        print('Please input the correct model name.')

    return model


def train_models_by_skl_xgb_lgb(dir_path):
    # 读取数据集
    data_all = pd.read_csv(dir_path + 'train_data_protein_symbol_features.csv')

    # 划分为4折交叉验证数据集
    train_data = data_all.drop(columns=['gene_name'])
    df_y = train_data['label']
    df_X = train_data.drop(columns=['label'])
    df_X = scale(df_X, axis=0)  # Z-Score: (X-mean)/std
    # 构建模型

    lr = LogisticRegression(random_state=2018, tol=1e-6)  # 逻辑回归模型
    dt = DecisionTreeClassifier(random_state=2018)  # 决策树模型
    svm = SVC(probability=True, random_state=2018, tol=1e-6)  # SVM模型
    rf = RandomForestClassifier(n_estimators=100, random_state=2018)  # 随机森林
    gbdt = GradientBoostingClassifier(random_state=2018)  # CBDT
    xgb = XGBClassifier(random_state=2018)  # Xgbc
    lgbm = lgb.LGBMClassifier(random_state=2018)  # lgb

    def muti_score(model):
        warnings.filterwarnings('ignore')
        accuracy = cross_val_score(model, df_X, df_y, scoring='accuracy', cv=4)
        precision = cross_val_score(model, df_X, df_y, scoring='precision', cv=4)
        recall = cross_val_score(model, df_X, df_y, scoring='recall', cv=4)
        f1_score = cross_val_score(model, df_X, df_y, scoring='f1', cv=4)
        auc = cross_val_score(model, df_X, df_y, scoring='roc_auc', cv=4)
        print("Accuracy:", accuracy.mean())
        print("Precision:", precision.mean())
        print("Recall:", recall.mean())
        print("F1_score:", f1_score.mean())
        print("AUC:", auc.mean())

    model_name = ["lr", "dt", "svm", "rf", "gbdt", "xgb", "lgbm"]
    for name in model_name:
        model = eval(name)
        print(name)
        muti_score(model)


def evaluation_models(dir_path,kfold):
    # 读取数据集
    data_all = pd.read_csv(dir_path + 'compare_models_results/train_data_protein_symbol_features.csv')

    train_data = data_all.drop(columns=['gene_name'])
    data_zero = train_data[train_data['label'] == 0]
    data_one = train_data[train_data['label'] == 1]

    kf = KFold(n_splits=kfold, shuffle=True, random_state=2)
    folds_zero = kf.split(data_zero)
    folds_one = kf.split(data_one)

    model_names = ["lr", "dt", "svm", "rf", "gbdt", "xgb", "lgbm"]
    indicators = ["accuracy", "precision", "recall", "f1_score"]
    model_name_indicators = OrderedDict()
    for model_name in model_names:
        model_name_indicators[model_name] = OrderedDict()
        for indicator in indicators:
            model_name_indicators[model_name][indicator] = []

    for i in range(kfold):

        fz = next(folds_zero, None)
        fo = next(folds_one, None)
        (train_data_0, test_data_0) = data_zero.iloc[fz[0]], data_zero.iloc[fz[1]]
        (train_data_1, test_data_1) = data_one.iloc[fo[0]], data_one.iloc[fo[1]]

        train_data_1 = pd.concat([train_data_1] * 3, ignore_index=True)
        test_data_1 = pd.concat([test_data_1] * 3, ignore_index=True)

        test_data = test_data_0.append(test_data_1)
        train_data = train_data_0.append(train_data_1)

        X_train, X_test = train_data.drop(columns=['label']), test_data.drop(columns=['label'])
        Y_train, Y_test = train_data['label'], test_data['label']

        X_train_normalized = normalize(X_train, norm='l2')
        X_test_normalized = normalize(X_test, norm='l2')

        acc_list = []
        pre_list = []
        rec_list = []
        f1_list = []
        print('model_name: ' + ', '.join(model_names))
        for model_name in model_names:
            model = build_model(model_name)
            model.fit(X_train_normalized, Y_train)
            y_pred_class = model.predict(X_test_normalized)
            accuracy = metrics.accuracy_score(Y_test, y_pred_class)
            acc_list.append(accuracy)
            model_name_indicators[model_name]['accuracy'].append(accuracy)
            precision = metrics.precision_score(Y_test, y_pred_class)
            pre_list.append(precision)
            model_name_indicators[model_name]['precision'].append(precision)
            recall = metrics.recall_score(Y_test, y_pred_class)
            rec_list.append(recall)
            model_name_indicators[model_name]['recall'].append(recall)
            f1_score = metrics.f1_score(Y_test, y_pred_class)
            f1_list.append(f1_score)
            model_name_indicators[model_name]['f1_score'].append(f1_score)

        print('acc_list: ' + ', '.join([str(x) for x in acc_list]))
        print('pre_list: ' + ', '.join([str(x) for x in pre_list]))
        print('rec_list: ' + ', '.join([str(x) for x in rec_list]))
        print('f1_list: ' + ', '.join([str(x) for x in f1_list]))

    for model_name in model_name_indicators.keys():
        k_fold_results_fp = dir_path + 'compare_models_results/skl_xgb_lgb/' + model_name.upper() + '_Evaluate_%d' %kfold +'folds_total_results.csv'
        with codecs.open(k_fold_results_fp, "a", "utf-8") as oFile_handler:
            acc_list = model_name_indicators[model_name]['accuracy']
            oFile_handler.write('accuracy' + '\t' + '\t'.join([str(x) for x in acc_list]) + '\n')
            pre_list = model_name_indicators[model_name]['precision']
            oFile_handler.write('precision' + '\t' + '\t'.join([str(x) for x in pre_list]) + '\n')
            rec_list = model_name_indicators[model_name]['recall']
            oFile_handler.write('recall' + '\t' + '\t'.join([str(x) for x in rec_list]) + '\n')
            f1_list = model_name_indicators[model_name]['f1_score']
            oFile_handler.write('f1_score' + '\t' + '\t'.join([str(x) for x in f1_list]) + '\n')



if __name__ == "__main__":

    dir_path = "F:/new_disease/4_phase1_training/"
    # dir_path = "F:/4_phase1_training/"
    # train_models_by_skl_xgb_lgb(dir_path)
    kfold = 9
    evaluation_models(dir_path, kfold)
