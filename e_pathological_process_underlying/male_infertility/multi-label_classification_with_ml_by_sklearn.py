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

import codecs
from itertools import islice
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle

from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import hamming_loss
from sklearn.metrics import classification_report
from sklearn.metrics import jaccard_similarity_score
from sklearn.metrics import precision_recall_fscore_support
from sklearn.utils.fixes import signature
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from skmultilearn.adapt import MLkNN
from skmultilearn.problem_transform import LabelPowerset
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from skmultilearn.problem_transform import ClassifierChain
from sklearn import svm
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import f1_score

from xgboost import XGBClassifier
# from sklearn.multiclass import OneVsRestClassifier # multiclass分类
from sklearn.multioutput import MultiOutputClassifier # multilabel分类
from sklearn.pipeline import Pipeline

import warnings
warnings.filterwarnings('ignore') # "error", "ignore", "always", "default", "module" or "once"

# import sys
#
# reload(sys)
# sys.setdefaultencoding('utf8')


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
    file_for_feature = embedding_data_dir + 'train_data_protein_symbol_features'
    protein_symbol_features = {}
    with codecs.open(file_for_feature, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split(',')
            if len(temp) < 2:
                continue
            if temp[0].strip() not in protein_symbol_multi_labeled:
                continue
            protein_symbol_features[temp[0].strip()] = temp[1:]
    print('protein_symbol_features: ' + str(len(protein_symbol_features.keys())))

    train_test_dataset_x = []
    train_test_dataset_y = []
    for protein_symbol in protein_symbol_features.keys():
        feature_float = [float(x) for x in protein_symbol_features[protein_symbol]]
        train_test_dataset_x.append(feature_float)
        train_test_dataset_y.append(protein_symbol_multi_labels[protein_symbol])

    return train_test_dataset_x,train_test_dataset_y

# Visual model output
def visualize_re_performance(Y_test,y_pred):
    precision = dict()
    recall = dict()
    average_precision = dict()

    # Compute micro-average ROC curve and ROC area
    precision_micro, recall_micro, _ = precision_recall_curve(Y_test.ravel(),y_pred.ravel()) # 把多维的数组降为1维
    average_precision_micro = average_precision_score(Y_test, y_pred, average="micro")

    precision['micro'] = precision_micro
    recall['micro'] = recall_micro
    average_precision['micro'] = average_precision_micro

    print("precision : "+str(precision_micro)+",recall :"+str(recall_micro)+",average precision :"+str(average_precision_micro))

    step_kwargs = ({'step': 'post'}
                   if 'step' in signature(plt.fill_between).parameters
                   else {})

    plt.figure()
    plt.step(recall['micro'], precision['micro'], color='b', alpha=0.2,
             where='post')
    plt.fill_between(recall["micro"], precision["micro"], alpha=0.2, color='b',
                     **step_kwargs)

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title(
        'Average precision score, micro-averaged over all classes: AP={0:0.2f}'
            .format(average_precision["micro"]))
    plt.show()

    # Compute Precision-Recall and plot curve
    n_classes = Y_test.shape[1]

    for i in range(n_classes):
        precision[i], recall[i], _ = precision_recall_curve(Y_test[:, i],
                                                            y_pred[:, i])
        # print(precision[i])
        average_precision[i] = average_precision_score(Y_test[:, i], y_pred[:, i])
        # print(average_precision[i])

    # Plot Precision-Recall curve for each class
    plt.clf()
    plt.plot(recall_micro, precision_micro,
             label='micro-average Precision-recall curve (area = {0:0.2f})'
                   ''.format(average_precision_micro))
    for i in range(n_classes):
        plt.plot(recall[i], precision[i],
                 label='Precision-recall curve of class {0} (area = {1:0.2f})'
                       ''.format(i, average_precision[i]))

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Extension of Precision-Recall curve to multi-class')
    plt.legend(loc="upper right")
    plt.show()

    # setup plot details
    colors = cycle(['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal','pink','green','purple'])

    plt.figure(figsize=(7, 8))
    f_scores = np.linspace(0.2, 0.8, num=8)
    lines = []
    labels = []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))

    lines.append(l)
    labels.append('iso-f1 curves')
    l, = plt.plot(recall["micro"], precision["micro"], color='gold', lw=2)
    lines.append(l)
    labels.append('micro-average Precision-recall (area = {0:0.2f})'
                  ''.format(average_precision["micro"]))

    for i, color in zip(range(n_classes), colors):
        l, = plt.plot(recall[i], precision[i], color=color, lw=2)
        lines.append(l)
        labels.append('Precision-recall for class {0} (area = {1:0.2f})'
                      ''.format(i, average_precision[i]))

    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Extension of Precision-Recall curve to multi-class')
    plt.legend(lines, labels, loc=(0.25, -.48), prop=dict(size=9))

    plt.show()

# training model
def train_model(dir_path,X_train, X_test, Y_train, y_test,algorithm):

    # Algorithm Adaptation approaches
    # a multi-label adapted kNN classifier with bayesian prior corrections
    if algorithm == 'MK':
        parameters = {'k': range(1, 15), 's': [0.5, 0.7, 1.0]}
        score = 'f1_macro' #accuracy,f1_macro
        clf = GridSearchCV(MLkNN(), parameters, scoring=score)
        clf.fit(X_train, Y_train)
        # print(clf.best_params_, clf.best_score_)
        y_probs = clf.predict_proba(X_test).toarray()  # scipy.sparse.lil.lil_matrix
    # Problem Transformation approaches
    # treats each label combination as a separate class with one multi-class classification problem
    elif algorithm == 'NB':
        clf = LabelPowerset(GaussianNB())
        clf.fit(X_train, Y_train)
        y_probs = clf.predict_proba(X_test).toarray() # scipy.sparse.lil.lil_matrix
    elif algorithm == 'LR':
        clf = LabelPowerset(LogisticRegression())
        clf.fit(X_train, Y_train)
        y_probs = clf.predict_proba(X_test).toarray() # scipy.sparse.lil.lil_matrix
    elif algorithm == 'SVM':
        clf = LabelPowerset(SVC(probability=True))
        clf.fit(X_train, Y_train)
        y_probs = clf.predict_proba(X_test).toarray()  # scipy.sparse.lil.lil_matrix
    elif algorithm == 'RF':
        parameters = [
            {
                'classifier': [RandomForestClassifier()],
                'classifier__criterion': ['gini', 'entropy'],
                'classifier__n_estimators': [30,50,100],
            }
        ]
        clf = GridSearchCV(LabelPowerset(), parameters, scoring='accuracy')
        clf.fit(X_train, Y_train)
        # print (clf.best_params_, clf.best_score_)
        y_probs = clf.predict_proba(X_test).toarray()  # scipy.sparse.lil.lil_matrix
    elif algorithm == 'XGB':
        # clf = OneVsRestClassifier(XGBClassifier())
        # clf.fit(X_train, Y_train)
        # y_probs = clf.predict(X_test)
        xgb_estimator = XGBClassifier(objective='binary:logistic')
        clf = MultiOutputClassifier(xgb_estimator)
        clf.fit(X_train, Y_train)
        # y_probs = clf.predict(X_test)
        y_probs = clf.predict_proba(X_test) # list of preds per class
        y_probs = np.array(y_probs)[:, :, 1].T # take the positive class
        print('Accuracy on test data: {:.1f}%'.format(accuracy_score(y_test, clf.predict(X_test)) * 100))
    else:
        print('Please input the correct algorithm name!')
    # multilabel classification metrics
    print(y_test.shape,y_probs.shape)
    average_precision_micro = average_precision_score(y_test, y_probs, average="micro")
    average_precision_macro = average_precision_score(y_test, y_probs, average="macro")
    auc_micro = roc_auc_score(y_test, y_probs, average='micro')
    auc_macro = roc_auc_score(y_test, y_probs, average='macro')

    print('average_precision_micro: %.4f' % average_precision_micro)
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
    hamming_loss_score = hamming_loss(y_test, y_pred)
    print('hamming_loss_score: %.4f' % hamming_loss_score)
    # 值越大模型越好
    jaccard_similarity = jaccard_similarity_score(y_test, y_pred)
    print('jaccard_similarity: %.4f' % jaccard_similarity)

    f1_score_micro = f1_score(y_test, y_pred, average='micro')
    print('f1_score_micro: %.4f' % f1_score_micro)
    f1_score_macro = f1_score(y_test, y_pred, average='macro')
    print('f1_score_macro: %.4f' % f1_score_macro)

    evaluate_scores = []
    evaluate_scores.append(hamming_loss_score)
    evaluate_scores.append(jaccard_similarity)
    evaluate_scores.append(auc_micro)
    evaluate_scores.append(auc_macro)
    evaluate_scores.append(average_precision_micro)
    evaluate_scores.append(average_precision_macro)
    evaluate_scores.append(f1_score_micro)
    evaluate_scores.append(f1_score_macro)
    scores_output_fp = dir_path + 'compare_models_results/'+str(kfold)+'-fold_compared_models_evaluation_indicators_results.csv'
    with open(scores_output_fp, 'a') as output_file:
        temp_values = [str(x) for x in evaluate_scores]
        output_file.write(algorithm + ',' + ','.join(temp_values) + '\n')

# evaluation models
def evaluation_models(dir_path,kfold,algorithm):

    train_test_dataset_x, train_test_dataset_y = prepare_dataset_for_train(dir_path)
    X = np.array(train_test_dataset_x)
    Y = np.array(train_test_dataset_y)
    kf = KFold(n_splits=kfold)

    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        Y_train, Y_test = Y[train_index], Y[test_index]
        train_model(dir_path, X_train, X_test, Y_train, Y_test, algorithm)


if __name__ == '__main__':

    dir_path = ''
    algorithms = ['MK','XGB','NB','LR','SVM','RF']
    # algorithms = ['XGB','LR']
    kfold = 5
    for algorithm in algorithms:
        evaluation_models(dir_path, kfold, algorithm)
