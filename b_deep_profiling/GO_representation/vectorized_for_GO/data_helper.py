#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 28/12/2017 下午 8:43
# @Author  : GUO Ganggang
# @email   : ganggangguo@csu.edu.cn
# @Site    : 
# @File    : protein_literature_2vec.py
# @Software: PyCharm

import warnings
warnings.filterwarnings(action='ignore', category=UserWarning, module='gensim')

import nltk
import nltk.data
from nltk.corpus import stopwords
import os
import codecs
import logging
import os.path
import sys
import multiprocessing
from gensim.models import Word2Vec
from gensim.models.word2vec import Text8Corpus
import numpy as np

# After the standardization of the text clause cleaning
def clean_literatures_corpus(dir_path):
    # Load the punkt tokenizer
    download_nltk_corpus_path = '/nltk_data/'
    if not os.path.exists(download_nltk_corpus_path):
        os.makedirs(download_nltk_corpus_path)
        nltk.download() # 下载到上面那个路径下
    english_stopwords = stopwords.words('english') # 自动检索上面那个路径
    english_punctuations = [',', '.', ':', ';', '?', '(', ')', '[', ']', '&', '!', '*', '@', '#', '$', '%',
                            '>', '<', '\'s', '\'', '-/-', '``', "''", '--', '-']

    raw_text_file_path = dir_path + 'text_corpus/corpus_literature_definition_for_go.csv'
    clean_text_file_path = dir_path + 'text_corpus/corpus_literature_definition_for_go_clean.csv'
    with open(clean_text_file_path, 'w') as output_file:
        with codecs.open(raw_text_file_path, "rb", "utf-8") as input_file:
            for line in input_file:
                corpus = line.strip()
                texts_tokenized = [nltk.word_tokenize(sentence) for sentence in nltk.sent_tokenize(corpus)]
                texts_filtered_stopwords = [[word for word in document if not word in english_stopwords] for document in \
                                            texts_tokenized]
                texts_filtered = [[word for word in document if not word in english_punctuations] for document in \
                                  texts_filtered_stopwords]
                # texts_stemmed = [[st.stem(word) for word in docment] for docment in texts_filtered]
                for text in texts_filtered:
                    if len(text) < 2:
                        continue
                    words = []
                    for word in text:
                        if word.strip() == '':
                            continue
                        if word.strip().istitle():
                            words.append(word.strip().lower())
                        else:
                            words.append(word.strip())
                    sent = ' '.join(words)
                    output_file.write(sent + '\n')

# After the standardization of the text clause cleaning
def clean_go_term_corpus(dir_path):
    # Load the punkt tokenizer
    download_nltk_corpus_path = '/nltk_data/'
    if not os.path.exists(download_nltk_corpus_path):
        os.makedirs(download_nltk_corpus_path)
        nltk.download() # 下载到上面那个路径下
    english_stopwords = stopwords.words('english') # 自动检索上面那个路径
    english_punctuations = [',', '.', ':', ';', '?', '(', ')', '[', ']', '&', '!', '*', '@', '#', '$', '%',
                            '>', '<', '\'s', '\'', '-/-', '``', "''", '--', '-']

    raw_text_file_path = dir_path + 'text_corpus/train_go_id_term.csv'
    clean_text_file_path = dir_path + 'text_corpus/train_go_id_term_clean.csv'
    with open(clean_text_file_path, 'w') as output_file:
        with codecs.open(raw_text_file_path, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                corpus = temp[1].strip()
                texts_tokenized = [nltk.word_tokenize(sentence) for sentence in nltk.sent_tokenize(corpus)]
                texts_filtered_stopwords = [[word for word in document if not word in english_stopwords] for document in \
                                            texts_tokenized]
                texts_filtered = [[word for word in document if not word in english_punctuations] for document in \
                                  texts_filtered_stopwords]
                # texts_stemmed = [[st.stem(word) for word in docment] for docment in texts_filtered]
                for text in texts_filtered:
                    if len(text) < 2:
                        continue
                    words = []
                    for word in text:
                        if word.strip() == '':
                            continue
                        if word.strip().istitle():
                            words.append(word.strip().lower())
                        else:
                            words.append(word.strip())
                    sent = ' '.join(words)
                    if sent.strip() == '':
                        continue
                    output_file.write(temp[0].strip() + '$' + sent + '\n')

# The Word2Vec model was used for training
def train_w2v_model(dir_path):
    print('data preparation ...')

    clean_text_corpus = dir_path + 'text_corpus/corpus_literature_definition_for_go_clean.csv'
    if (not os.path.isfile(clean_text_corpus)) or (os.path.getsize(clean_text_corpus) < 100):
        print ('INFORM : There is no corpus file. Generate Corpus file from clean_text_corpus...')
        clean_literatures_corpus(dir_path)
    else:
        print("INFORM : File's Existence is confirmed")

    clean_go_term = dir_path + 'text_corpus/train_go_id_term_clean.csv'
    if (not os.path.isfile(clean_go_term)) or (os.path.getsize(clean_go_term) < 100):
        print ('INFORM : There is no corpus file. Generate Corpus file from clean_go_term_corpus...')
        clean_go_term_corpus(dir_path)
    else:
        print("INFORM : File's Existence is confirmed")

    print('train model ...')
    program = os.path.basename(sys.argv[0])
    logger = logging.getLogger(program)

    logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s')
    logging.root.setLevel(level=logging.INFO)
    logger.info("running %s" % ' '.join(sys.argv))

    # check and process input arguments

    out_model = dir_path + 'w2v_model/protein_go_term_corpus_sg1_s128_w5_m1.bin'
    out_vec = dir_path + 'w2v_model/protein_go_term_corpus_sg1_s128_w5_m1.vector'

    if os.path.isfile(out_model) and os.path.isfile(out_vec):
        print(os.path.getsize(out_model),os.path.getsize(out_vec))
    else:
        model = Word2Vec(Text8Corpus(clean_text_corpus), sg=1, size=128, window=5, min_count=1,workers=int(multiprocessing.cpu_count()-2))

        # trim unneeded model memory = use(much) less RAM
        # model.init_sims(replace=True)
        model.save(out_model)
        model.wv.save_word2vec_format(out_vec, binary=False)


if __name__ == '__main__':

    dir_path = ''
    train_w2v_model(dir_path)
