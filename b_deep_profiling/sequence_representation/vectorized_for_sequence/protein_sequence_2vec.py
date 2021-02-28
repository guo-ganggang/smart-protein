#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 11/11/2018 4:08 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : protein_sequence_2vec.py
# @Software: Mining from a specific set of proteins in human sperm


import warnings
warnings.filterwarnings(action='ignore', category=UserWarning, module='gensim')
import platform
import logging
import os.path
import sys
import multiprocessing
from gensim.models import Word2Vec
from gensim.models.word2vec import Text8Corpus
import numpy as np
import codecs
from os import listdir


def to_vecs(n,size,protein_seq,ngram_vectors):

    """
    convert sequence to three n-length vectors
    e.g. 'AGAMQSASM' => [ array([  ... * 100 ], array([  ... * 100 ], array([  ... * 100 ] ]
    """

    # ngram_patterns = data_preparation.split_ngrams(protein_seq, n)
    # protvecs = []
    # for ngrams in ngram_patterns:
    #     ngram_vecs = []
    #     for ngram in ngrams:
    #         try:
    #             ngram_vecs.append(ngram_vectors[ngram])
    #         except:
    #             raise Exception("Model has never trained this n-gram: " + ngram)
    #     protvecs.append(sum(ngram_vecs))
    # return protvecs

    """
    convert sequence to one n-length vectors
    e.g. 'AGAMQSASM' => 'AGA GAM AMQ MQS QSA SAS ASM' +=> array([  ... * 100 ]
    """
    protvec = np.zeros(size, dtype=np.float32)
    for index in xrange(len(protein_seq) + 1 - n):
        ngram = protein_seq[index:index + n]
        if ngram in ngram_vectors:
            ngram_vector = ngram_vectors[ngram]
            protvec += ngram_vector
    return protvec

def get_ngram_vectors(ngram_corpus_vector):
    ngram_vectors = {}
    with open(ngram_corpus_vector) as infile:
        for line in infile:
            line_parts = line.rstrip().split()
            # skip first line with metadata in word2vec text file format
            if len(line_parts) > 2:
                ngram, vector_values = line_parts[0], line_parts[1:]
                ngram_vectors[ngram] = np.array(map(float, vector_values), dtype=np.float32)
    return ngram_vectors

def split_ngrams(seq, n):
    """
    'AGAMQSASM' => [['AGA', 'MQS', 'ASM'], ['GAM','QSA'], ['AMQ', 'SAS']]
    """
    a, b, c = zip(*[iter(seq)]*n), zip(*[iter(seq[1:])]*n), zip(*[iter(seq[2:])]*n)
    str_ngrams = []
    for ngrams in [a,b,c]:
        x = []
        for ngram in ngrams:
            x.append("".join(ngram))
        str_ngrams.append(x)
    return str_ngrams

def generate_corpusfile(dir_path, n, out_corpus,path_type):

    protein_sequence_entry_name = set()
    protein_sequence_fp = dir_path + 'sequence_corpus'+path_type + 'protein_sequence_including_entry_ensemble.csv'
    with open(out_corpus, 'a') as output_file:
        with codecs.open(protein_sequence_fp, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split(',')
                if len(temp) != 2:
                    continue
                ngram_patterns = split_ngrams(temp[1].strip().upper(), n)
                for ngram_pattern in ngram_patterns:
                    output_file.write(" ".join(ngram_pattern) + '\n')
                protein_sequence_entry_name.add(temp[0].strip())
    print(len(protein_sequence_entry_name))

def train_w2v_model(file_path,path_type):
    print('data preparation ...')
    n = 3
    out_corpus = file_path + 'sequence_corpus'+path_type+'protein_sequence_ngram_corpus.csv'

    if (not os.path.isfile(out_corpus)) or (os.path.getsize(out_corpus) < 5):
        print ('INFORM : There is no corpus file. Generate Corpus file from fasta file...')
        generate_corpusfile(file_path, n, out_corpus,path_type)
    else:
        print("INFORM : File's Existence is confirmed")

    print('train model ...')
    program = os.path.basename(sys.argv[0])
    logger = logging.getLogger(program)

    logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s')
    logging.root.setLevel(level=logging.INFO)
    logger.info("running %s" % ' '.join(sys.argv))

    # check and process input arguments

    out_model = file_path + 'bio2vec_model'+path_type+'protein_sequence_corpus_3letter_sg1_s128_w25_m1.bin'
    out_vec = file_path + 'bio2vec_model'+path_type+'protein_sequence_corpus_3letter_sg1_s128_w25_m1.vector'

    if os.path.isfile(out_model) and os.path.isfile(out_vec):
        print(os.path.getsize(out_model),os.path.getsize(out_vec))
    else:
        model = Word2Vec(Text8Corpus(out_corpus), sg=1, size=128, window=25, min_count=1,workers=multiprocessing.cpu_count()-1)

        # trim unneeded model memory = use(much) less RAM
        # model.init_sims(replace=True)
        model.save(out_model)
        model.wv.save_word2vec_format(out_vec, binary=False)

    return out_model,out_vec

def make_protein_vector_for_uniprot(file_path,path_type):

    ngram_protein_vector_fname = file_path + 'train_protein_seq_ngram_vector.csv'
    n = 3
    size = 128

    _,ngram_corpus_vector_fname = train_w2v_model(file_path,path_type)

    print('load model ...')
    ngram_vectors = get_ngram_vectors(ngram_corpus_vector_fname)

    print('make protein vector ...')

    train_protein_sequence_fp = file_path + 'sequence_corpus' + path_type + 'train_protein_symbol_equences.csv'
    with open(ngram_protein_vector_fname, 'a') as output_file:
        with codecs.open(train_protein_sequence_fp, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                protein_symbol = temp[0].strip()
                protein_sequences = temp[1].strip().split(';')
                sum_protein_vectors = np.zeros((128))
                temp_count = 0
                for protein_sequence in protein_sequences:
                    if protein_sequence.strip() == '':
                        continue
                    protein_vector = to_vecs(n,size,protein_sequence,ngram_vectors)
                    sum_protein_vectors += protein_vector
                    temp_count += 1
                avg_protein_vectors = sum_protein_vectors / float(temp_count)
                output_file.write(protein_symbol+'$'+ ','.join(map(str, avg_protein_vectors))+'\n')


if __name__ == '__main__':
    print(platform.system())
    file_path = ''
    path_type = '/'
    if platform.system() == 'Windows':
        file_path = 'vector_protein_seq/'
        path_type = '/'

    make_protein_vector_for_uniprot(file_path,path_type)