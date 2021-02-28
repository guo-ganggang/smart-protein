#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 17/11/2018 9:13 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : data_helper.py
# @Software: Mining from a specific set of proteins in human sperm

import codecs
from itertools import islice
import numpy as np
from os import listdir

import sys

# reload(sys)
# sys.setdefaultencoding('utf8')


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

# Protein sequences associated with candidate sets
def popular_organisms_protein_sequence_3gram():

    corpus_fp = '/vector_protein_seq/raw_data/'

    protein_sequence_3gram = set()

    file_name_feature = '.csv'
    file_names = [f for f in listdir(corpus_fp) if f.endswith(file_name_feature)]
    print('The number of the files in this dir is %d' % len(file_names))

    for file_name in file_names:
        if file_name == 'other_organisms_protein_uniprot_entry_geneid_sequence_pfam.csv':
            continue
        protein_sequence_fp = corpus_fp + file_name
        with codecs.open(protein_sequence_fp, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                if len(temp) != 3:
                    continue
                ngram_seq = split_ngrams(temp[2].strip().upper(),3)
                for ngrams in ngram_seq:
                    for ngram in ngrams:
                        protein_sequence_3gram.add(ngram)
    print(len(protein_sequence_3gram))

    return protein_sequence_3gram


def get_ngram_vectors(ngram_corpus_vector):
    # protein_sequence_3gram = popular_organisms_protein_sequence_3gram()
    ngram_vectors = {}
    with open(ngram_corpus_vector) as infile:
        for line in infile:
            line_parts = line.rstrip().split()
            # skip first line with metadata in word2vec text file format
            if len(line_parts) > 2:
                ngram, vector_values = line_parts[0], line_parts[1:]
                # if ngram not in protein_sequence_3gram:
                #     continue
                # print(vector_values)
                ngram_vectors[ngram] = np.array(map(float, vector_values), dtype=np.float32)
    return ngram_vectors
