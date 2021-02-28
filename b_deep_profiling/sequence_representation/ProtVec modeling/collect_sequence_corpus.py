#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 27/12/2018 6:13 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : collect_sequence_corpus.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import gzip
from Bio import SeqIO
import codecs
from itertools import islice
import data_helper as hhgt
from os import listdir
import re
import os

import sys

reload(sys)
sys.setdefaultencoding('utf8')

species_dict = {'human': '9606', 'zebrafish': '7955', 'fruit_fly': '7227','rat': '10116', 'mouse': '10090'}


# Prepare the SEQUENCE training corpus, focusing on proteins reviewed on UNIPROT
def generate_sequence_corpus_file(dir_path):

    corpus_fp = dir_path + 'seq_corpus/raw_data/'
    seq_corpus_fn = corpus_fp + 'protein_sequence_including_entry_ensemble.csv'
    first_batch_uniprot_entries = set()

    file_names = [f for f in listdir(corpus_fp) if f.endswith('fasta.gz')]
    f = open(seq_corpus_fn, "a")
    for file_name in file_names:
        species_tax_id = ''
        # if 'sprot' in file_name:
        #     continue
        for tax_name in species_dict.keys():
            if tax_name in file_name:
                species_tax_id = species_dict[tax_name]
                break
        uniprot_human_corpus_fname = corpus_fp + file_name
        uniprot_entries_set = set()
        # _, uniprot_other_species_homologous_gene_entry_human_gene_entry = hhgt.other_species_homologous_gene_names_ids(species_tax_id)
        # species_protein_symbol_entries = set()
        # for protein_symbol in uniprot_other_species_homologous_gene_entry_human_gene_entry.keys():
        #     for entry in uniprot_other_species_homologous_gene_entry_human_gene_entry[protein_symbol]:
        #         species_protein_symbol_entries.add(entry)
        with gzip.open(uniprot_human_corpus_fname, 'rb', 'utf-8') as gzipped_file:
            # first_record = SeqIO.parse(gzipped_file, "fasta").next()
            # print(first_record)
            for record in SeqIO.parse(gzipped_file, "fasta"):
                # lstrip：用来去除开头字符、空白符(包括\n、\r、\t、' '，即：换行、回车、制表符、空格)
                if len(record.name.split('|')) < 2:
                    continue
                uniprot_entry = record.name.split('|')[1].strip()
                # if uniprot_entry not in species_protein_symbol_entries:
                #     continue
                if uniprot_entry in uniprot_entries_set:
                    continue
                uniprot_entries_set.add(uniprot_entry)
                first_batch_uniprot_entries.add(uniprot_entry)
                protein_seq = str(record.seq).strip()
                if not protein_seq.isalpha():
                    print(uniprot_entry+'------------------')
                f.write(uniprot_entry + ',' + protein_seq + "\n")
        print(species_tax_id + ': ' + str(len(uniprot_entries_set)))
    f.close()

    f = open(seq_corpus_fn, "a")
    uniprot_corpus_fname = corpus_fp + 'uniprot_sprot.fa.gz'
    with gzip.open(uniprot_corpus_fname, 'rb', 'utf-8') as gzipped_file:
        # first_record = SeqIO.parse(gzipped_file, "fasta").next()
        # print(first_record)
        for record in SeqIO.parse(gzipped_file, "fasta"):
            # lstrip：用来去除开头字符、空白符(包括\n、\r、\t、' '，即：换行、回车、制表符、空格)
            if len(record.name.split('|')) < 2:
                continue
            uniprot_entry = record.name.split('|')[1].strip()
            if uniprot_entry in first_batch_uniprot_entries:
                continue
            first_batch_uniprot_entries.add(uniprot_entry)
            protein_seq = str(record.seq).strip()
            if not protein_seq.isalpha():
                print(uniprot_entry + '------------------')
            f.write(uniprot_entry + ',' + protein_seq + "\n")
    print(len(first_batch_uniprot_entries))
    f.close()

    generate_sequence_corpus_file_supplyment_by_geneid(dir_path, seq_corpus_fn)

# Retrieve only NCBI geneid protein sequence, filter out the string with *
def generate_sequence_corpus_file_supplyment_by_geneid(dir_path,seq_corpus_fn):

    subdir_path = dir_path + '/species_protein_pair_score/'
    dir_file_names = [f for f in listdir(subdir_path) if f.startswith('protein_ensembl')]
    protein_ensemble_ids_supplyment_by_geneid = set()
    for file_name in dir_file_names:
        if file_name == 'protein_ensembl_identifier_homologous_gene_human_symbol':
            continue
        with codecs.open(subdir_path+file_name, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                if len(temp) < 3:
                    continue
                if temp[1].strip() != '':
                    continue
                if ';' not in temp[2].strip():
                    protein_ensemble_ids_supplyment_by_geneid.add(temp[0].strip())
                else:
                    for ensemble_id in temp[2].strip().split(';'):
                        if ensemble_id.strip() == '':
                            continue
                        protein_ensemble_ids_supplyment_by_geneid.add(ensemble_id.strip())
    print('protein_ensemble_ids_supplyment_by_geneid: ' + str(len(protein_ensemble_ids_supplyment_by_geneid)))

    # for ensemble_id in protein_ensemble_ids_supplyment_by_geneid:
    #     print(ensemble_id)

    including_ensemble_ids = set()
    ensemble_corpus_fp = dir_path + 'seq_corpus/raw_data/'
    dir_file_names = [f for f in listdir(ensemble_corpus_fp) if f.startswith('ensemble')]

    f = open(seq_corpus_fn, "a")
    for file_name in dir_file_names:
        with gzip.open(ensemble_corpus_fp+file_name, 'rb', 'utf-8') as gzipped_file:
            # first_record = SeqIO.parse(gzipped_file, "fasta").next()
            # print(first_record)
            for record in SeqIO.parse(gzipped_file, "fasta"):
                # lstrip：用来去除开头字符、空白符(包括\n、\r、\t、' '，即：换行、回车、制表符、空格)
                ensemble_id = record.id.strip().split('.')[0].strip()
                if ensemble_id not in protein_ensemble_ids_supplyment_by_geneid:
                    continue
                including_ensemble_ids.add(ensemble_id)
                protein_seq = str(record.seq).strip()
                if '*' in protein_seq:
                    continue
                if not protein_seq.isalpha(): #.isalpha(),.isupper()
                    print(ensemble_id+'+++++++++++++++')
                f.write(ensemble_id + ',' + protein_seq + "\n")
        print(file_name+' : '+str(len(including_ensemble_ids)))
    f.close()

    print('including_ensemble_ids: ' + str(len(including_ensemble_ids)))

# Prepare training data sequence vectors
def prepare_proein_sequence_vector_for_train(dir_path):

    subdir_path = dir_path + 'seq_corpus/protein_symbols_for_train/'
    dir_file_names = [f for f in listdir(subdir_path) if f.startswith('protein_ensembl')]
    ncbi_id_human_protein_ensemble_ids = {}
    uniprot_entry_human_protein_ensemble_ids = {}
    ncbi_ids = set()
    uniprot_entries = set()
    with codecs.open(subdir_path + 'protein_ensembl_identifier_homologous_gene_human', "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) < 3:
                continue
            if temp[1].strip() != '':
                if ';' not in temp[1].strip():
                    if temp[1].strip() not in uniprot_entries:
                        uniprot_entry_human_protein_ensemble_ids[temp[1].strip()] = set()
                        uniprot_entries.add(temp[1].strip())
                    uniprot_entry_human_protein_ensemble_ids[temp[1].strip()].add(temp[0].strip())
                else:
                    for uniprot_entry in temp[1].strip().split(';'):
                        if uniprot_entry.strip() == '':
                            continue
                        if uniprot_entry.strip() not in uniprot_entries:
                            uniprot_entry_human_protein_ensemble_ids[uniprot_entry.strip()] = set()
                            uniprot_entries.add(uniprot_entry.strip())
                        uniprot_entry_human_protein_ensemble_ids[uniprot_entry.strip()].add(temp[0].strip())
            if temp[2].strip() != '':
                if ';' not in temp[2].strip():
                    if temp[2].strip() not in ncbi_ids:
                        ncbi_id_human_protein_ensemble_ids[temp[2].strip()] = set()
                        ncbi_ids.add(temp[2].strip())
                    ncbi_id_human_protein_ensemble_ids[temp[2].strip()].add(temp[0].strip())
                else:
                    for ncbi_id in temp[2].strip().split(';'):
                        if ncbi_id.strip() == '':
                            continue
                        if ncbi_id.strip() not in ncbi_ids:
                            ncbi_id_human_protein_ensemble_ids[ncbi_id.strip()] = set()
                            ncbi_ids.add(ncbi_id.strip())
                        ncbi_id_human_protein_ensemble_ids[ncbi_id.strip()].add(temp[0].strip())
    print('ncbi_id_human_protein_ensemble_ids: ' + str(
        len(ncbi_id_human_protein_ensemble_ids.keys())))
    print('uniprot_entry_human_protein_ensemble_ids: ' + str(
        len(uniprot_entry_human_protein_ensemble_ids.keys())))

    human_ensemble_id_other_species_ensemble_ids = {}
    human_ensemble_ids = set()
    for file_name in dir_file_names:
        if file_name == 'protein_ensembl_identifier_homologous_gene_human':
           continue
        with codecs.open(subdir_path + file_name, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                if len(temp) < 3:
                    continue
                if temp[1].strip() != '':
                    if ';' not in temp[1].strip():
                        if temp[1].strip() in uniprot_entries:
                            temp_human_ensemble_ids = uniprot_entry_human_protein_ensemble_ids[temp[1].strip()]
                            for temp_human_ensemble_id in temp_human_ensemble_ids:
                                if temp_human_ensemble_id not in human_ensemble_ids:
                                    human_ensemble_id_other_species_ensemble_ids[temp_human_ensemble_id] = set()
                                    human_ensemble_ids.add(temp_human_ensemble_id)
                                human_ensemble_id_other_species_ensemble_ids[temp_human_ensemble_id].add(temp[0].strip())
                    else:
                        for uniprot_entry in temp[1].strip().split(';'):
                            if uniprot_entry.strip() == '':
                                continue
                            if uniprot_entry.strip() in uniprot_entries:
                                temp_human_ensemble_ids = uniprot_entry_human_protein_ensemble_ids[uniprot_entry.strip()]
                                for temp_human_ensemble_id in temp_human_ensemble_ids:
                                    if temp_human_ensemble_id not in human_ensemble_ids:
                                        human_ensemble_id_other_species_ensemble_ids[temp_human_ensemble_id] = set()
                                        human_ensemble_ids.add(temp_human_ensemble_id)
                                    human_ensemble_id_other_species_ensemble_ids[temp_human_ensemble_id].add(temp[0].strip())

                if temp[2].strip() != '':
                    if ';' not in temp[2].strip():
                        if temp[2].strip() in ncbi_ids:
                            temp_human_ensemble_ids = ncbi_id_human_protein_ensemble_ids[temp[2].strip()]
                            for temp_human_ensemble_id in temp_human_ensemble_ids:
                                if temp_human_ensemble_id not in human_ensemble_ids:
                                    human_ensemble_id_other_species_ensemble_ids[temp_human_ensemble_id] = set()
                                    human_ensemble_ids.add(temp_human_ensemble_id)
                                human_ensemble_id_other_species_ensemble_ids[temp_human_ensemble_id].add(temp[0].strip())

                    else:
                        for ncbi_id in temp[2].strip().split(';'):
                            if ncbi_id.strip() == '':
                                continue
                            if ncbi_id.strip() in ncbi_ids:
                                temp_human_ensemble_ids = ncbi_id_human_protein_ensemble_ids[ncbi_id.strip()]
                                for temp_human_ensemble_id in temp_human_ensemble_ids:
                                    if temp_human_ensemble_id not in human_ensemble_ids:
                                        human_ensemble_id_other_species_ensemble_ids[temp_human_ensemble_id] = set()
                                        human_ensemble_ids.add(temp_human_ensemble_id)
                                    human_ensemble_id_other_species_ensemble_ids[temp_human_ensemble_id].add(temp[0].strip())

    print('human_ensemble_id_other_species_ensemble_ids: ' + str(len(human_ensemble_id_other_species_ensemble_ids.keys())))

    _, human_protein_official_symbol_uniprot_entries,_,_ = hhgt.prepare_human_protein_symbols_ids()

    train_protein_symbol_ensemble_ids = {}
    miss_train_protein_symbol_geneid_entry = {}
    miss_train_protein_symbols = set()
    mapped_protein_ensemble_id_by_hand = {'C12orf80':'ENSG00000257137','FTH1P18':'ENSG00000243048','RPSAP58':'ENSG00000225178'}
    train_protein_symbol_vector_fp = subdir_path + 'train_protein_official_symbols_ensemble_ids.csv'
    with codecs.open(train_protein_symbol_vector_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if temp[0].strip() in mapped_protein_ensemble_id_by_hand.keys():
                train_protein_symbol_ensemble_ids[temp[0].strip()] = set()
                train_protein_symbol_ensemble_ids[temp[0].strip()].add(mapped_protein_ensemble_id_by_hand[temp[0].strip()])
                continue
            if temp[2].strip() != '':
                ensemble_ids = temp[2].strip().split(';')
                train_protein_symbol_ensemble_ids[temp[0].strip()] = set()
                for ensemble_id in ensemble_ids:
                    train_protein_symbol_ensemble_ids[temp[0].strip()].add(ensemble_id)
                    if ensemble_id in human_ensemble_id_other_species_ensemble_ids.keys():
                        train_protein_symbol_ensemble_ids[temp[0].strip()] = train_protein_symbol_ensemble_ids[temp[0].strip()] | human_ensemble_id_other_species_ensemble_ids[ensemble_id]
            else:
                if temp[0].strip() in human_protein_official_symbol_uniprot_entries.keys():
                    miss_train_protein_symbol_geneid_entry[temp[0].strip()] = set(human_protein_official_symbol_uniprot_entries[temp[0].strip()])
                else:
                    miss_train_protein_symbols.add(temp[0].strip())
                    print(temp[0].strip())
    print('train_protein_symbol_ensemble_ids: ' + str(len(train_protein_symbol_ensemble_ids.keys())))
    print('miss_train_protein_symbol_geneid_entry: ' + str(len(miss_train_protein_symbol_geneid_entry.keys())))
    print('miss_train_protein_symbols: ' + str(len(miss_train_protein_symbols)))

    return train_protein_symbol_ensemble_ids,miss_train_protein_symbol_geneid_entry

# Prepare sequence correspondence based on Ensemble_id and Protein_Entry respectively
def mapped_ensemble_uniprot_protein_id_sequence(dir_path):

    ensemble_protein_id_sequence = {}
    uniprot_protein_id_sequence = {}
    corpus_fp = dir_path + 'seq_corpus/raw_data/'
    dir_file_names = [f for f in listdir(corpus_fp) if f.endswith('fasta.gz')]
    print('dir_file_names: ' + str(len(dir_file_names)))
    for file_name in dir_file_names:
        uniprot_human_corpus_fname = corpus_fp + file_name
        with gzip.open(uniprot_human_corpus_fname, 'rb', 'utf-8') as gzipped_file:
            # first_record = SeqIO.parse(gzipped_file, "fasta").next()
            # print(first_record)
            for record in SeqIO.parse(gzipped_file, "fasta"):
                # lstrip：用来去除开头字符、空白符(包括\n、\r、\t、' '，即：换行、回车、制表符、空格)
                if len(record.name.split('|')) < 2:
                    continue
                uniprot_entry = record.name.split('|')[1].strip()
                protein_seq = re.sub('\n','',str(record.seq)).strip()
                uniprot_protein_id_sequence[uniprot_entry] = protein_seq
    print('uniprot_protein_id_sequence: ' + str(len(uniprot_protein_id_sequence.keys())))

    dir_file_names = [f for f in listdir(corpus_fp) if f.startswith('ensemble')]
    print('dir_file_names: '+ str(len(dir_file_names)))

    for file_name in dir_file_names:
        ensemble_corpus_fname = corpus_fp + file_name
        with gzip.open(ensemble_corpus_fname, 'rb', 'utf-8') as gzipped_file:
            # first_record = SeqIO.parse(gzipped_file, "fasta").next()
            # print(first_record)
            for record in SeqIO.parse(gzipped_file, "fasta"):
                # lstrip：用来去除开头字符、空白符(包括\n、\r、\t、' '，即：换行、回车、制表符、空格)
                ensemble_id = record.id.strip().split('.')[0].strip()
                protein_seq = re.sub('\n','',str(record.seq)).strip()
                ensemble_protein_id_sequence[ensemble_id] = protein_seq
    print('ensemble_protein_id_sequence: ' + str(len(ensemble_protein_id_sequence.keys())))

    return uniprot_protein_id_sequence,ensemble_protein_id_sequence

# Generate protein names and sequence lists
def generate_protein_symbol_sequences(dir_path):
    train_protein_symbol_equences = {}

    uniprot_protein_id_sequence, ensemble_protein_id_sequence = mapped_ensemble_uniprot_protein_id_sequence(dir_path)
    train_protein_symbol_ensemble_ids, miss_train_protein_symbol_geneid_entry = prepare_proein_sequence_vector_for_train(dir_path)

    uniprot_protein_ids = set(uniprot_protein_id_sequence.keys())
    ensemble_protein_ids = set(ensemble_protein_id_sequence.keys())

    for train_protein_symbol in train_protein_symbol_ensemble_ids.keys():
        train_protein_symbol_equences[train_protein_symbol] = set()
        ensemble_ids = train_protein_symbol_ensemble_ids[train_protein_symbol]
        for ensemble_id in ensemble_ids:
            if ensemble_id.strip() in ensemble_protein_ids:
                temp_sequence = ensemble_protein_id_sequence[ensemble_id.strip()]
                if temp_sequence == '':
                    continue
                if (not temp_sequence.isalpha()) or (not temp_sequence.isupper()):
                    continue
                train_protein_symbol_equences[train_protein_symbol].add(temp_sequence)
        if len(train_protein_symbol_equences[train_protein_symbol]) == 0:
            print(train_protein_symbol+' ********************')
    print('train_protein_symbol_equences: ' + str(len(train_protein_symbol_equences.keys())))

    for train_protein_symbol in miss_train_protein_symbol_geneid_entry.keys():
        train_protein_symbol_equences[train_protein_symbol] = set()
        uniprot_entries = miss_train_protein_symbol_geneid_entry[train_protein_symbol]
        for uniprot_entry in uniprot_entries:
            if uniprot_entry.strip() in uniprot_protein_ids:
                temp_sequence = uniprot_protein_id_sequence[uniprot_entry.strip()]
                if temp_sequence == '':
                    continue
                if (not temp_sequence.isalpha()) or (not temp_sequence.isupper()):
                    continue
                train_protein_symbol_equences[train_protein_symbol].add(temp_sequence)
        if len(train_protein_symbol_equences[train_protein_symbol]) == 0:
            print(train_protein_symbol + ' ********************')
    print('train_protein_symbol_equences: ' + str(len(train_protein_symbol_equences.keys())))

    train_protein_symbol_equences_fname = dir_path + 'seq_corpus/train_protein_symbol_equences.csv'
    if not os.path.exists(train_protein_symbol_equences_fname):
        with open(train_protein_symbol_equences_fname, 'w') as output_file:
            for key in train_protein_symbol_equences.keys():
                if len(train_protein_symbol_equences[key]) == 0:
                    continue
                sequences = ';'.join(train_protein_symbol_equences[key])
                output_file.write(key+'$'+sequences+'\n')




if __name__ == '__main__':

    dir_path = ''
    # generate_sequence_corpus_file(dir_path)
    generate_protein_symbol_sequences(dir_path)