#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 21/12/2018 3:18 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : data_helper.py
# @Software: Mining from a specific set of proteins in human sperm

from itertools import islice
import os
import codecs
import pandas as pd
import re


species_dict = {'human': '9606', 'zebrafish': '7955', 'fruit_fly': '7227','rat': '10116', 'mouse': '10090'}


'''
    Prepare the human protein list, prepare the training set human protein list, 
    prepare the four species homologous gene list of human genes
'''


# Read in a list of human protein names and their abbreviations,
# and at the same time, read in a list of protein names and
# their abbreviations related to the training set
def prepare_human_protein_symbols_ids():
    dir_path = '/text_corpus/'
    raw_data_dir = '/protein_entities_standardized/'
    raw_data_fp = raw_data_dir + 'human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    human_protein_official_symbol_names = {}
    with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            human_protein_official_symbol_names[temp[0]] = set(temp[1:])
    print('human_protein_official_symbol_names: ' + str(len(human_protein_official_symbol_names.keys())))
    human_protein_official_symbols = set(human_protein_official_symbol_names.keys())

    train_data_fp = dir_path + 'prepare_gene_symbols_abbreviation_for_train_model.csv'
    train_protein_official_symbols_aliases = {}
    special_train_protein_symbols = {'NAT6': 'NAA80', 'EFTUD1': 'EFL1', 'LNP': 'LNPK', 'APITD1': 'CENPS-CORT'}
    del_train_protein_symbols = ['FAM58BP']
    with codecs.open(train_data_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if temp[0].strip() in del_train_protein_symbols:
                continue
            if temp[0].strip() in special_train_protein_symbols.keys():
                train_protein_official_symbols_aliases[special_train_protein_symbols[temp[0].strip()]] = \
                human_protein_official_symbol_names[special_train_protein_symbols[temp[0].strip()]]
                continue
            if temp[0].strip() in human_protein_official_symbols:
                train_protein_official_symbols_aliases[temp[0].strip()] = human_protein_official_symbol_names[temp[0].strip()]
            else:
                flag = 0
                for symbol in human_protein_official_symbol_names.keys():
                    if temp[0].strip() in human_protein_official_symbol_names[symbol]:
                        if (temp[0].strip() == symbol) and (len(human_protein_official_symbol_names[symbol]) == 1):
                            continue
                        train_protein_official_symbols_aliases[symbol] = human_protein_official_symbol_names[symbol]
                        flag += 1
                if flag > 1:
                    print(temp[0].strip() + ': ' + str(flag))

    print('train_protein_official_symbols_aliases: ' + str(len(train_protein_official_symbols_aliases.keys())))

    raw_data_mapped_ids_dbs_fp = dir_path + 'human_protein_official_symbol_ncbi_ids_uniprot_entries_v2.csv'
    human_protein_official_symbol_ncbi_ids = {}
    human_protein_official_symbol_uniprot_entries = {}

    # # 核对匹配上uniprot 与 ncbi id 的蛋白质集合
    # human_protein_official_symbols_mapped = set()
    identify_problems = ['HLA-DRB3$0301', 'HLA-DRB4$0101']
    # special_mapped_protein_symbols = {'hCG_2039718': 'RAD51L3-RFFL', 'hCG_1807616': 'LINC02218',
    #                                   'hCG_1809904': 'TAF11L7','LOC105371242': 'PPIAL4H',
    #                                   'WIPI3': 'WDR45B', 'hCG_2002594': 'SEPT5',
    #                                   'TMEM27': 'CLTRN', 'KIAA1024L': 'MINAR2',
    #                                   'C9orf84': 'SHOC1', 'IL-21': 'IL21',
    #                                   'hCG_1796489': 'LOC101059948'}
    # special_mapped_protein_symbols_for_check = dict([val, key] for key, val in special_mapped_protein_symbols.items())

    with codecs.open(raw_data_mapped_ids_dbs_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            for special_protein_name in identify_problems:
                if special_protein_name in line.strip():
                    temp_line = (re.sub(special_protein_name, '', line.strip())).strip().split('$')
                    if temp_line[0].strip() != '':
                        human_protein_official_symbol_ncbi_ids[special_protein_name] = temp_line[0].strip()
                    if temp_line[1].strip() != '':
                        human_protein_official_symbol_uniprot_entries[special_protein_name] = temp_line[1].strip()
                    continue
            temp = line.strip().split('$')
            if temp[1].strip() != '':
                human_protein_official_symbol_ncbi_ids[temp[0].strip()] = temp[1].strip()
            if temp[2].strip() != '':
                human_protein_official_symbol_uniprot_entries[temp[0].strip()] = temp[2].strip()
    print('human_protein_official_symbol_ncbi_ids: ' + str(len(human_protein_official_symbol_ncbi_ids.keys())))
    print('human_protein_official_symbol_uniprot_entries: ' + str(
        len(human_protein_official_symbol_uniprot_entries.keys())))
    human_protein_official_symbols_mapped = set(human_protein_official_symbol_ncbi_ids.keys()).union(
        human_protein_official_symbol_uniprot_entries.keys())
    train_protein_official_symbols_mapped = set(train_protein_official_symbols_aliases.keys()).intersection(
        human_protein_official_symbols_mapped)
    print('train_protein_official_symbols_mapped: ' + str(len(train_protein_official_symbols_mapped)))

    return human_protein_official_symbol_ncbi_ids, human_protein_official_symbol_uniprot_entries, train_protein_official_symbols_aliases,human_protein_official_symbol_names


# List of homologous genes of the four major species and human genes
def other_species_homologous_gene_names_ids(tax_id):

    ncbi_other_species_homologous_gene_id_human_gene_id = {}
    uniprot_other_species_homologous_gene_entry_human_gene_entry = {}

    homologene_results = '/homologene/results/'
    ncbi_homologous_gene_out_fp = homologene_results + 'human_homologous_gene_by_ncbi_id.csv'
    uniprot_homologous_gene_out_fp = homologene_results + 'human_homologous_gene_by_uniprot_entry.csv'

    with codecs.open(ncbi_homologous_gene_out_fp, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if len(temp) < 3:
                continue
            if temp[0].strip() != tax_id:
                continue
            if ';' not in temp[2].strip():
                ncbi_other_species_homologous_gene_id_human_gene_id[temp[2].strip()] = temp[1].strip()
            else:
                other_species_homologous_gene_ids = temp[2].strip().split(';')
                for other_species_homologous_gene_id in other_species_homologous_gene_ids:
                    if other_species_homologous_gene_id.strip() == '':
                        continue
                    ncbi_other_species_homologous_gene_id_human_gene_id[other_species_homologous_gene_id.strip()] = \
                    temp[1].strip()

    with codecs.open(uniprot_homologous_gene_out_fp, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if len(temp) < 3:
                continue
            if temp[0].strip() != tax_id:
                continue
            if ';' not in temp[2].strip():
                uniprot_other_species_homologous_gene_entry_human_gene_entry[temp[2].strip()] = temp[1].strip()
            else:
                other_species_homologous_gene_entries = temp[2].strip().split(';')
                for other_species_homologous_gene_entry in other_species_homologous_gene_entries:
                    if other_species_homologous_gene_entry.strip() == '':
                        continue
                    uniprot_other_species_homologous_gene_entry_human_gene_entry[
                        other_species_homologous_gene_entry.strip()] = \
                        temp[1].strip()
    return ncbi_other_species_homologous_gene_id_human_gene_id,uniprot_other_species_homologous_gene_entry_human_gene_entry