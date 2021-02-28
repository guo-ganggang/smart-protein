#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 3/12/2018 9:58 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : collect_ppi_score_network.py
# @Software: Mining from a specific set of proteins in human sperm

import gzip
import shutil
import codecs
from os import listdir
from itertools import islice
import os
import xlrd
import re
from Bio import Entrez
import urllib2
import socket
from Bio import Medline
import ssl
from bs4 import BeautifulSoup

###############
###Intermediate process code, for user reference only
###############

import sys

reload(sys)
sys.setdefaultencoding('utf-8')


Entrez.email = "ggguo@smu.edu.sg"
socket.setdefaulttimeout(600)

species_dict = {'human': '9606', 'mouse': '10090', 'zebrafish': '7955', 'rat': '10116', 'fruit_fly': '7227'}

'''
    1. Clean-up data of association between genes and literatures downloaded 
    from two databases, including five species; meanwhile, collect literatures 
    related to mouse genes from MGI database
'''

# UnZip File
def un_gz(file_name):
    """ungz zip file"""
    f_name = file_name.replace(".gz", "")
    with gzip.open(file_name, 'rb') as f_in:
        with open(f_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return f_name

# NCBI Literature collection
def get_gene_pubmed_list_from_ncbi(dir_path,species):

    # gene id 与相关文献的对应关系
    species_dict = {'human':'9606', 'mouse':'10090', 'zebrafish':'7955', 'rat':'10116', 'fruit_fly':'7227'}
    ncbi_gene_pmid_fp = dir_path + 'ncbi/'
    gz_file_names = [f for f in listdir(ncbi_gene_pmid_fp) if f.endswith('.gz')]
    for gz_file_name in gz_file_names:
        file_name = gz_file_name.split('.')[0].strip()
        if not os.path.exists(ncbi_gene_pmid_fp + file_name):
            un_gz(ncbi_gene_pmid_fp + gz_file_name)
    gene_pubmed_fp = dir_path + 'ncbi/gene2pubmed'
    gene_info_fp = dir_path + 'ncbi/gene_info'
    for key in species_dict.keys():
        if key != species:
            continue
        tax_id = species_dict[key]
        gene_id_pubmed_list = {}
        ncbi_gene_ids = set()
        gene_id_symbol = {}
        with codecs.open(gene_pubmed_fp, "rb", "utf-8") as input_file:
            for line in input_file:
                if '#' == line.strip()[0]:
                    continue
                temp = line.strip().split('	')
                if len(temp) != 3:
                    continue
                if temp[0] != tax_id:
                    continue
                if temp[1].strip() not in ncbi_gene_ids:
                    gene_id_pubmed_list[temp[1].strip()] = []
                    ncbi_gene_ids.add(temp[1].strip())
                gene_id_pubmed_list[temp[1].strip()].append(temp[2].strip())

        print('gene_id_pubmed_list: '+str(len(gene_id_pubmed_list.keys())))

        # 各个物种基因ID 与 基因名称之间的对应关系
        with codecs.open(gene_info_fp, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('	')
                if len(temp) < 3:
                    continue
                if temp[0].strip() != tax_id:
                    continue
                if temp[1].strip() not in ncbi_gene_ids:
                    continue
                gene_id_symbol[temp[1].strip()] = temp[2].strip()

        print('gene_id_symbol: ' + str(len(gene_id_symbol.keys())))

        return gene_id_pubmed_list,gene_id_symbol

# Uniprot Literature collection
def get_gene_pubmed_list_from_uniprot(dir_path,species):

    uniprot_gene_pmid_fp = dir_path + 'uniprot/'
    gz_file_names = [f for f in listdir(uniprot_gene_pmid_fp) if f.endswith('.gz')]
    uniprot_species_pmids = {}
    for gz_file_name in gz_file_names:
        file_name = gz_file_name.split('.')[0].strip()
        if species not in file_name:
            continue
        if not os.path.exists(uniprot_gene_pmid_fp + file_name):
            un_gz(uniprot_gene_pmid_fp + gz_file_name)
        uniprot_species_pmids[file_name] = set()
        file_name_path = un_gz(uniprot_gene_pmid_fp + gz_file_name)
        try:
            data = xlrd.open_workbook(file_name_path)
        except Exception, e:
            print str(e)
        table = data.sheet_by_name('Sheet0')
        uniprot_gene_entries = set()
        gene_pubmed_list = {}
        gene_entry_symbol = {}
        uniprot_entry_list = table.col_values(0)
        pubmed_ids_str_list = table.col_values(1)
        gene_symbol_list = table.col_values(2)
        for i in range(1,len(uniprot_entry_list)):
            if gene_symbol_list[i].strip() == '':
                continue
            gene_entry = uniprot_entry_list[i].strip()
            pubmed_ids = pubmed_ids_str_list[i].strip().split(';')
            gene_symbol_str = re.sub(' ',';',gene_symbol_list[i].strip())
            gene_symbol_str = eval(repr(gene_symbol_str).replace('/', ';'))
            gene_symbol = [gene_symbol_str.split(';')[0].strip()]
            if gene_entry not in uniprot_gene_entries:
                gene_pubmed_list[gene_entry] = []
                gene_entry_symbol[gene_entry] = []
                uniprot_gene_entries.add(gene_entry)
            gene_pubmed_list[gene_entry] += pubmed_ids
            gene_entry_symbol[gene_entry] += gene_symbol
        print('gene_pubmed_list: ' + str(len(gene_pubmed_list.keys())))
        print('gene_entry_symbol: ' + str(len(gene_entry_symbol.keys())))
        return gene_pubmed_list, gene_entry_symbol

# MGI Literature collection,Only mouse genes
def get_gene_pubmed_list_from_mgi(dir_path):

    mgi_homolo_gene_fp = dir_path + 'mgi/HMD_HumanPhenotype.rpt'
    mouse_gene_id_mgi_maker_accession_id = {}
    mouse_gene_ids = set()
    mgi_maker_accession_ids_homolo = set()
    with codecs.open(mgi_homolo_gene_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) < 6:
                continue
            if temp[5].strip() not in mouse_gene_ids:
                mouse_gene_id_mgi_maker_accession_id[temp[2].strip()] = set()
                mouse_gene_ids.add(temp[2].strip())
            mouse_gene_id_mgi_maker_accession_id[temp[2].strip()].add(temp[5].strip())
            mgi_maker_accession_ids_homolo.add(temp[5].strip())

    print('mouse_gene_id_mgi_maker_accession_id: ' + str(len(mouse_gene_id_mgi_maker_accession_id.keys())))

    gene_pubmed_fp = dir_path + 'mgi/MRK_Reference.rpt'
    mgi_maker_accession_id_pubmed_ids = {}
    mgi_maker_accession_ids = set()
    with codecs.open(gene_pubmed_fp, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) < 5:
                continue
            if temp[0].strip() not in mgi_maker_accession_ids_homolo:
                continue
            pubmed_ids = temp[4].strip().split('|')
            if temp[0].strip() not in mgi_maker_accession_ids:
                mgi_maker_accession_id_pubmed_ids[temp[0].strip()] = []
            mgi_maker_accession_id_pubmed_ids[temp[0].strip()] += pubmed_ids
    print('mgi_maker_accession_id_pubmed_ids: ' + str(len(mgi_maker_accession_id_pubmed_ids.keys())))

    mgi_homolo_gene_id_pubmed_ids = {}

    for key in mouse_gene_id_mgi_maker_accession_id.keys():
        temp_pubmed_ids = []
        for mgi_maker_accession_id in mouse_gene_id_mgi_maker_accession_id[key]:
            if mgi_maker_accession_id in mgi_maker_accession_id_pubmed_ids.keys():
                temp_pubmed_ids += mgi_maker_accession_id_pubmed_ids[mgi_maker_accession_id]
        while '' in temp_pubmed_ids:
            temp_pubmed_ids.remove('')
        if len(temp_pubmed_ids) == 0:
            continue
        mgi_homolo_gene_id_pubmed_ids[key] = temp_pubmed_ids
    print('mgi_homolo_gene_id_pubmed_ids: ' + str(len(mgi_homolo_gene_id_pubmed_ids.keys())))

    return mgi_homolo_gene_id_pubmed_ids

'''
    2、Literature related to human genes/proteins
'''

# Read in a list of human protein names and their abbreviations,
# and at the same time, read in a list of protein names and
# their abbreviations related to the training set
def prepare_human_protein_symbols_ids(dir_path):
    raw_data_fp = 'C:/ggguo/1_data_help/protein_entities_standardized/human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    human_protein_official_symbol_names = {}
    with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            human_protein_official_symbol_names[temp[0]] = set(temp[1:])
    print('human_protein_official_symbol_names: '+str(len(human_protein_official_symbol_names.keys())))
    human_protein_official_symbols = set(human_protein_official_symbol_names.keys())

    train_data_fp = dir_path + 'prepare_gene_symbols_abbreviation_for_train_model.csv'
    train_protein_official_symbols = {}
    special_train_protein_symbols = {'NAT6':'NAA80','EFTUD1':'EFL1','LNP':'LNPK','APITD1':'CENPS-CORT'}
    del_train_protein_symbols = ['FAM58BP']
    with codecs.open(train_data_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if temp[0].strip() in del_train_protein_symbols:
                continue
            if temp[0].strip() in special_train_protein_symbols.keys():

                train_protein_official_symbols[special_train_protein_symbols[temp[0].strip()]] = human_protein_official_symbol_names[special_train_protein_symbols[temp[0].strip()]]
                continue
            if temp[0].strip() in human_protein_official_symbols:
                train_protein_official_symbols[temp[0].strip()] = human_protein_official_symbol_names[temp[0].strip()]
            else:
                flag = 0
                for symbol in human_protein_official_symbol_names.keys():
                    if temp[0].strip() in human_protein_official_symbol_names[symbol]:
                        if (temp[0].strip() == symbol) and (len(human_protein_official_symbol_names[symbol]) == 1):
                            continue
                        train_protein_official_symbols[symbol] = human_protein_official_symbol_names[symbol]
                        flag += 1
                if flag > 1:
                    print(temp[0].strip()+': '+str(flag))

    print('train_protein_official_symbols: ' + str(len(train_protein_official_symbols.keys())))

    raw_data_mapped_ids_dbs_fp = dir_path + 'human_protein_official_symbol_ncbi_ids_uniprot_entries_v2.csv'
    human_protein_official_symbol_ncbi_ids = {}
    human_protein_official_symbol_uniprot_entries = {}

    # # 核对匹配上uniprot 与 ncbi id 的蛋白质集合
    # human_protein_official_symbols_mapped = set()
    # identify_problems = ['HLA-DRB3$0301','HLA-DRB4$0101']
    # special_mapped_protein_symbols = {'hCG_2039718': 'RAD51L3-RFFL', 'hCG_1807616': 'LINC02218',
    #                                   'hCG_1809904': 'TAF11L7','LOC105371242': 'PPIAL4H',
    #                                   'WIPI3': 'WDR45B', 'hCG_2002594': 'SEPT5',
    #                                   'TMEM27': 'CLTRN', 'KIAA1024L': 'MINAR2',
    #                                   'C9orf84': 'SHOC1', 'IL-21': 'IL21',
    #                                   'hCG_1796489': 'LOC101059948'}
    # special_mapped_protein_symbols_for_check = dict([val, key] for key, val in special_mapped_protein_symbols.items())
    # with codecs.open(raw_data_mapped_ids_dbs_fp, 'rb', 'utf-8') as input_file:
    #     for line in islice(input_file.readlines(), 0, None):
    #         for special_protein_name in identify_problems:
    #             if special_protein_name in line.strip():
    #                 human_protein_official_symbols_mapped.add(special_protein_name)
    #                 continue
    #         temp = line.strip().split('$')
    #         if temp[0].strip() in special_mapped_protein_symbols_for_check.keys():
    #             human_protein_official_symbols_mapped.add(special_mapped_protein_symbols_for_check[temp[0].strip()])
    #         human_protein_official_symbols_mapped.add(temp[0].strip())
    # print('human_protein_official_symbols_mapped: ' + str(len(human_protein_official_symbols_mapped)))
    #
    # skip_special_protein_names = set(human_protein_official_symbol_names.keys()).difference(
    #     human_protein_official_symbols_mapped)
    # print(skip_special_protein_names)
    # print(train_protein_official_symbols.difference(human_protein_official_symbols_mapped))

    identify_problems = ['HLA-DRB3$0301', 'HLA-DRB4$0101']
    with codecs.open(raw_data_mapped_ids_dbs_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            for special_protein_name in identify_problems:
                if special_protein_name in line.strip():
                    temp_line = (re.sub(special_protein_name,'',line.strip())).strip().split('$')
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
    print('human_protein_official_symbol_uniprot_entries: ' + str(len(human_protein_official_symbol_uniprot_entries.keys())))
    human_protein_official_symbols_mapped = set(human_protein_official_symbol_ncbi_ids.keys()).union(human_protein_official_symbol_uniprot_entries.keys())
    train_protein_official_symbols_mapped = set(train_protein_official_symbols.keys()).intersection(human_protein_official_symbols_mapped)
    print('train_protein_official_symbols_mapped: '+ str(len(train_protein_official_symbols_mapped)))

    return human_protein_official_symbol_ncbi_ids,human_protein_official_symbol_uniprot_entries,train_protein_official_symbols

def collect_human_gene_related_pubmed_corpus(dir_path):
    human_protein_official_symbol_ncbi_ids, human_protein_official_symbol_uniprot_entries, _ = prepare_human_protein_symbols_ids(dir_path)

    species = 'human'
    ncbi_gene_id_pubmed_list, _ = get_gene_pubmed_list_from_ncbi(dir_path, species)
    uniprot_gene_id_pubmed_list, _ = get_gene_pubmed_list_from_uniprot(dir_path, species)

    human_protein_official_symbol_pubmed_ids = {}
    human_protein_official_symbols_temp = set()
    for human_protein_official_symbol in human_protein_official_symbol_ncbi_ids.keys():
        ncbi_id_str = human_protein_official_symbol_ncbi_ids[human_protein_official_symbol]

        if ';' not in ncbi_id_str:
            ncbi_id = ncbi_id_str
            if ncbi_id.strip() in ncbi_gene_id_pubmed_list.keys():
                if human_protein_official_symbol not in human_protein_official_symbols_temp:
                    human_protein_official_symbols_temp.add(human_protein_official_symbol)
                    human_protein_official_symbol_pubmed_ids[human_protein_official_symbol] = ncbi_gene_id_pubmed_list[ncbi_id.strip()]
                else:
                    human_protein_official_symbol_pubmed_ids[human_protein_official_symbol] += ncbi_gene_id_pubmed_list[
                        ncbi_id.strip()]
        else:
            ncbi_ids = ncbi_id_str.split(';')
            for ncbi_id in ncbi_ids:
                if ncbi_id.strip() == '':
                    continue
                if ncbi_id.strip() in ncbi_gene_id_pubmed_list.keys():
                    if human_protein_official_symbol not in human_protein_official_symbols_temp:
                        human_protein_official_symbols_temp.add(human_protein_official_symbol)
                        human_protein_official_symbol_pubmed_ids[human_protein_official_symbol] = ncbi_gene_id_pubmed_list[
                            ncbi_id.strip()]
                    else:
                        human_protein_official_symbol_pubmed_ids[human_protein_official_symbol] += ncbi_gene_id_pubmed_list[
                            ncbi_id.strip()]

    for human_protein_official_symbol in human_protein_official_symbol_uniprot_entries.keys():
        uniprot_entry_str = human_protein_official_symbol_uniprot_entries[human_protein_official_symbol]

        if ';' not in uniprot_entry_str:
            uniprot_entry = uniprot_entry_str
            if uniprot_entry.strip() in uniprot_gene_id_pubmed_list.keys():
                if human_protein_official_symbol not in human_protein_official_symbols_temp:
                    human_protein_official_symbols_temp.add(human_protein_official_symbol)
                    human_protein_official_symbol_pubmed_ids[human_protein_official_symbol] = uniprot_gene_id_pubmed_list[
                        uniprot_entry.strip()]
                else:
                    human_protein_official_symbol_pubmed_ids[human_protein_official_symbol] += uniprot_gene_id_pubmed_list[
                        uniprot_entry.strip()]
        else:
            uniprot_entries = uniprot_entry_str.split(';')
            for uniprot_entry in uniprot_entries:
                if uniprot_entry.strip() == '':
                    continue
                if uniprot_entry.strip() in uniprot_gene_id_pubmed_list.keys():
                    if human_protein_official_symbol not in human_protein_official_symbols_temp:
                        human_protein_official_symbols_temp.add(human_protein_official_symbol)
                        human_protein_official_symbol_pubmed_ids[human_protein_official_symbol] = uniprot_gene_id_pubmed_list[
                            uniprot_entry.strip()]
                    else:
                        human_protein_official_symbol_pubmed_ids[human_protein_official_symbol] += uniprot_gene_id_pubmed_list[
                            uniprot_entry.strip()]
    print('human_protein_official_symbol_pubmed_ids: '+ str(len(human_protein_official_symbol_pubmed_ids.keys())))
    print('human_protein_official_symbols_temp: '+ str(len(human_protein_official_symbols_temp)))

    out_fp = dir_path + 'human_protein_official_symbol_pubmed_ids.csv'
    if not os.path.exists(out_fp):
        with open(out_fp, 'a') as output_file:
            for human_protein_official_symbol in human_protein_official_symbol_pubmed_ids.keys():
                pubmed_ids = list(human_protein_official_symbol_pubmed_ids[human_protein_official_symbol])
                while '' in pubmed_ids:
                    pubmed_ids.remove('')
                output_file.write(human_protein_official_symbol + '$' + ';'.join(list(pubmed_ids)) + '\n')

def collect_other_species_gene_related_pubmed_corpus(dir_path):

    # 获取其他物种对应的人类基因名称以及ID
    human_protein_official_symbol_ncbi_ids, human_protein_official_symbol_uniprot_entries, _ = prepare_human_protein_symbols_ids(dir_path)

    # 转换以ID为健，基因名称为值
    human_protein_ncbi_ids_official_symbol = {}
    human_protein_uniprot_entries_official_symbol = {}
    for key in human_protein_official_symbol_ncbi_ids.keys():
        if ';' not in human_protein_official_symbol_ncbi_ids[key]:
            ncbi_id = human_protein_official_symbol_ncbi_ids[key].strip()
            human_protein_ncbi_ids_official_symbol[ncbi_id] = key
        else:
            ncbi_ids = [x.strip() for x in human_protein_official_symbol_ncbi_ids[key].strip().split(';')]
            while '' in ncbi_ids:
                ncbi_ids.remove('')
            for ncbi_id in ncbi_ids:
                human_protein_ncbi_ids_official_symbol[ncbi_id] = key

    for key in human_protein_official_symbol_uniprot_entries.keys():
        if ';' not in human_protein_official_symbol_uniprot_entries[key]:
            uniprot_entry = human_protein_official_symbol_uniprot_entries[key].strip()
            human_protein_uniprot_entries_official_symbol[uniprot_entry] = key
        else:
            uniprot_entries = [x.strip() for x in human_protein_official_symbol_uniprot_entries[key].strip().split(';')]
            while '' in uniprot_entries:
                uniprot_entries.remove('')
            for uniprot_entry in uniprot_entries:
                human_protein_uniprot_entries_official_symbol[uniprot_entry] = key

    ncbi_id_out_fp = dir_path + 'other_species_ncbi_id_pubmed_ids.csv'
    uniprot_entry_out_fp = dir_path + 'other_species_uniprot_entry_pubmed_ids.csv'

    for species_name in species_dict.keys()[1:]:
        homologene_results = 'D:/1_data_help/text_corpus/homologene/results/'
        ncbi_homologous_gene_out_fp = homologene_results + 'human_homologous_gene_by_ncbi_id.csv'
        uniprot_homologous_gene_out_fp = homologene_results + 'human_homologous_gene_by_uniprot_entry.csv'

        ncbi_other_species_homologous_gene_id_human_gene_id = {}
        with codecs.open(ncbi_homologous_gene_out_fp, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if len(temp) < 3:
                    continue
                if temp[0].strip() != species_dict[species_name]:
                    continue
                if ';' not in temp[2].strip():
                    ncbi_other_species_homologous_gene_id_human_gene_id[temp[2].strip()] = temp[1].strip()
                else:
                    other_species_homologous_gene_ids = temp[2].strip().split(';')
                    for other_species_homologous_gene_id in other_species_homologous_gene_ids:
                        if other_species_homologous_gene_id.strip() == '':
                            continue
                        ncbi_other_species_homologous_gene_id_human_gene_id[other_species_homologous_gene_id.strip()] = temp[1].strip()

        uniprot_other_species_homologous_gene_entry_human_gene_entry = {}
        with codecs.open(uniprot_homologous_gene_out_fp, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('$')
                if len(temp) < 3:
                    continue
                if temp[0].strip() != species_dict[species_name]:
                    continue
                if ';' not in temp[2].strip():
                    uniprot_other_species_homologous_gene_entry_human_gene_entry[temp[2].strip()] = temp[1].strip()
                else:
                    other_species_homologous_gene_entries = temp[2].strip().split(';')
                    for other_species_homologous_gene_entry in other_species_homologous_gene_entries:
                        if other_species_homologous_gene_entry.strip() == '':
                            continue
                        uniprot_other_species_homologous_gene_entry_human_gene_entry[other_species_homologous_gene_entry.strip()] = \
                        temp[1].strip()

        ncbi_gene_id_pubmed_list, _ = get_gene_pubmed_list_from_ncbi(dir_path, species_name)
        uniprot_gene_id_pubmed_list, _ = get_gene_pubmed_list_from_uniprot(dir_path, species_name)

        with open(ncbi_id_out_fp, 'a') as output_file:
            for ncbi_gene_id in ncbi_gene_id_pubmed_list.keys():
                if ncbi_gene_id not in ncbi_other_species_homologous_gene_id_human_gene_id.keys():
                    continue
                human_gene_id = ncbi_other_species_homologous_gene_id_human_gene_id[ncbi_gene_id]
                if human_gene_id not in human_protein_ncbi_ids_official_symbol.keys():
                    continue
                human_protein_official_symbol = human_protein_ncbi_ids_official_symbol[human_gene_id]
                pubmed_ids = ncbi_gene_id_pubmed_list[ncbi_gene_id]
                while '' in pubmed_ids:
                    pubmed_ids.remove('')
                output_file.write(species_dict[species_name] + '$' + human_protein_official_symbol + '$' + \
                                  human_gene_id + '$' + ncbi_gene_id + '$' + ';'.join(pubmed_ids) + '\n')

        with open(uniprot_entry_out_fp, 'a') as output_file:
            for uniprot_gene_entry in uniprot_gene_id_pubmed_list.keys():
                if uniprot_gene_entry not in uniprot_other_species_homologous_gene_entry_human_gene_entry.keys():
                    continue
                human_gene_entry = uniprot_other_species_homologous_gene_entry_human_gene_entry[uniprot_gene_entry]
                if human_gene_entry not in human_protein_uniprot_entries_official_symbol.keys():
                    continue
                human_protein_official_symbol = human_protein_uniprot_entries_official_symbol[human_gene_entry]
                pubmed_ids = uniprot_gene_id_pubmed_list[uniprot_gene_entry]
                while '' in pubmed_ids:
                    pubmed_ids.remove('')
                output_file.write(species_dict[species_name] + '$' + human_protein_official_symbol + '$' + \
                                  human_gene_entry + '$' + uniprot_gene_entry + '$' + ';'.join(pubmed_ids) + '\n')

def collect_mgi_mouse_gene_phenotype_related_pubmed_corpus(dir_path):

    human_protein_official_symbol_ncbi_ids, _, _ = prepare_human_protein_symbols_ids(dir_path)
    human_protein_ncbi_ids_official_symbol = {}
    for key in human_protein_official_symbol_ncbi_ids.keys():
        if ';' not in human_protein_official_symbol_ncbi_ids[key]:
            ncbi_id = human_protein_official_symbol_ncbi_ids[key].strip()
            human_protein_ncbi_ids_official_symbol[ncbi_id] = key
        else:
            ncbi_ids = [x.strip() for x in human_protein_official_symbol_ncbi_ids[key].strip().split(';')]
            while '' in ncbi_ids:
                ncbi_ids.remove('')
            for ncbi_id in ncbi_ids:
                human_protein_ncbi_ids_official_symbol[ncbi_id] = key
    species_name = 'mouse'
    ncbi_gene_id_pubmed_list = get_gene_pubmed_list_from_mgi(dir_path)
    ncbi_homologous_gene_out_fp = dir_path + 'homologene/results/human_homologous_gene_by_ncbi_id.csv'
    ncbi_other_species_homologous_gene_id_human_gene_id = {}
    with codecs.open(ncbi_homologous_gene_out_fp, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if len(temp) < 3:
                continue
            if temp[0].strip() != species_dict[species_name]:
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

    ncbi_id_out_fp = dir_path + 'mouse_from_mgi_ncbi_id_pubmed_ids.csv'

    with open(ncbi_id_out_fp, 'a') as output_file:
        for ncbi_gene_id in ncbi_gene_id_pubmed_list.keys():
            if ncbi_gene_id not in ncbi_other_species_homologous_gene_id_human_gene_id.keys():
                continue
            human_gene_id = ncbi_other_species_homologous_gene_id_human_gene_id[ncbi_gene_id]
            if human_gene_id not in human_protein_ncbi_ids_official_symbol.keys():
                continue
            human_protein_official_symbol = human_protein_ncbi_ids_official_symbol[human_gene_id]
            pubmed_ids = ncbi_gene_id_pubmed_list[ncbi_gene_id]
            while '' in pubmed_ids:
                pubmed_ids.remove('')
            output_file.write(species_dict[species_name] + '$' + human_protein_official_symbol + '$' + \
                              human_gene_id + '$' + ncbi_gene_id + '$' + ';'.join(pubmed_ids) + '\n')

def collect_clinvar_human_gene_disease_related_pubmed_corpus(dir_path):

    human_protein_official_symbol_ncbi_ids, _, _ = prepare_human_protein_symbols_ids(dir_path)
    human_protein_ncbi_ids = set()
    for key in human_protein_official_symbol_ncbi_ids.keys():
        if ';' not in human_protein_official_symbol_ncbi_ids[key]:
            ncbi_id = human_protein_official_symbol_ncbi_ids[key].strip()
            human_protein_ncbi_ids.add(ncbi_id)
        else:
            ncbi_ids = [x.strip() for x in human_protein_official_symbol_ncbi_ids[key].strip().split(';')]
            while '' in ncbi_ids:
                ncbi_ids.remove('')
            for ncbi_id in ncbi_ids:
                human_protein_ncbi_ids.add(ncbi_id)
    print('human_protein_ncbi_ids: ' + str(len(human_protein_ncbi_ids)))

    allele_gene_id_ncbi_geneids = {}
    allele_gene_ids = set()
    mapping_allele_gene_id_ncbi_geneid_fname = dir_path + 'protein_symbol_pmid/disease_pubmed_raw_files/allele_gene.txt'
    with codecs.open(mapping_allele_gene_id_ncbi_geneid_fname, "rb", "utf-8") as input_file:
        for line in input_file:
            if line.strip()[0] == '#':
                continue
            temp = line.strip().split('	')
            if len(temp) < 2:
                continue
            if temp[0].strip() not in allele_gene_ids:
                allele_gene_id_ncbi_geneids[temp[0].strip()] = set()
                allele_gene_ids.add(temp[0].strip())
            allele_gene_id_ncbi_geneids[temp[0].strip()].add(temp[1].strip())
    print('allele_gene_id_ncbi_geneids: ' + str(len(allele_gene_id_ncbi_geneids.keys())))

    ncbi_geneid_pubmed_ids = {}
    ncbi_ids_from_clinvar = set()
    allele_gene_id_pmid_fname = dir_path + 'protein_symbol_pmid/disease_pubmed_raw_files/var_citations.txt'
    with codecs.open(allele_gene_id_pmid_fname, "rb", "utf-8") as input_file:
        for line in input_file:
            if line.strip()[0] == '#':
                continue
            temp = line.strip().split('	')
            if not temp[-1].strip().isdigit():
                continue
            if temp[0].strip() not in allele_gene_ids:
                continue
            geneids = allele_gene_id_ncbi_geneids[temp[0].strip()]
            for geneid in geneids:
                if geneid.strip() == '':
                    continue
                if geneid.strip() not in ncbi_ids_from_clinvar:
                    ncbi_geneid_pubmed_ids[geneid.strip()] = set()
                    ncbi_ids_from_clinvar.add(geneid.strip())
                ncbi_geneid_pubmed_ids[geneid.strip()].add(temp[-1].strip())
    print('ncbi_geneid_pubmed_ids: ' + str(len(ncbi_geneid_pubmed_ids.keys())))

    protein_symbol_pmids = {}
    protein_symbols = set()
    for ncbi_gene_id in ncbi_geneid_pubmed_ids.keys():
        pubmed_ids = ncbi_geneid_pubmed_ids[ncbi_gene_id]
        while '' in pubmed_ids:
            pubmed_ids.remove('')
        for protein_symbol in human_protein_official_symbol_ncbi_ids.keys():
            if protein_symbol.strip() == '':
                continue
            if ncbi_gene_id in human_protein_official_symbol_ncbi_ids[protein_symbol]:
                if protein_symbol not in protein_symbols:
                    protein_symbol_pmids[protein_symbol] = set()
                    protein_symbols.add(protein_symbol)
                protein_symbol_pmids[protein_symbol] = protein_symbol_pmids[protein_symbol] | pubmed_ids
    print('protein_symbol_pmids: ' + str(len(protein_symbol_pmids.keys())))

    ncbi_id_out_fp = dir_path + 'protein_symbol_pmid/human_geneid_from_clinvar_disease_pubmed_ids_supplement.csv'

    with open(ncbi_id_out_fp, 'a') as output_file:
        for protein_symbol in protein_symbol_pmids.keys():
            pubmed_ids = list(protein_symbol_pmids[protein_symbol])
            output_file.write(protein_symbol + '$' + ';'.join(pubmed_ids) + '\n')


'''
    3、Collect relevant literature
'''

# Generating pmid set
def generate_pubmed_ids_set(dir_path):

    file_name_feature = '_pubmed_ids.csv'
    raw_data_dir = dir_path + 'raw_data/'
    file_names = [f for f in listdir(raw_data_dir) if f.endswith(file_name_feature)]
    print('The number of the file_names: '+str(len(file_names)))
    pubmed_ids_all = set()
    for file_name in file_names:
        raw_data_fp = raw_data_dir + file_name
        with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                if temp[-1].strip() == '':
                    continue
                pubmed_ids_line = temp[-1].strip().split(';')
                for pubmed_id in pubmed_ids_line:
                    if not pubmed_id.strip().isdigit():
                        continue
                    pubmed_ids_all.add(pubmed_id.strip())
    print('pubmed_ids_all: '+str(len(pubmed_ids_all)))

    # raw_data_fp = dir_path + 'mouse_from_mgi_ncbi_id_pubmed_ids.csv'
    # pubmed_ids_supplement = set()
    # with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
    #     for line in islice(input_file.readlines(), 0, None):
    #         temp = line.strip().split('$')
    #         if temp[-1].strip() == '':
    #             continue
    #         pubmed_ids_line = temp[-1].strip().split(';')
    #         for pubmed_id in pubmed_ids_line:
    #             if not pubmed_id.strip().isdigit():
    #                 continue
    #             pubmed_ids_supplement.add(pubmed_id.strip())
    # print('pubmed_ids_supplement: ' + str(len(pubmed_ids_supplement)))

    # raw_data_fp = dir_path + 'go_literature_reference_pubmed_ids_from_Pfam_supplement.csv'
    # pubmed_ids_supplement = set()
    # with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
    #     for line in islice(input_file.readlines(), 0, None):
    #         temp = line.strip().split('	')
    #         if len(temp) < 3:
    #             continue
    #         pubmed_ids_supplement.add(temp[1].strip())
    # print('pubmed_ids_supplement: ' + str(len(pubmed_ids_supplement)))

    # raw_data_fp = dir_path + 'gene2go'
    # pubmed_ids_supplement = set()
    # with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
    #     for line in input_file:
    #         if line.strip()[0] == '#':
    #             continue
    #         temp = line.strip().split('	')
    #         if len(temp) < 7:
    #             continue
    #         if temp[6].strip() == '-':
    #             continue
    #         pubmed_ids_supplement= pubmed_ids_supplement | set(temp[6].strip().split('|'))
    # print('pubmed_ids_supplement: ' + str(len(pubmed_ids_supplement)))

    # 解析go是出现bug ,漏掉以“|”分割的pmids
    # raw_data_fp = dir_path + 'go_literature_reference_pubmed_ids_supplement.csv'
    # pubmed_ids_supplement = set()
    # with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
    #     for line in input_file:
    #         temp = line.strip().split("$")
    #         if len(temp) < 2:
    #             continue
    #         if ';' not in temp[1].strip():
    #             continue
    #         pubmed_ids_supplement= pubmed_ids_supplement | set(temp[1].strip().split(';'))
    # print('pubmed_ids_supplement: ' + str(len(pubmed_ids_supplement)))

    raw_data_fp = dir_path + 'human_geneid_from_clinvar_disease_pubmed_ids.csv'
    pubmed_ids_supplement = set()
    with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if len(temp) < 2:
                continue
            if temp[1].strip() == '':
                continue
            if ';' not in temp[1].strip():
                pubmed_ids_supplement.add(temp[1].strip())
            else:
                for pmid in temp[1].strip().split(';'):
                    if pmid.strip() == '':
                        continue
                    pubmed_ids_supplement.add(pmid.strip())
    print('pubmed_ids_supplement: ' + str(len(pubmed_ids_supplement)))

    pubmed_ids_supplement_difference = pubmed_ids_supplement.difference(pubmed_ids_all)
    print('pubmed_ids_supplement_difference: ' + str(len(pubmed_ids_supplement_difference)))

    sub_list_size = 60
    pubmed_ids_all_list = list(pubmed_ids_supplement_difference)
    pubmed_ids_all_list_chunks = [pubmed_ids_all_list[i:i + sub_list_size] for i in range(0, len(pubmed_ids_all_list), sub_list_size)]
    print('pubmed_ids_all_list_chunks: ' + str(len(pubmed_ids_all_list_chunks)))

    return pubmed_ids_all_list_chunks

# If there is a problem with batch collection, single collection
def efetch_literatures_by_individual(pmid_list,dir_path,download_fn):

    # 每次fetch一条文献
    file_name = dir_path + download_fn + ".txt"
    with open(file_name, 'a') as output_file:
        for pmid in pmid_list:
            if not os.path.isfile(file_name):
                # Downloading...
                net_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
                try:
                    output_file.write(net_handle.read())
                except AttributeError as e:
                    print(e)
                    continue

# Standardized literature data were obtained according to PMID
def efetch_literatures_by_batch(dir_path):

    # 按蛋白质名称读入所有PMIDs

    pubmed_ids_all_list_chunks = generate_pubmed_ids_set(dir_path)
    download_fp = dir_path + "download_literatures_supplement/"
    if not os.path.exists(download_fp):
        os.mkdir(download_fp)
    completed_file_names = [f for f in listdir(download_fp) if f.endswith('txt')]
    completed_serial_numbers = set()
    for completed_file_name in completed_file_names:
        serial_number = re.sub('\.txt','',completed_file_name)
        completed_serial_numbers.add(serial_number)

    for i in range(len(pubmed_ids_all_list_chunks)):
        pmid_list = pubmed_ids_all_list_chunks[i]
        download_fn = str(i+len(completed_file_names))
        if download_fn in completed_serial_numbers:
            continue
        if len(pmid_list) == 0:
            continue
        print(download_fn)
        # 每次fetch一批文献
        search_results = Entrez.epost("pubmed", id=",".join(pmid_list)).read()
        soup = BeautifulSoup(search_results, 'lxml')
        try:
            webenv = soup.find("webenv").text.strip()
            query_key = soup.find("querykey").text.strip()
        except AttributeError as e:
            efetch_literatures_by_individual(pmid_list, dir_path,download_fn)
            print(e)
            continue
        count = len(pmid_list)
        print "Found %i results" % count
        batch_size = 10
        out_fp = download_fp + download_fn + ".txt"
        out_handle = open(out_fp, "w")
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            try:
                fetch_handle = Entrez.efetch(db="pubmed",
                                             rettype="medline", retmode="text",
                                             retstart=start, retmax=batch_size,
                                             webenv=webenv,
                                             query_key=query_key
                                            )
            except urllib2.URLError as err:
                print err
                continue
            except socket.error:
                errno, errstr = sys.exc_info()[:2]
                if errno == socket.timeout:
                    print "There was a timeout"
                else:
                    print "There was some other socket error"
                continue
            print "Going to download record %i to %i" % (start + 1, end)
            try:
                data = fetch_handle.read()
            except ssl.SSLError as e:
                print e
                continue
            fetch_handle.close()
            out_handle.write(data)
        out_handle.close()

# Parse literature standardization data
def parse_file(filePath):

    # 解析已下载的文件
    download_fp = filePath + "download_literatures_supplement/"
    parse_fp = filePath + 'parse_literatures_supplement/'
    if not os.path.exists(parse_fp):
        os.mkdir(parse_fp)
    download_file_names = [f for f in listdir(download_fp) if f.endswith('.txt')]
    print 'download_file_names num: ' + str(len(download_file_names))

    completed_file_names = [f for f in listdir(parse_fp) if f.endswith('.txt')]
    print 'completed_file_names num: ' + str(len(completed_file_names))

    for download_file_name in download_file_names:
        if download_file_name in completed_file_names:
            continue
        input = open(download_fp + download_file_name)

        # 一次解析多条文献数据
        records = Medline.parse(input)
        records = list(records)
        with open(parse_fp + download_file_name, 'a+') as output_file:
            for record in records:
                parse_results = []
                pmid = record.get("PMID", "?")
                parse_results.append(pmid)
                date = record.get("DCOM", "?")
                parse_results.append(date)
                title = record.get("TI", "?")
                parse_results.append(re.sub('\n','',title))
                abstract = record.get("AB", "?")
                parse_results.append(re.sub('\n','',abstract))
                country = record.get("PL", "?")
                parse_results.append(re.sub('\n', '', country))
                journal = record.get("JT", "?")
                parse_results.append(re.sub('\n', '', journal))
                topic = record.get("MH", "?")
                parse_results.append(re.sub('\n', '', ';;'.join(topic)))
                output_file.write('$'.join(parse_results) + '\n')

# Perform literature collection procedures
def exec_fetch_parse(filePath):
    efetch_literatures_by_batch(filePath)
    parse_file(filePath)

'''
    4、Prepare the protein definition text corpus
'''

# Summarize data from three databases
def prepare_protein_summary_definition_text(dir_path):
    human_protein_official_symbol_ncbi_ids, human_protein_official_symbol_uniprot_entries,train_protein_official_symbols = prepare_human_protein_symbols_ids(dir_path)

    # human_protein_official_symbols_set = set(human_protein_official_symbol_ncbi_ids.keys()).union(human_protein_official_symbol_uniprot_entries.keys())

    raw_data_fp = dir_path + 'protein_summary/protein_summary_from_genecards.csv'
    human_protein_official_symbol_text = {}
    human_protein_official_symbols = set()
    with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) < 2:
                continue
            if temp[1].strip() == '':
                continue
            text_clean = ' '.join([x.strip() for x in '$'.join(temp[1:]).strip().split(' ')])
            if temp[0].strip() in human_protein_official_symbols:
                if len(text_clean) > len(human_protein_official_symbol_text[temp[0].strip()]):
                    human_protein_official_symbol_text[temp[0].strip()] = text_clean
            else:
                human_protein_official_symbol_text[temp[0].strip()] = text_clean
    print('human_protein_official_symbol_text: ' + str(len(human_protein_official_symbol_text.keys())))

    proteins_summary_mapped = set()
    # 注意并不用相应的基因名称进行替换，只是用基因别称对应的定义
    special_processing_proteins_for_train = {'TSBP1':'C6orf10','HEXD':'HEXDC','DOP1B':'DOPEY2','PRXL2B':'FAM213B','PRXL2A':'FAM213A','SHOC1':'C9orf84'}
    skip_proteins_for_train = ['PRED58']
    for protein_official_symbol in train_protein_official_symbols.keys():
        if protein_official_symbol in skip_proteins_for_train:
            continue
        if protein_official_symbol in human_protein_official_symbol_text.keys():
            proteins_summary_mapped.add(protein_official_symbol)
        else:
            for protein in train_protein_official_symbols[protein_official_symbol]:
                if protein in human_protein_official_symbol_text.keys():
                    if protein in special_processing_proteins_for_train.values():
                        if protein in proteins_summary_mapped:
                            print(protein_official_symbol,protein)
                        proteins_summary_mapped.add(protein_official_symbol)
    print('proteins_summary_mapped: ' + str(len(proteins_summary_mapped)))

    # print(set(proteins_summary_mapped).difference(human_protein_official_symbols_set))

    out_fp = dir_path + 'protein_summary/human_protein_official_symbol_in_train_dataset_text_clean.csv'
    if not os.path.exists(out_fp):
        with open(out_fp, 'a') as output_file:
            for protein in proteins_summary_mapped:
                if protein in special_processing_proteins_for_train.keys():
                    protein_text = human_protein_official_symbol_text[special_processing_proteins_for_train[protein]]
                else:
                    protein_text = human_protein_official_symbol_text[protein]
                protein_text = ' '.join(protein_text.split('$'))
                protein_text = re.sub('      ',' ',protein_text)
                output_file.write(protein + '$' + protein_text + '\n')

# Get the mesh for each gene/protein
def extract_literature_mesh_words_by_pmid(dir_path):
    parse_data_dir = dir_path + 'parse_literatures'
    parse_files_path = parse_data_dir + '/'
    for i in range(2):
        if i == 1:
            parse_files_path = parse_data_dir + '_supplement' + '/'
        parse_file_names = [f for f in listdir(parse_files_path) if f.endswith('txt')]
        pubmed_id_mesh_words = {}
        for file_name in parse_file_names:
            file_path = parse_files_path + file_name
            with codecs.open(file_path, 'rb', 'utf-8') as input_file:
                for line in input_file:
                    temp = line.strip().split('$')
                    if (temp[-1].strip() == '') or (temp[-1].strip() == '?'):
                        continue
                    if temp[0].strip().isdigit() is False:
                        continue
                    pubmed_id_mesh_words[temp[0].strip()] = temp[-1].strip()

        write_out_fp = dir_path + 'pubmed_id_mesh_words.csv'
        # if not os.path.exists(write_out_fp):
        with open(write_out_fp, 'a') as output_file:
            for pubmed_id in pubmed_id_mesh_words.keys():
                output_file.write(pubmed_id + '$' + pubmed_id_mesh_words[pubmed_id] + '\n')

# Get the title and abstract of each pmid
def extract_literature_date_title_abstract_by_pmid(dir_path):
    parse_data_dir = dir_path + 'parse_literatures'
    parse_files_path = parse_data_dir + '/'
    for i in range(2):
        if i == 1:
            parse_files_path = parse_data_dir + '_supplement' + '/'
        parse_file_names = [f for f in listdir(parse_files_path) if f.endswith('txt')]
        pubmed_id_date_title_abstract_journal = {}
        for file_name in parse_file_names:
            file_path = parse_files_path + file_name
            with codecs.open(file_path, 'rb', 'utf-8') as input_file:
                for line in input_file:
                    temp = line.strip().split('$')
                    if temp[0].strip().isdigit() is False:
                        continue
                    pubmed_id_date_title_abstract_journal[temp[0].strip()] = '$'.join(temp[1:-1])

        write_out_fp = dir_path + 'pubmed_id_date_title_abstract_journal.csv'
        # if not os.path.exists(write_out_fp):
        with open(write_out_fp, 'a') as output_file:
            for pubmed_id in pubmed_id_date_title_abstract_journal.keys():
                output_file.write(pubmed_id + '$' + pubmed_id_date_title_abstract_journal[pubmed_id] + '\n')

# The PMID set corresponding to each protein is obtained
def obtain_protein_pmid_set(dir_path):
    file_name_feature = 'csv'
    raw_data_dir = dir_path + 'protein_symbol_pmid/'
    file_names = [f for f in listdir(raw_data_dir) if f.endswith(file_name_feature)]
    print('The number of the file_names: ' + str(len(file_names)))
    human_protein_symbol_pmids = {}
    human_protein_symbols = set()
    for file_name in file_names:
        if 'go_' in file_name:
            continue
        raw_data_fp = raw_data_dir + file_name
        with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                if temp[-1].strip() == '':
                    continue
                pubmed_ids_line = temp[-1].strip().split(';')
                while '' in pubmed_ids_line:
                    pubmed_ids_line.remove('')
                if 'human' in file_name:
                    if temp[0].strip() not in human_protein_symbols:
                        human_protein_symbol_pmids[temp[0].strip()] = []
                        human_protein_symbols.add(temp[0].strip())
                    human_protein_symbol_pmids[temp[0].strip()] += pubmed_ids_line
                else:
                    if temp[1].strip() not in human_protein_symbols:
                        human_protein_symbol_pmids[temp[1].strip()] = []
                        human_protein_symbols.add(temp[1].strip())
                    human_protein_symbol_pmids[temp[1].strip()] += pubmed_ids_line
    print('human_protein_symbol_pmids: ' + str(len(human_protein_symbol_pmids.keys())))

    return human_protein_symbol_pmids

# MeSH terms are processed and matched by protein
def processing_mesh_words_by_protein(dir_path):

    input_fp = dir_path + 'extract_literatures/pubmed_id_mesh_words.csv'
    pubmed_id_mesh_words_set = {}
    with codecs.open(input_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if len(temp) < 2:
                continue
            if temp[1].strip() == '?':
                continue
            pubmed_id_mesh_words_set[temp[0].strip()] = temp[1].strip().split(';;')
    print('pubmed_id_mesh_words_set: ' + str(len(pubmed_id_mesh_words_set.keys())))
    pubmed_id_set = set(pubmed_id_mesh_words_set.keys())

    human_protein_symbol_pmids = obtain_protein_pmid_set(dir_path)

    _, _, train_protein_official_symbols = prepare_human_protein_symbols_ids(dir_path)

    protein_official_symbols_for_train_mesh = set(train_protein_official_symbols.keys()).intersection(set(human_protein_symbol_pmids.keys()))
    protein_official_symbol_pmids_mesh_words_for_train = {}
    for protein_official_symbol in protein_official_symbols_for_train_mesh:
        pubmed_ids = human_protein_symbol_pmids[protein_official_symbol]
        temp_mesh_words =set()
        for pubmed_id in pubmed_ids:
            if pubmed_id in pubmed_id_set:
                for mesh_word in pubmed_id_mesh_words_set[pubmed_id]:
                    temp_mesh_words.add(mesh_word.strip())
        protein_official_symbol_pmids_mesh_words_for_train[protein_official_symbol] = temp_mesh_words

    print('protein_official_symbol_pmids_mesh_words_for_train: ' + str(len(protein_official_symbol_pmids_mesh_words_for_train.keys())))

    mesh_words_phrases = set()
    for key in protein_official_symbol_pmids_mesh_words_for_train.keys():
        for phrase in protein_official_symbol_pmids_mesh_words_for_train[key]:
            mesh_words_phrases.add(phrase)
    print('mesh_words_phrases: ' + str(
        len(mesh_words_phrases)))

    write_out_fp = dir_path + 'MeSH_terms/protein_official_symbol_pmids_mesh_words_for_train_raw_mesh_terms.csv'
    if not os.path.exists(write_out_fp):
        with open(write_out_fp, 'w') as output_file:
            for protein_official_symbol in protein_official_symbol_pmids_mesh_words_for_train.keys():
                mesh_words = ';;'.join(list(protein_official_symbol_pmids_mesh_words_for_train[protein_official_symbol]))
                output_file.write(protein_official_symbol + '$' + mesh_words + '\n')

    write_out_fp = dir_path + 'MeSH_terms/protein_official_symbol_pmids_mesh_words_for_train_modify_mesh_terms.csv'
    if not os.path.exists(write_out_fp):
        with open(write_out_fp, 'w') as output_file:
            for protein_official_symbol in protein_official_symbol_pmids_mesh_words_for_train.keys():
                mesh_words = set()
                for terms in protein_official_symbol_pmids_mesh_words_for_train[protein_official_symbol]:
                    temp_terms = re.sub('\*','',terms.strip())
                    temp_terms = re.sub('/', ' ', temp_terms.strip())
                    mesh_words.add(temp_terms.strip())
                mesh_words_str = ' '.join(list(mesh_words))
                output_file.write(protein_official_symbol + '$' + mesh_words_str + '\n')

# Extract the title and summary
def processing_title_abstract(dir_path):

    pmid_title_abstract = {}
    pmids = set()
    pubmed_data_fname = dir_path + 'extract_literatures/pubmed_id_date_title_abstract_journal.csv'
    with codecs.open(pubmed_data_fname, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if len(temp) < 6:
                continue
            title_abstrct = ''.join(temp[2:4])
            if (temp[0].strip() in pmids) or (title_abstrct == '??'):
                continue
            pmid_title_abstract[temp[0].strip()] = title_abstrct
            pmids.add(temp[0].strip())
    print('pmid_title_abstract: '+str(len(pmid_title_abstract.keys())))

    human_protein_symbol_pmids = obtain_protein_pmid_set(dir_path)
    write_out_fp =  dir_path + 'extract_literatures/corpus_literatures_title_abstract_by_protein.csv'
    with open(write_out_fp, 'w') as output_file:
        for human_protein_symbol in human_protein_symbol_pmids.keys():
            if len(human_protein_symbol_pmids[human_protein_symbol]) == 0:
                continue
            temp_pmids = set(human_protein_symbol_pmids[human_protein_symbol])
            for pmid in temp_pmids:
                if pmid.strip() == '':
                    continue
                if pmid.strip() not in pmids:
                    continue
                output_file.write(pmid.strip()+'$'+pmid_title_abstract[pmid.strip()] + '\n')

# Read the MeSH terms interpretation
def parse_xml_for_mesh_terms_info(dir_path):

    mesh_dir_path = dir_path + 'MeSH_terms/desc2019.xml'
    handle = open(mesh_dir_path)
    records = Entrez.parse(handle)
    write_out_fp = dir_path + 'MeSH_terms/parse_MeSH_terms_info.csv'
    with open(write_out_fp, 'w') as output_file:
        for record in records:
            mesh_info = []
            mesh_term_name = ''
            mesh_term_id = ''
            mesh_term_concept = ''
            mesh_term_tree = ''
            try:
                mesh_term_name = re.sub('\n','',record['ConceptList'][0]['ConceptName']['String'].strip()).strip()
            except KeyError as e:
                print e
                continue
            try:
                mesh_term_id = re.sub('\n','',record['DescriptorUI'].strip()).strip()
            except KeyError as e:
                print e
                continue
            try:
                mesh_term_concept = re.sub('\n', '', record['ConceptList'][0]['ScopeNote'].strip()).strip()
            except KeyError as e:
                print e
            try:
                mesh_term_tree = ';'.join([x.strip() for x in record['TreeNumberList']])
            except KeyError as e:
                print e
            mesh_info.append(mesh_term_name)
            mesh_info.append(mesh_term_id)
            mesh_info.append(mesh_term_concept)
            mesh_info.append(mesh_term_tree)
            output_file.write('$'.join(mesh_info)+ '\n')


def temp_read_write_function(dir_path):
    read_in_fp = dir_path + 'protein_symbol_pmid/go_pubmed_raw_files/gene2go_pubmed_ids'
    write_out_fp = dir_path + 'protein_symbol_pmid/go_literature_reference_pubmed_ids_supplement.csv'
    with open(write_out_fp, 'w') as output_file:
        with codecs.open(read_in_fp, 'rb', 'utf-8') as input_file:
            for line in input_file:
                if line.strip()[0] == '#':
                    continue
                temp = line.strip().split('	')
                if len(temp) < 7:
                    continue
                if temp[6].strip() == '-':
                    continue
                pmids = temp[6].strip().split('|')
                output_file.write(temp[1].strip()+'$'+';'.join(pmids)+'\n')



if __name__ == '__main__':

    dir_path = ''
    collect_human_gene_related_pubmed_corpus(dir_path)
    collect_other_species_gene_related_pubmed_corpus(dir_path)
    collect_mgi_mouse_gene_phenotype_related_pubmed_corpus(dir_path)
    collect_clinvar_human_gene_disease_related_pubmed_corpus(dir_path)
    prepare_human_protein_symbols_ids(dir_path)
    generate_pubmed_ids_set(dir_path)
    exec_fetch_parse(dir_path)
    prepare_protein_summary_definition_text(dir_path)
    extract_literature_mesh_words_by_pmid(dir_path)
    extract_literature_date_title_abstract_by_pmid(dir_path)
    processing_mesh_words_by_protein(dir_path)
    processing_title_abstract(dir_path)
    parse_xml_for_mesh_terms_info(dir_path)
    temp_read_write_function(dir_path)




