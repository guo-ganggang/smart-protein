#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 22/11/2018 3:36 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : data_processing.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

from collections import OrderedDict
import gzip
import codecs,os
from itertools import islice
import scipy.stats as stats
import re
from sklearn.decomposition import PCA
import numpy as np

'''
    Data were downloaded from GTEX database, that is, 
    mRNA expression amount in each human tissue, 
    and relevant protein expression data were screened out
'''


# Spearman rank correlation coefficient was calculated
def spearman_rank(X,Y):
    n = len(X)
    Xrank = stats.rankdata(X, method='average')
    # N minus 2 is the penultimate element
    Yrank = stats.rankdata(Y, method='average')
    # This is a place where the order of Xrank and Yrank doesn't matter, because you're squaring them
    diffs = Xrank - Yrank
    r_s = 1 - 6*sum(diffs*diffs)/(n*(n**2 - 1))

    return r_s

# alculate the median value of the list
def get_median(data):
    data = sorted(data)
    size = len(data)
    if size % 2 == 0: # 判断列表长度为偶数
        median = (data[size//2]+data[size//2-1])/2
        data[0] = median
    if size % 2 == 1: # 判断列表长度为奇数
        median = data[(size-1)//2]
        data[0] = median
    return data[0]

# Read in a list of human protein names and their abbreviations,
# and at the same time, read in a list of protein names
# and their abbreviations related to the training set
def prepare_human_protein_symbols_ids():
    dir_path = '/text_corpus/'
    raw_data_dir = '/protein_entities_standardized/'
    raw_data_fp = raw_data_dir + 'human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    human_protein_official_symbol_names = {}
    del_train_protein_symbols = ['FAM58BP']
    with codecs.open(raw_data_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if temp[0].strip() in del_train_protein_symbols:
                continue
            human_protein_official_symbol_names[temp[0].strip()] = set(temp[1:])
    print('human_protein_official_symbol_names: ' + str(len(human_protein_official_symbol_names.keys())))
    human_protein_official_symbols = set(human_protein_official_symbol_names.keys())

    train_data_fp = dir_path + 'prepare_gene_symbols_abbreviation_for_train_model.csv'
    train_protein_official_symbols_aliases = {}
    special_train_protein_symbols = {'NAT6': 'NAA80', 'EFTUD1': 'EFL1', 'LNP': 'LNPK', 'APITD1': 'CENPS-CORT'}

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


# The sample number corresponding to each organization is counted
def prepare_tissue_sample_ids(dir_path):
    gtex_expression_tissue_fp = dir_path + 'GTEx/GTEx_v7_Annotations_SampleAttributesDS.txt'

    flag_terms = ['Blood	Whole Blood', 'Brain	Brain', 'Lung	Lung', 'Muscle	Muscle', 'Bone Marrow	Cells',
                  'Heart	Heart', 'Skin	Skin', 'Nerve	Nerve', 'Thyroid	Thyroid', 'Esophagus	Esophagus',
                  'Adipose Tissue', 'Blood Vessel', 'Blood	Cells', 'Pituitary	Pituitary', 'Testis	Testis',
                  'Pancreas	Pancreas', 'Prostate	Prostate', 'Skin	Cells']

    tissue_sample_ids = OrderedDict()
    sample_attribute_ids = set()
    with codecs.open(gtex_expression_tissue_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split('	')
            if len(temp) == 63:
                tissue = temp[6].strip()
                if tissue not in tissue_sample_ids:
                    tissue_sample_ids[tissue] = set()
                tissue_sample_ids[tissue].add(temp[0].strip())
                sample_attribute_ids.add(temp[0].strip())
            else:
                for flag_term in flag_terms:
                    if flag_term in line.strip():
                        tissue = ''
                        for i in range(len(temp)):
                            if temp[i].strip() == flag_term.split('	')[0].strip():
                                tissue = temp[i + 1].strip()
                                break
                        if tissue not in tissue_sample_ids:
                            tissue_sample_ids[tissue] = set()
                        tissue_sample_ids[tissue].add(temp[0].strip())
                        sample_attribute_ids.add(temp[0].strip())
                        break
    print(len(tissue_sample_ids.keys()))
    print(list(tissue_sample_ids.keys())[-1])
    print(len(sample_attribute_ids))

    return tissue_sample_ids

# Median value was used to obtain tissue data related to each gene
def processing_GTEx_expression_data_base_median(dir_path,first_complete):

    tissue_sample_ids = prepare_tissue_sample_ids(dir_path)

    # for key in tissue_sample_ids.keys():
    #     print(key,str(len(tissue_sample_ids[key])))

    # gene_symbols_fp = dir_path + 'gene_symbols_mgi_links.csv'
    # gene_symbols = []
    # with codecs.open(gene_symbols_fp, "rb", "utf-8") as input_file:
    #     for line in islice(input_file.readlines(), 0, None):
    #         temp = line.strip().split(',')
    #         gene_symbols.append(temp[0].strip())
    # print(len(gene_symbols))

    _,_,train_protein_official_symbols_aliases,_ = prepare_human_protein_symbols_ids()

    total_relative_protein_entries = set()
    for key in train_protein_official_symbols_aliases.keys():
        total_relative_protein_entries.add(key)
        for protein_entry in train_protein_official_symbols_aliases[key]:
            total_relative_protein_entries.add(protein_entry)
    print('total_relative_protein_entries: '+str(len(total_relative_protein_entries)))

    gtex_expression_values_fp = dir_path + 'GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz'
    # gtex_expression_values_fp = dir_path + 'GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz'
    skip_lines = 0
    sample_ids = []
    gene_name_values = {}
    with gzip.open(gtex_expression_values_fp, 'rb', 'utf-8') as gzipped_file:
        for line in gzipped_file:
            if skip_lines < 2:
                skip_lines += 1
                continue
            temp = line.strip().split('	')
            if skip_lines == 2:
                for sample_id in temp[2:]:
                    sample_ids.append(sample_id.strip())
                skip_lines += 1
                print('sample_ids: '+str(len(sample_ids)))
                continue
            temp = line.strip().split('	')
            if (temp[1].strip() in total_relative_protein_entries) or (temp[1].strip() in first_complete):
                gene_name_values[temp[1].strip()] = temp[2:]
    print('gene_name_values: '+str(len(gene_name_values.keys())))

    gene_symbol_vector = {}
    gene_expression_in_all_tissues_fp = dir_path + 'gene_expression_in_all_tissues_by_tpm_median.csv'

    complete_gene_symbols = set()
    if os.path.isfile(gene_expression_in_all_tissues_fp):
        with codecs.open(gene_expression_in_all_tissues_fp, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 1, None):
                temp = line.strip().split(',')
                if len(temp) != len(tissue_sample_ids.keys()):
                    continue
                complete_gene_symbols.add(temp[0].strip())
        print('Complete gene symbols: '+str(len(complete_gene_symbols)))

    with open(gene_expression_in_all_tissues_fp, 'a') as output_file:
        if os.path.getsize(gene_expression_in_all_tissues_fp) == 0:
            output_file.write('gene_symbol' + ',' + ','.join(tissue_sample_ids.keys()[:-1])+'\n')
        for gene_name in gene_name_values.keys():
            if gene_name not in first_complete:
                continue
            if gene_name in complete_gene_symbols:
                continue
            tissue_sample_average_values = []
            for tissue in tissue_sample_ids.keys()[:-1]:
                temp_value = []
                for sample_id in tissue_sample_ids[tissue]:
                    if sample_id not in sample_ids:
                        # print(sample_id)
                        continue
                    temp_value.append(float(gene_name_values[gene_name][sample_ids.index(sample_id)]))
                tissue_sample_average_values.append(get_median(temp_value))
            gene_symbol_vector[gene_name] = tissue_sample_average_values
            output_file.write(gene_name + ',' + ','.join([str(x) for x in tissue_sample_average_values])+'\n')
            print(gene_name + ',' + ','.join([str(x) for x in tissue_sample_average_values]))
    print(len(gene_symbol_vector.keys()))

# Median value was used to obtain tissue data related to each gene
def processing_GTEx_expression_data(dir_path):

    tissue_sample_ids = prepare_tissue_sample_ids(dir_path)

    _,_,train_protein_official_symbols_aliases, human_protein_official_symbol_names = prepare_human_protein_symbols_ids()

    total_human_protein_entries = set()
    for key in human_protein_official_symbol_names.keys():
        total_human_protein_entries.add(key.strip())
        for aliase in human_protein_official_symbol_names[key]:
            total_human_protein_entries.add(aliase.strip())

    print('total_human_protein_entries: '+str(len(total_human_protein_entries)))

    gtex_expression_values_fp = dir_path + 'GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz'
    # gtex_expression_values_fp = dir_path + 'GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz'
    skip_lines = 0
    sample_ids = []
    gene_name_values = {}
    with gzip.open(gtex_expression_values_fp, 'rb', 'utf-8') as gzipped_file:
        for line in gzipped_file:
            if skip_lines < 2:
                skip_lines += 1
                continue
            temp = line.strip().split('	')
            if skip_lines == 2:
                for sample_id in temp[2:]:
                    sample_ids.append(sample_id.strip())
                skip_lines += 1
                print('sample_ids: '+str(len(sample_ids)))
                continue
            temp = line.strip().split('	')
            if temp[1].strip() in total_human_protein_entries:
                gene_name_values[temp[1].strip()] = temp[2:]
    print('gene_name_values: '+str(len(gene_name_values.keys())))
    gene_names = list(gene_name_values.keys())

    # gene_symbol_vector = {}
    gene_expression_in_all_tissues_fp = dir_path + 'gene_expression_in_all_tissues_by_total_tpm_values.csv'
    if not os.path.exists(gene_expression_in_all_tissues_fp):

        with open(gene_expression_in_all_tissues_fp, 'w') as output_file:
            output_file.write('gene_symbol' + '$' + '$'.join(tissue_sample_ids.keys()[:-1]) + '\n')

        # 按组织匹配样本索引号
        tissue_sample_indexes = []
        for tissue in tissue_sample_ids.keys()[:-1]:
            temp_indexes = []
            for sample_id in tissue_sample_ids[tissue]:
                if sample_id not in sample_ids:
                    # print(sample_id)
                    continue
                temp_indexes.append(sample_ids.index(sample_id))
            tissue_sample_indexes.append(temp_indexes)

        tissue_sample_values = []
        for tissue_sample_index in tissue_sample_indexes:
            tissue_sample_value = []
            for gene_name in gene_names:
                gene_tissue_sample_values = gene_name_values[gene_name]
                temp_values = []
                for sample_index in tissue_sample_index:
                    temp_values.append(gene_tissue_sample_values[sample_index])
                tissue_sample_value.append(map(float, temp_values))
            tissue_sample_value_pca = PCA(n_components=5).fit_transform(np.array(tissue_sample_value, dtype=np.float32).reshape(len(gene_names),len(tissue_sample_index)))
            tissue_sample_values.append(tissue_sample_value_pca.tolist())
        print('tissue_sample_values: ' + str(len(tissue_sample_values)))

        with open(gene_expression_in_all_tissues_fp, 'a') as output_file:
            for i in range(len(gene_names)):
                gene_name = gene_names[i]
                temp_tissue_sample_value = []
                for tissue_sample_value in tissue_sample_values:
                    temp_tissue_sample_value.append(','.join((map(str, tissue_sample_value[i]))))
                output_file.write(gene_name + '$' + '$'.join(temp_tissue_sample_value) +'\n')

    human_gene_symbol_tissues_expression_values = {}
    human_gene_symbols = set()
    with codecs.open(gene_expression_in_all_tissues_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            expression_values = ','.join(temp[1:])
            human_gene_symbol_tissues_expression_values[temp[0].strip()] = expression_values
            human_gene_symbols.add(temp[0].strip())
    print('human_gene_symbol_tissues_expression_values: ' + str(len(human_gene_symbol_tissues_expression_values.keys())))

    train_gene_expression_in_all_tissues_fp = dir_path + 'train_gene_expression_in_all_tissues_by_total_tpm_values.csv'
    with open(train_gene_expression_in_all_tissues_fp, 'a') as output_file:
        for protein_symbol in train_protein_official_symbols_aliases.keys():
            if protein_symbol in human_gene_symbols:
                output_file.write(protein_symbol + '$' + human_gene_symbol_tissues_expression_values[protein_symbol] + '\n')
            else:
                for protein in train_protein_official_symbols_aliases[protein_symbol]:
                    if protein.strip() in human_gene_symbols:
                        output_file.write(
                            protein_symbol + '$' + human_gene_symbol_tissues_expression_values[protein.strip()] + '\n')
                        break

# The protein names that need to be collected are screened out
def supplementary_GTEx_expression_data(dir_path,first_complete):

    # 补充采集表达数据
    _, _, train_protein_official_symbols_aliases, _ = prepare_human_protein_symbols_ids()

    total_relative_protein_entries = set()
    for key in train_protein_official_symbols_aliases.keys():
        total_relative_protein_entries.add(key)
        for protein_entry in train_protein_official_symbols_aliases[key]:
            total_relative_protein_entries.add(protein_entry)
    print('total_relative_protein_entries: ' + str(len(total_relative_protein_entries)))

    completed_protein_entries = set()
    # total_supplementary_protein_entries = []
    completed_gene_expression_entries_fp = dir_path + 'gene_expression_in_all_tissues_by_tpm_median.csv'
    with codecs.open(completed_gene_expression_entries_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split(',')
            completed_protein_entries.add(temp[0].strip())

    print('completed_protein_entries: ' + str(len(completed_protein_entries)))

    supplementary_protein_entries = {}

    print('difference: ' + str(len(first_complete.difference(completed_protein_entries)))) # difference,union

    for entry in first_complete.difference(completed_protein_entries):
        supplementary_protein_entries[entry] = []
        if entry not in train_protein_official_symbols_aliases.keys():
            for entry_abbreviation in train_protein_official_symbols_aliases.keys():
                if entry in train_protein_official_symbols_aliases[entry_abbreviation]:
                    # if entry_abbreviation in first_complete:
                    #     print(entry)
                    for ea in train_protein_official_symbols_aliases[entry_abbreviation]:
                        if ea in completed_protein_entries:
                            supplementary_protein_entries[entry].append(ea)
        else:
            for e in train_protein_official_symbols_aliases[entry]:
                if e in completed_protein_entries:
                    supplementary_protein_entries[entry].append(e)

    print('supplementary_protein_entries: ' + str(len(supplementary_protein_entries.keys())))

    return supplementary_protein_entries

# Read the gene expression data on the HPA
def processing_hpa_expression_data(dir_path,supplementary_gene):

    gene_expression_fp = dir_path + 'the_human_protein_atlas/rna_tissue.tsv'
    gene_tissue_values = {}
    gene_names = set()
    with codecs.open(gene_expression_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split('	')
            if temp[1].strip() not in gene_names:
                gene_tissue_values[temp[1].strip()] = []
                gene_names.add(temp[1].strip())
            gene_tissue_values[temp[1].strip()].append([temp[2].strip(),temp[3].strip()])
    print('gene_tissue_values: '+str(len(gene_tissue_values.keys())))
    tissue_names_set = set()
    for gene in gene_tissue_values.keys():
        for tv in gene_tissue_values[gene]:
            tissue_names_set.add(tv[0].strip())
    print('tissue_names: ' + str(len(tissue_names_set)))

    tissue_names = list(tissue_names_set)
    gene_tissues_values = {}
    for gene in gene_tissue_values.keys():
        if gene not in supplementary_gene:
            continue
        tissue_values = []
        for i in range(len(tissue_names)):
            tissue_values.append(0.0)
        if gene in supplementary_gene:
            for tv in gene_tissue_values[gene]:
                tissue_values[tissue_names.index(tv[0])] = tv[1]
        gene_tissues_values[gene] = tissue_values

    return tissue_names,gene_tissues_values

# Supplement relevant expression data from HPA
def supplementary_expression_data_from_hpa(protein_entries_labeled):

    vp_exp_dir_path = '/vector_protein_exp/raw_data/'
    supplementary_gene = set()
    supplementary_protein_entries = supplementary_GTEx_expression_data(vp_exp_dir_path,protein_entries_labeled)
    for entry in supplementary_protein_entries.keys():
        if len(set(supplementary_protein_entries[entry]))== 0:
            supplementary_gene.add(entry)
    tissue_names, gene_tissues_values = processing_hpa_expression_data(vp_exp_dir_path,supplementary_gene)
    print(','.join(tissue_names))
    for gene in gene_tissues_values.keys():
        print(gene)
        print(gene_tissues_values[gene])
    dnd1_fp = vp_exp_dir_path + 'hpa_gtex_tissue_name_mapping.txt'
    gene_expression_in_all_tissues_fp = vp_exp_dir_path + 'gene_expression_in_all_tissues_by_tpm_median.csv'
    with open(gene_expression_in_all_tissues_fp, 'a') as output_file:
        for gene in gene_tissues_values.keys():
            hpa_values = gene_tissues_values[gene]
            gtex_values = []
            with codecs.open(dnd1_fp, "rb", "utf-8") as input_file:
                for line in islice(input_file.readlines(), 0, None):
                    temp = line.strip().split(',')
                    if temp[1].strip() == 'null':
                        gtex_values.append('0.0')
                    else:
                        gtex_values.append(hpa_values[tissue_names.index(temp[1].strip())])
            output_file.write(gene + ',' + ','.join(gtex_values))













