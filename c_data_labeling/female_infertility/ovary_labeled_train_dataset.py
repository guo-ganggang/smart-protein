#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 26/1/2021 8:38 PM
# @Author  : ggguo
# @Site    : SYM PROJECT
# @File    : ovary_labeled_train_dataset.py
# @Software: SYM application case

###############
###Intermediate process code, for user reference only
###############

import codecs,re,os
from itertools import islice
import xlrd
import xlwt
from os import listdir
import Levenshtein

import sys

reload(sys)
sys.setdefaultencoding('utf-8')


################################################################
## 1. According to the gene expression in different tissues of GTEX, the differentially expressed gene combinations were selected
## 2. Read the high confidence non-redundant human egg protein collection
################################################################

def read_human_oocyte_genes(dir_path):

    oocyte_gene_fp = dir_path + 'High confidence non-redundant human egg protein collection.xls'
    try:
        data = xlrd.open_workbook(oocyte_gene_fp)
    except Exception, e:
        print str(e)
    table = data.sheet_by_name('Hoja1')
    uniprot_entry = table.col_values(2)[3:1379]
    # print(uniprot_entry[0],uniprot_entry[-1])
    # print('uniprot_entry: ' + str(len(uniprot_entry)))
    gene_id = table.col_values(3)[3:1379]
    # print(gene_id[0], gene_id[-1])
    print('high_confidence_nonredundancy_human_oocyte_genes: ' + str(len(gene_id)))

    out_file_path = dir_path + 'high_confidence_nonredundancy_human_oocyte_genes.csv'
    with open(out_file_path, 'w') as output_file:
        output_file.write('uniprot_entry'+'$'+'gene_id' + '\n')
        for i in range(len(uniprot_entry)):
            output_file.write(uniprot_entry[i] + '$' + gene_id[i] + '\n')
    return gene_id

# A list of highly expressed genes obtained based on RNASEQ gene expression data
def prepare_gene_list_base_hight_expression_in_gtex(dir_path):

    tissue_expression_from_gtex_fp = dir_path + 'ovary_hyperexpression/ovary.csv'
    human_gene_from_gtex = set()
    with codecs.open(tissue_expression_from_gtex_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split(',')
            human_gene_from_gtex.add(temp[0].strip())
    print('human_gene_from_gtex: ' + str(len(human_gene_from_gtex)))

    files_path = dir_path + 'ovary_hyperexpression/results_fc2/'
    file_names = [f for f in listdir(files_path) if f.startswith('ovary_')]
    total_gene_name_list_by_high_expression = human_gene_from_gtex
    for file_name in file_names:
        temp_gene_name = set()
        with codecs.open(files_path + file_name, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 1, None):
                temp = line.strip().split(',')
                gene_name = re.sub('"','',temp[0]).strip()
                temp_gene_name.add(gene_name)
        # print(len(temp_gene_name),len(total_gene_name_list_by_high_expression))
        total_gene_name_list_by_high_expression = total_gene_name_list_by_high_expression.intersection(temp_gene_name)
    print('total_gene_name_list_by_high_expression: ' + str(len(total_gene_name_list_by_high_expression)))

    out_file_path_total = dir_path + 'total_gene_names_with_overexpressed_in_ovary_by_deseq2.csv'
    with open(out_file_path_total, 'w') as output_file:
        for gene_name in total_gene_name_list_by_high_expression:
            output_file.write(gene_name + '\n')

    return total_gene_name_list_by_high_expression


################################################################
## 2. The potential disease names were matched
# from the three disease databases according to the disease keywords,
# and the positive gene sets were selected according to the disease names
##
################################################################

def get_disorder_keywords(dir_path):
    disease_gene_fp = dir_path + 'label_by_disorder/Key words female infertility.xlsx'
    try:
        data = xlrd.open_workbook(disease_gene_fp)
    except Exception, e:
        print str(e)
    table = data.sheet_by_name(u'keywords')
    disorder_keywords = table.col_values(0)
    print('disorder_keywords: ' + str(len(disorder_keywords)))

    return disorder_keywords

def data_write(file_path, datas):
    f = xlwt.Workbook()
    sheet1 = f.add_sheet('Disease Name', cell_overwrite_ok=True)  # 创建sheet

    # 将数据写入第 i 行，第 j 列
    i = 0
    for data in datas:
        sheet1.write(i, 0, data)
        i = i + 1
    f.save(file_path)  # 保存文件

def select_diorder_names_by_keywords_from_disease_db(dir_path):

    # 打标签的关键词
    disorder_keywords = get_disorder_keywords(dir_path)
    disorder_keywords_clean = set()
    for disorder_z in disorder_keywords:
        disorder_keywords_clean.add(disorder_z.lower())

    # 根据关键词匹配出三个数据库中的疾病名称
    disease_db_summary_fp = dir_path + 'label_by_disorder/several_gene_disease_db_summary.csv'
    disease_names_selected = set()
    with codecs.open(disease_db_summary_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if (temp[2].strip() == '') or (len(temp[2].strip())==1):
                continue
            disease_name = temp[2].strip()
            for keyword in disorder_keywords_clean:
                if keyword in disease_name:
                    disease_names_selected.add(disease_name)
                    break
    print('disease_names_selected: ' + str(len(disease_names_selected)))

    #  将数据写入文件,excel 文件与csv 文件
    disease_out_file_path = dir_path + 'label_by_disorder/matching_disorder_symbols/disease_names_by_disorder_labeled'
    data_write(disease_out_file_path+'.xls', disease_names_selected)
    with open(disease_out_file_path+'.csv', 'w') as output_file:
        for diease in disease_names_selected:
            output_file.write(diease.strip() + '\n')

def label_positive_samples_by_diorder_names(dir_path):

    disease_gene_fp = dir_path + 'label_by_disorder/disease_names_from_keywords.xls'
    try:
        data = xlrd.open_workbook(disease_gene_fp)
    except Exception, e:
        print str(e)
    table = data.sheet_by_name('Disease Name')
    target_disease_names = set(table.col_values(0))
    print('target_disease_names: '+str(len(target_disease_names)))

    disease_db_summary_fp = dir_path + 'label_by_disorder/several_gene_disease_db_summary.csv'
    genes_selected_from_disease_db = {}
    # extended_disease_names = set()
    with codecs.open(disease_db_summary_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if (temp[2].strip() == '') or (len(temp[2].strip()) == 1):
                continue
            disease_name = temp[2].strip()
            for target_disease_name in target_disease_names:
                if disease_name == target_disease_name:
                    if temp[1].strip() not in genes_selected_from_disease_db.keys():
                        genes_selected_from_disease_db[temp[1].strip()] = []
                    genes_selected_from_disease_db[temp[1].strip()].append(disease_name)
                    break
                # else:
                #     similarity_ratio = Levenshtein.ratio(disease_name, target_disease_name)
                #     if similarity_ratio > 0.9:
                #         disease_genes_selected.add(temp[1].strip())
                #         if disease_name not in target_disease_names:
                #             extended_disease_names.add(disease_name)

    print('genes_selected_from_disease_db: ' + str(len(genes_selected_from_disease_db)))
    # print('extended_disease_names: ' + str(len(extended_disease_names)))

    # #  将数据写入文件,excel 文件与csv 文件
    # extended_disease_names_fp = dir_path + 'label_by_disorder/extended_disease_names'
    # data_write(extended_disease_names_fp + '.xls', extended_disease_names)

    disease_genes_selected_fp = dir_path + 'label_by_disorder/genes_selected_from_disease_db.csv'
    with open(disease_genes_selected_fp, 'w') as output_file:
        for gene_name in genes_selected_from_disease_db:
            output_file.write(gene_name +'$'+ '$'.join(genes_selected_from_disease_db[gene_name]) + '\n')


################################################################
## 3. According to the phenotype name set,
# the gene set corresponding to human homologous mouse phenotype
# was selected as the positive gene set
##
################################################################

def get_phenotype_terms(dir_path):
    phenotype_terms_fp = dir_path + 'label_by_phenotype/Phenotype of female infertile mice.xlsx'
    try:
        data = xlrd.open_workbook(phenotype_terms_fp)
    except Exception, e:
        print str(e)
    table = data.sheet_by_name('names of the phenotypic')
    phenotype_terms = set(table.col_values(0))
    print('phenotype_terms: ' + str(len(phenotype_terms)))

    return phenotype_terms

def labeled_positive_example_by_phenotype(dir_path):

    target_phenotype_terms = get_phenotype_terms(dir_path)
    protein_entries_labeled_by_phenotype = {}
    gene_phenoes_fp = dir_path + 'label_by_phenotype/mgi_gene_phenoes_db_summary.csv'
    with codecs.open(gene_phenoes_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) < 3:
                continue
            line_phenotype_terms = temp[2].strip().split(';')
            human_protein_name = temp[0].strip().split('_')[0].strip()
            for phenotype_term in target_phenotype_terms:
                if phenotype_term in line_phenotype_terms:
                    if human_protein_name not in protein_entries_labeled_by_phenotype.keys():
                        protein_entries_labeled_by_phenotype[human_protein_name] = []
                    protein_entries_labeled_by_phenotype[human_protein_name].append(temp[2].strip())
                    break
    print('protein_entries_labeled_by_phenotype: '+str(len(protein_entries_labeled_by_phenotype)))

    out_file_path = dir_path + 'label_by_phenotype/protein_entries_by_phenotype_labeled.csv'
    with open(out_file_path, 'w') as output_file:
        for protein_symbol in protein_entries_labeled_by_phenotype.keys():
            output_file.write(protein_symbol.strip() + '$' + '$'.join(protein_entries_labeled_by_phenotype[protein_symbol]) + '\n')


################################################################
## 4. Negative labels were obtained from mouse phenotypes according to positive labels
##
################################################################

# Add a few manual screening positive cases
def read_some_positive_label_genes(dir_path):
    phenotype_terms_fp = dir_path + 'manual_selecting_positive_label_genes.xls'
    try:
        data = xlrd.open_workbook(phenotype_terms_fp)
    except Exception, e:
        print str(e)
    table = data.sheet_by_name(u'阳性标签的基因')
    some_positive_label_genes = set(table.col_values(0))
    print('some_positive_label_genes: ' + str(len(some_positive_label_genes)))

    return some_positive_label_genes

def prepare_human_protein_official_symbol_names(dir_path):

    official_protein_symbol_fp = dir_path + 'human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    human_protein_official_symbol_names = {}
    with codecs.open(official_protein_symbol_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            human_protein_official_symbol_names[temp[0]] = set(temp[1:])
    print('human_protein_official_symbol_names: ' + str(len(human_protein_official_symbol_names.keys())))

    return human_protein_official_symbol_names

# Combined all positive samples, the name of the unified treatment
def prepare_human_official_protein_names(dir_path):

    some_positive_label_genes = read_some_positive_label_genes(dir_path)

    human_protein_official_symbol_names = prepare_human_protein_official_symbol_names(dir_path)

    genes_labeled_positive_samples = set()
    gene_phenoes_fp = dir_path + 'label_by_phenotype/protein_entries_by_phenotype_labeled.csv'
    with codecs.open(gene_phenoes_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            genes_labeled_positive_samples.add(temp[0].strip())
    gene_disease_fp = dir_path + 'label_by_disorder/genes_selected_from_disease_db.csv'
    with codecs.open(gene_disease_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            genes_labeled_positive_samples.add(temp[0].strip())

    genes_labeled_positive_samples = genes_labeled_positive_samples.union(some_positive_label_genes)
    print('genes_labeled_positive_samples: ' + str(len(genes_labeled_positive_samples)))

    gene_names_for_positive_samples_standardization = set()
    for gene_name in genes_labeled_positive_samples:
        if gene_name in human_protein_official_symbol_names.keys():
            gene_names_for_positive_samples_standardization.add(gene_name)
            continue
        for key in human_protein_official_symbol_names.keys():
            if gene_name in human_protein_official_symbol_names[key]:
                gene_names_for_positive_samples_standardization.add(gene_name)
                break
    print('gene_names_for_positive_samples_standardization: ' + str(len(gene_names_for_positive_samples_standardization)))

    out_file_path = dir_path + 'gene_names_for_positive_samples_standardization.csv'
    with open(out_file_path, 'w') as output_file:
        for protein_symbol in gene_names_for_positive_samples_standardization:
            output_file.write(protein_symbol.strip() + '\n')

    return gene_names_for_positive_samples_standardization

# Negative samples were screened
def labeled_train_data_negative_example_by_mammalian_phenotype(dir_path):

    genes_labeled_positive_samples_standardization = prepare_human_official_protein_names(dir_path)

    # 匹配是否条件敲除
    file_path_gene_knockout = dir_path + 'label_by_phenotype/mgi_mouse_gene_knockout_protein_info.csv'

    gene_names_from_cko = set()
    total_samples = set()
    with codecs.open(file_path_gene_knockout, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            total_samples.add(temp[0].strip())
            results = temp[4].strip().split('!!')
            for re in results:
                if re[:2] == 'cn':
                    if temp[0].strip() not in genes_labeled_positive_samples_standardization:
                        gene_names_from_cko.add(temp[0].strip())
                        break
    print('total_samples: ' + str(len(total_samples)))
    print('gene_names_from_cko: '+str(len(gene_names_from_cko)))

    gene_names_for_negative_example = (total_samples.difference(genes_labeled_positive_samples_standardization)).difference(gene_names_from_cko)

    human_protein_official_symbol_names = prepare_human_protein_official_symbol_names(dir_path)

    gene_names_for_negative_example_standardization = set()
    for gene_name in gene_names_for_negative_example:
        if gene_name in human_protein_official_symbol_names.keys():
            gene_names_for_negative_example_standardization.add(gene_name)
            continue
        for key in human_protein_official_symbol_names.keys():
            if gene_name in human_protein_official_symbol_names[key]:
                gene_names_for_negative_example_standardization.add(gene_name)
                break

    print('gene_names_for_negative_example_standardization: '+str(len(gene_names_for_negative_example_standardization)))

    out_file_path = dir_path + 'gene_names_for_negative_example_standardization.csv'
    with open(out_file_path, 'w') as output_file:
        for protein_symbol in gene_names_for_negative_example_standardization:
            output_file.write(protein_symbol.strip() + '\n')

    return gene_names_for_negative_example_standardization


################################################################
## 5. All gene data were merged to generate the disease gene domain,
# and each gene in the disease domain was represented by matching vector
##
################################################################

# All gene data were combined to generate disease gene domains
def merge_selected_gene_symbols_for_seed_universe(dir_path):

    high_confidence_nonredundancy_human_oocyte_genes = set(read_human_oocyte_genes(dir_path))
    total_gene_name_list_by_high_expression = prepare_gene_list_base_hight_expression_in_gtex(dir_path)
    genes_labeled_positive_samples_standardization = prepare_human_official_protein_names(dir_path)
    gene_names_for_negative_example_standardization = labeled_train_data_negative_example_by_mammalian_phenotype(dir_path)

    seed_universe_gene_symbols = high_confidence_nonredundancy_human_oocyte_genes.union(total_gene_name_list_by_high_expression.union(genes_labeled_positive_samples_standardization.union(gene_names_for_negative_example_standardization)))

    out_file_path = dir_path + 'seed_universe_gene_symbols_female_infertility.csv'
    with open(out_file_path, 'w') as output_file:
        for protein_symbol in seed_universe_gene_symbols:
            output_file.write(protein_symbol.strip() + '\n')

    print('seed_universe_gene_symbols: '+str(len(seed_universe_gene_symbols)))



if __name__ == '__main__':

    dir_path = ''
    # read_human_oocyte_genes(dir_path)
    # prepare_gene_list_base_hight_expression_in_gtex(dir_path)
    # select_diorder_names_by_keywords_from_disease_db(dir_path)
    # label_positive_samples_by_diorder_names(dir_path)
    # labeled_positive_example_by_phenotype(dir_path)
    # labeled_train_data_negative_example_by_mammalian_phenotype(dir_path)
    merge_selected_gene_symbols_for_seed_universe(dir_path)
