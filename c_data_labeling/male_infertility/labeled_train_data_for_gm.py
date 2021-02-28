#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 17/1/2019 11:34 AM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : labeled_train_data_for_gm.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import codecs
from itertools import islice
import pandas as pd


# Read in the disease name and phenotype name to match the protein name
def read_phenotype_terms_disorder_keywords(dir_path):

    # 读入表型类别名称及其对应表型短语
    classification_phenotype_terms = {}
    phenotype_gene_fp = dir_path + 'classification_label_for_graph_mining/Phenotypic classification of mice.xlsx'
    data_phenotype = pd.read_excel(phenotype_gene_fp,None)
    sheets = data_phenotype.keys()

    for sheet in sheets:
        sheet_data = data_phenotype[sheet]
        # columns_dict = {col:sheet_data[col].tolist() for col in sheet_data.columns}
        col = sheet_data.columns[0]
        terms = sheet_data[col].tolist()
        classification_phenotype_terms[sheet] = [col] + terms

    for key in classification_phenotype_terms.keys():
        print(key)
        print(len(classification_phenotype_terms[key]))
    print('classification_phenotype_terms: ' + str(len(classification_phenotype_terms.keys())))


    # 读入疾病类别名称及其对应疾病名称
    classification_disease_terms = {}
    disease_gene_fp = dir_path + 'classification_label_for_graph_mining/Disease label classification.xlsx'
    data_disease = pd.read_excel(disease_gene_fp, None)
    sheets = data_disease.keys()

    for sheet in sheets:
        sheet_data = data_disease[sheet]
        # columns_dict = {col: sheet_data[col].tolist() for col in sheet_data.columns}
        col = sheet_data.columns[0]
        terms = sheet_data[col].tolist()
        classification_disease_terms[sheet] = [col] + terms

    for key in classification_disease_terms.keys():
        print(key)
        print(len(classification_disease_terms[key]))
    print('classification_disease_terms: ' + str(len(classification_disease_terms.keys())))

    print('*************************************')

    return classification_phenotype_terms,classification_disease_terms


# Match the name of the gene to which each class belongs
def mapping_gene_set_for_classification(dir_path):

    # 读入训练集中打上阳性标签的基因名称
    genes_from_train_positive_samples = set()
    genes_from_train_positive_samples_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/total_train_dataset_with_positive_samples_by_disease_phenotype.csv'
    with codecs.open(genes_from_train_positive_samples_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            if temp == '':
                continue
            genes_from_train_positive_samples.add(temp)
    print('genes_from_train_positive_samples: ' + str(len(genes_from_train_positive_samples)))

    positive_samples_gene_all_aliases = {}
    positive_samples_gene_all_aliases_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/prepare_gene_symbols_abbreviation_for_train_model.csv'
    with codecs.open(positive_samples_gene_all_aliases_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            if line.strip() == '':
                continue
            temp = line.strip().split(',')
            if temp[0].strip() in genes_from_train_positive_samples:
                positive_samples_gene_all_aliases[temp[0].strip()] = temp[1:]
    print('positive_samples_gene_all_aliases: ' + str(len(positive_samples_gene_all_aliases.keys())))

    all_gene_aliases_from_train_positive_samples = set()
    for key in positive_samples_gene_all_aliases.keys():
        for gene in positive_samples_gene_all_aliases[key]:
            if gene.strip() == '':
                continue
            all_gene_aliases_from_train_positive_samples.add(gene.strip())
    print('all_gene_aliases_from_train_positive_samples: ' + str(len(all_gene_aliases_from_train_positive_samples)))

    # 读入表型或者疾病名称与基因的对应关系
    diorder_gene_fp = dir_path + 'label_related_dataset_final/several_gene_disease_db_summary.csv'
    phenotype_gene_fp = dir_path + 'label_related_dataset_final/mgi_gene_phenoes_db_summary.csv'

    gene_disorders = {}
    temp_genes_for_disorder = set()
    with codecs.open(diorder_gene_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) < 3:
                continue
            human_gene = temp[1].strip()
            if '_' in temp[1].strip():
                human_gene = temp[1].strip().split('_')[0].strip()
            elif '   ' in temp[1].strip():
                human_gene = temp[1].strip().split('   ')[0].strip()
            if human_gene not in all_gene_aliases_from_train_positive_samples:
                continue
            if temp[2].strip() == '':
                continue
            if human_gene not in temp_genes_for_disorder:
                gene_disorders[human_gene] = set()
                temp_genes_for_disorder.add(human_gene)
            gene_disorders[human_gene].add(temp[2].strip())
    print('gene_disorders: ' + str(len(gene_disorders.keys())))

    gene_phenotypes = {}
    temp_genes_for_phenotypes = set()
    with codecs.open(phenotype_gene_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) < 3:
                continue
            human_gene = temp[0].strip()
            if '_' in temp[0].strip():
                human_gene = temp[0].strip().split('_')[0].strip()
            elif '   ' in temp[0].strip():
                human_gene = temp[0].strip().split('   ')[0].strip()
            if human_gene not in all_gene_aliases_from_train_positive_samples:
                continue
            if temp[2].strip() == '':
                continue
            if human_gene not in temp_genes_for_phenotypes:
                gene_phenotypes[human_gene] = set()
                temp_genes_for_phenotypes.add(human_gene)
            phenotypes = set(temp[2].strip().split(';'))
            gene_phenotypes[human_gene] = gene_phenotypes[human_gene].union(phenotypes)
    print('gene_phenotypes: ' + str(len(gene_phenotypes.keys())))
    # print(len(set(gene_disorders.keys()).union(gene_phenotypes.keys())))

    print('-----------------------------------------')

    # 读入疾病与表型term信息

    classification_protein_symbols = {}
    classification_phenotype_terms, classification_disease_terms = read_phenotype_terms_disorder_keywords(dir_path)

    # classification_for_positive_samples = [u'精子发生',u'成熟精子及异常',u'受精及受精过程',u'胚胎',u'除睾丸外生殖疾病',u'睾丸相关生殖疾病',
    #                                        u'生育力',u'精子病理类型',u'生殖肿瘤',u'性别',u'其他疾病伴随生殖异常']

    gene_symbols_abbreviation = set()
    for classifi in classification_phenotype_terms.keys():
        phenotype_terms = classification_phenotype_terms[classifi]
        temp_protein_symbols = set()
        for phenotype_term in phenotype_terms:
            flag = 0
            for protein_symbol in gene_phenotypes:
                if phenotype_term in gene_phenotypes[protein_symbol]:
                    temp_protein_symbols.add(protein_symbol)
                    gene_symbols_abbreviation.add(protein_symbol)
                    flag = 1
            if flag == 0:
                print(classifi)
                print(phenotype_term)
        classification_protein_symbols[classifi] = temp_protein_symbols
    print('classification_protein_symbols: ' + str(len(classification_protein_symbols.keys())))

    print('++++++++++++++++++++++++++++++++++++++++++')

    for classifi in classification_disease_terms.keys():
        disease_terms = classification_disease_terms[classifi]
        temp_protein_symbols = set()
        for disease_term in disease_terms:
            flag = 0
            for protein_symbol in gene_disorders:
                if disease_term in gene_disorders[protein_symbol]:
                    temp_protein_symbols.add(protein_symbol)
                    gene_symbols_abbreviation.add(protein_symbol)
                    flag = 1
            if flag == 0:
                print(classifi)
                print(disease_term)
        if classifi in classification_protein_symbols.keys():
            classification_protein_symbols[classifi] = classification_protein_symbols[classifi].union(temp_protein_symbols)
        else:
            classification_protein_symbols[classifi] = temp_protein_symbols
    print('classification_protein_symbols: ' + str(len(classification_protein_symbols.keys())))
    print('gene_symbols_abbreviation: ' + str(len(gene_symbols_abbreviation)))

    gene_symbols_abbreviation_difference = gene_symbols_abbreviation.difference(genes_from_train_positive_samples)
    print(len(gene_symbols_abbreviation_difference))

    genes_from_train_positive_samples_difference = genes_from_train_positive_samples.difference(gene_symbols_abbreviation)
    print(len(genes_from_train_positive_samples_difference))

    belong_gene_mapping_offical_gene = {}
    for gene in genes_from_train_positive_samples_difference:
        for protein in positive_samples_gene_all_aliases[gene]:
            if protein in gene_symbols_abbreviation:
                belong_gene_mapping_offical_gene[protein] = gene
                break
    print(len(belong_gene_mapping_offical_gene.keys()))

    protein_symbols_selected = set()
    sample_classifications_protein_symbols = {}
    for key in classification_protein_symbols.keys():
        print(key)
        temp_proteins = set()
        protein_symbols_selected = protein_symbols_selected.union(classification_protein_symbols[key])
        for protein in classification_protein_symbols[key]:
            if protein in belong_gene_mapping_offical_gene.keys():
                temp_proteins.add(belong_gene_mapping_offical_gene[protein])
                protein_symbols_selected.add(belong_gene_mapping_offical_gene[protein])
            elif protein in genes_from_train_positive_samples:
                temp_proteins.add(protein)
                protein_symbols_selected.add(protein)
            else:
                continue
        sample_classifications_protein_symbols[key] = temp_proteins
        print(len(temp_proteins))
    print('protein_symbols_selected: ' + str(len(protein_symbols_selected)))

    print(len(genes_from_train_positive_samples.difference(protein_symbols_selected)))

    return sample_classifications_protein_symbols,gene_disorders,gene_phenotypes


# Calculate the multi-label vector for each protein
def compute_multi_labels_vector(dir_path):

    sample_classifications_protein_symbols,_,_ = mapping_gene_set_for_classification(dir_path)

    merge_classifi_for_positive_samples = {
        u'精子发生': [u'精子发生'],
        u'成熟精子及异常': [u'成熟精子及异常',u'精子病理类型'],
        u'受精过程及胚胎发育': [u'受精及受精过程', u'胚胎'],
        u'睾丸相关生殖疾病': [u'睾丸相关生殖疾病'],
        u'除睾丸外生殖疾病': [u'除睾丸外生殖疾病', u'性别'],
        u'生育力': [u'生育力'],
        u'生殖肿瘤': [u'生殖肿瘤'],
        u'其他疾病伴随生殖异常': [u'其他疾病伴随生殖异常']
    }

    count_sum = 0
    sample_classifications_protein_symbols_with_group = {}
    for key in merge_classifi_for_positive_samples.keys():
        temp_protein_symbols = set()
        for classifi in merge_classifi_for_positive_samples[key]:
            temp_protein_symbols = temp_protein_symbols.union(sample_classifications_protein_symbols[classifi])
        sample_classifications_protein_symbols_with_group[key] = temp_protein_symbols
        count_sum += len(temp_protein_symbols)
        print(key)
        print(len(temp_protein_symbols))
    print(count_sum)

    for out_key in sample_classifications_protein_symbols_with_group.keys():
        for inner_key in sample_classifications_protein_symbols_with_group.keys():
            if inner_key == out_key:
                continue
            ratio_coincide = float(len(sample_classifications_protein_symbols_with_group[out_key].intersection(sample_classifications_protein_symbols_with_group[inner_key]))) / len(sample_classifications_protein_symbols_with_group[out_key])
            if ratio_coincide < 0.1:
                continue
            print(out_key+' --- '+inner_key +': %f' % ratio_coincide)

    multi_labels_terms = [u'精子发生',u'成熟精子及异常',u'受精过程及胚胎发育',u'睾丸相关生殖疾病',
                                        u'除睾丸外生殖疾病',u'生育力',u'生殖肿瘤',u'其他疾病伴随生殖异常']
    # multi_labels_terms = list(merge_classifi_for_positive_samples.keys())
    print(multi_labels_terms)
    protein_multi_labels = {}
    temp_initi_protein_labels = set()
    for label_term in sample_classifications_protein_symbols_with_group.keys():
        for protein_symbol in sample_classifications_protein_symbols_with_group[label_term]:
            if protein_symbol not in temp_initi_protein_labels:
                protein_multi_labels[protein_symbol] = []
                for i in range(len(multi_labels_terms)):
                    protein_multi_labels[protein_symbol].append(0)
                temp_initi_protein_labels.add(protein_symbol)
            protein_multi_labels[protein_symbol][multi_labels_terms.index(label_term)] = 1


    out_pf = dir_path + 'classification_label_for_graph_mining/protein_multi_labels.csv'
    with open(out_pf, 'w') as output_file:
        for protein_symbol in protein_multi_labels.keys():
            print(protein_symbol, protein_multi_labels[protein_symbol])
            values = [str(x) for x in protein_multi_labels[protein_symbol]]
            output_file.write(protein_symbol.strip()+'$'+ ','.join(values) + '\n')


# Outputs a set of phenotypes and diseases corresponding to the specified gene/protein
def obtain_phenotype_and_disease_terms_by_protein_symbol(dir_path):

    sample_classifications_protein_symbols, gene_disorders, gene_phenotypes = mapping_gene_set_for_classification(dir_path)
    classification_phenotype_terms, classification_disease_terms = read_phenotype_terms_disorder_keywords(dir_path)
    intersection_two_class = sample_classifications_protein_symbols[u'生殖肿瘤'].intersection(sample_classifications_protein_symbols[u'直接生殖疾病'])

    phenotypes_terms = set()
    for classifi in classification_phenotype_terms.keys():
        phenotype_terms = classification_phenotype_terms[classifi]
        for phenotype_term in phenotype_terms:
            for protein_symbol in gene_phenotypes:
                if phenotype_term in gene_phenotypes[protein_symbol]:
                    if protein_symbol in intersection_two_class:
                        phenotypes_terms.add(phenotype_term)
    print('phenotypes_terms: ' + str(len(phenotypes_terms)))

    disease_terms = set()
    for classifi in classification_disease_terms.keys():
        disease_terms = classification_disease_terms[classifi]
        for disease_term in disease_terms:
            for protein_symbol in gene_disorders:
                if disease_term in gene_disorders[protein_symbol]:
                    if protein_symbol in intersection_two_class:
                        disease_terms.add(disease_term)
    print('disease_terms: ' + str(len(disease_terms)))

    # out_pf = dir_path + 'classification_label_for_graph_mining/phenotypes_gene_for_coincide_classifi.csv'
    out_pf = dir_path + 'classification_label_for_graph_mining/disease_gene_for_coincide_classifi.csv'
    with open(out_pf, 'w') as output_file:
        for term in disease_terms:
            output_file.write(term.strip() + '\n')


if __name__ == '__main__':

    dir_path = ''
    compute_multi_labels_vector(dir_path)
