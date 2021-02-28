#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 20/10/2018 12:01 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : data_process_for_semsim.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import codecs,re,os
from itertools import islice
import gzip
import xlrd
import pronto
from collections import OrderedDict
from os import listdir
import Levenshtein

import sys

reload(sys)
sys.setdefaultencoding('utf-8')


################################################################
## 1. Training set labeling based on human sperm specific protein sets
##
################################################################

def prepare_human_protein_names():
    raw_data_fp = 'C:/ggguo/1_data_help/protein_entities_standardized/'
    out_file_path = raw_data_fp + 'human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    human_protein_official_symbol_names = {}
    with codecs.open(out_file_path, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            human_protein_official_symbol_names[temp[0]] = set(temp[1:])
    print('human_protein_official_symbol_names: '+str(len(human_protein_official_symbol_names.keys())))
    return human_protein_official_symbol_names

def prepare_human_sperm_protein_name():
    gene_symbols_fp = 'C:/ggguo/1_data_help/protein_entities_standardized/'
    human_sperm_protein_aliases_names_fp = gene_symbols_fp + 'human_sperm_protein_name_all_abbreviation_aliases.csv'
    all_abbreviation_names = {}
    with codecs.open(human_sperm_protein_aliases_names_fp, 'rb', 'utf-8') as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            all_abbreviation_names[temp[0]] = temp[1:]
    print('human_sperm_protein_name_all_abbreviation_aliases: '+str(len(all_abbreviation_names.keys())))

    return all_abbreviation_names

def labeled_train_dataset_by_phenotype(dir_path,phenotype_terms):

    protein_entries_labeled_by_phenotype = set()
    total_protein_entries_labeled_by_phenotype = set()
    knockout_gene_fp = dir_path + 'raw_data/mgi_mouse_gene_knockout_protein_info.csv'
    with codecs.open(knockout_gene_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if temp[3].strip() in phenotype_terms:
                protein_entries_labeled_by_phenotype.add(temp[0].strip())
            total_protein_entries_labeled_by_phenotype.add(temp[0].strip())
    print('protein_entries_labeled_by_phenotype: '+str(len(protein_entries_labeled_by_phenotype)))
    print('total_protein_entries_labeled_by_phenotype: ' + str(len(total_protein_entries_labeled_by_phenotype)))
    return protein_entries_labeled_by_phenotype,total_protein_entries_labeled_by_phenotype

def phenotype_terms_disorder_keywords(dir_path):
    disease_gene_fp = dir_path + 'raw_data/label_for_phenotype_disorder.xlsx'
    try:
        data = xlrd.open_workbook(disease_gene_fp)
    except Exception, e:
        print str(e)
    phenotype_terms = []
    disorder_keywords = []
    sheets = [u'受精及受精过程', u'生育力', u'成熟精子及异常', u'精子病理类型', u'精子发生', u'表观遗传学', u'男性生殖系统', u'胚胎', u'关键词']
    for sheet in sheets:
        if sheet != u'关键词':
            table = data.sheet_by_name(sheet)
            phenotype_terms += table.col_values(0)
        else:
            table = data.sheet_by_name(sheet)
            disorder_keywords = table.col_values(0)
    print('phenotype_terms: ' + str(len(phenotype_terms)))
    print('disorder_keywords: ' + str(len(disorder_keywords)))

    return phenotype_terms,disorder_keywords

def generate_mammalian_phenotype_label_tree():

    dir_path = 'C:/ggguo/1_data_help/mgi/'
    count = 25

    # 读入基因名称
    file_path = dir_path + 'raw_data/mgi_mouse_gene_knockout_protein_info.csv'
    various_labels = {}
    with codecs.open(file_path, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if temp[2].strip() not in various_labels.keys():
                various_labels[temp[2].strip()] = []
            various_labels[temp[2].strip()].append(temp[3].strip())
    print(len(various_labels))

    # 生成数据格式
    # {
    #     "name": "Openi 数据库(图片)",
    #     "children": [
    #         {
    #             "name": "第一阶段 采集 Openi 图片 ID", "value": 21,
    #             "children": [
    #                 {"name": "爬虫主体采集", "value": 211},
    #                 {"name": "API 补充采集", "value": 212}
    #             ]
    #         },
    #         {
    #             "name": "第二阶段 采集 图片自动摘要信息 ", "value": 21,
    #             "children": [
    #                 {"name": "PMID 与 PMCID 的转换", "value": 211},
    #                 {"name": "采集文献图片相关信息", "value": 212}
    #             ]
    #         }
    #     ]
    # }

    tree_dict = []
    for key in various_labels.keys():
        temp_dict = {}
        for label in various_labels[key]:
            if label not in temp_dict.keys():
                temp_dict[label] = 1
            else:
                temp_dict[label] += 1
        print(key,str(len(temp_dict)))

        second_class_list = []
        if len(temp_dict) <= count:
            temp_third_class = {}
            temp_third_class['name'] = 'Whole label'
            temp_third_class['children'] = []
            for last_label in temp_dict.keys():
                temp_third_class['children'].append({"name": str(last_label), "value": temp_dict[last_label]})
            second_class_list.append(temp_third_class)
        else:
            for i in range(len(temp_dict)//count):
                temp_third_class = {}
                temp_third_class['name'] = 'Part '+str(i)
                temp_third_class['children'] = []
                for last_label in temp_dict.keys()[(i*count):((i+1)*count)]:
                    temp_third_class['children'].append({"name": str(last_label), "value": temp_dict[last_label]})
                second_class_list.append(temp_third_class)
            temp_third_class = {}
            temp_third_class['name'] = 'Part ' + str(len(temp_dict)//count)
            temp_third_class['children'] = []
            for last_label in temp_dict.keys()[-(len(temp_dict)%count):]:
                temp_third_class['children'].append({"name": str(last_label), "value": temp_dict[last_label]})
            second_class_list.append(temp_third_class)

        tree_dict.append({"name": str(key),"children": second_class_list})
        # print(str(key), str(len(second_class_list)))

    print('---------------------------')
    standard_label_entries = set()
    with open('D:/3_label_properties/standard_label_entries.csv', 'w') as output_file:
        for tissue in tree_dict:
            for part in tissue['children']:
                for name_value in part['children']:
                    # print(name_value['name'])
                    if name_value['name'].strip() not in standard_label_entries:
                        standard_label_entries.add(name_value['name'].strip())
                        output_file.write(name_value['name'].strip()+'\n')

    filter_label_entries = set()
    with codecs.open('D:/3_label_properties/labels.txt', "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            if temp not in standard_label_entries:
                print(temp)
            filter_label_entries.add(temp)
    return filter_label_entries

def generate_human_disorder_label_tree(dir_path):

    # 读入GENECARDS 疾病与基因对应关系，疾病有很多别称
    file_path = dir_path + 'label_by_disorders/protein_disorder_from_genecards.csv'
    various_labels = {}
    protein_symbol_disorder = []
    with codecs.open(file_path, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if temp[1].strip() not in various_labels.keys():
                if temp[2].strip() == '':
                    various_labels[temp[1].strip().lower()] = [temp[1].strip().lower()]
                else:
                    various_labels[temp[1].strip().lower()] = [temp[1].strip().lower()]+[x.lower() for x in temp[2].strip().split(';')]
            else:
                various_labels[temp[1].strip().lower()] += [x.lower() for x in temp[2].strip().split(';')]
            if [temp[0].strip(),temp[1].strip().lower()] not in protein_symbol_disorder:
                protein_symbol_disorder.append([temp[0].strip(),temp[1].strip().lower()])
    print(len(various_labels.keys()),len(protein_symbol_disorder))

    # 读入humsavar 与clinvar 上基因与疾病的对应关系对
    other_protein_symbol_disorder = []

    file_path_humsavar = dir_path + 'label_by_disorders/humsavar_gene_disease.txt'
    with codecs.open(file_path_humsavar, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 49, None):
            if 'Disease' not in line.strip():
                continue
            temp = line.strip().split('    ')
            disease = temp[-1].strip().split(' [MIM:')[0].strip().lower()
            if disease.startswith('rs'):
                disease = disease[12:].strip()
            if [temp[0].strip(),disease] not in other_protein_symbol_disorder:
                other_protein_symbol_disorder.append([temp[0].strip(),disease])
    print('HUMSAVAR: ' + str(len(other_protein_symbol_disorder)))

    file_path_clinvar = dir_path + 'label_by_disorders/clinvar_gene_condition_source_id.txt'
    with codecs.open(file_path_clinvar, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split('	')
            if len(temp) != 9:
                continue
            if [temp[1].strip(), temp[3].strip().lower()] not in other_protein_symbol_disorder:
                other_protein_symbol_disorder.append([temp[1].strip(), temp[3].strip().lower()])
        print('clinvar + HUMSAVAR: ' + str(len(other_protein_symbol_disorder)))

    out_file_path = dir_path + 'label_by_disorders/statistic_disorder_aliases.csv'
    if not os.path.isfile(out_file_path):
        with open(out_file_path, 'a') as output_file:
            for key in various_labels.keys():
                aliases = list(set(various_labels[key]))
                output_file.write(key + '$' + ';'.join(aliases) + '\n')

    return protein_symbol_disorder, various_labels,other_protein_symbol_disorder

# mgi
def labeled_train_data_positive_example_by_mammalian_phenotype(dir_path):

    phenotype_terms, _ = phenotype_terms_disorder_keywords(dir_path)
    protein_entries_labeled_by_phenotype,_ = labeled_train_dataset_by_phenotype(dir_path,phenotype_terms)
    all_abbreviation_names = prepare_human_sperm_protein_name()

    out_file_path = dir_path + 'label_by_phenotypes/protein_entries_by_phenotype_labeled.csv'
    is_protein_names = set()
    if not os.path.exists(out_file_path):
        with open(out_file_path, 'w') as output_file:
            for phenotype in protein_entries_labeled_by_phenotype:
                flag = 0
                if (phenotype.strip() in all_abbreviation_names.keys()) and (phenotype.strip() not in is_protein_names):
                    output_file.write(phenotype.strip() + '\t' + '1' + '\n')
                    is_protein_names.add(phenotype.strip())
                    continue
                else:
                    for key in all_abbreviation_names.keys():
                        other_names = all_abbreviation_names[key]
                        if (phenotype.strip() in other_names) and (key.strip() not in is_protein_names):
                            output_file.write(key.strip() + '\t' + '1' + '\n')
                            is_protein_names.add(key.strip())
                            flag = 1
                            break
                if flag == 0:
                    output_file.write(phenotype.strip() + '\t' + '0' + '\n')
    print('protein_names: ' + str(len(is_protein_names)))

# mgi,genecards(...),humsavar,clinvar
def labeled_train_data_positive_example_by_human_disease(dir_path):

    _, disorder_keywords = phenotype_terms_disorder_keywords(dir_path)

    disorder_keywords_clean = set()
    for disorder_z in disorder_keywords:
        disorder_keywords_clean.add(disorder_z.lower())

    disease_gene_fp = dir_path + 'label_by_disorders/mgi_gene_disease.xlsx'
    try:
        data = xlrd.open_workbook(disease_gene_fp)
    except Exception, e:
        print str(e)

    mgi_gene_disease = []
    table = data.sheet_by_name('Sheet0')
    # colnames = table.row_values(0)
    # print(colnames)
    diseases = table.col_values(1)
    genes = table.col_values(3)

    print(len(diseases),len(genes))
    # print(diseases[0],genes[0])
    for i in range(1,len(diseases)):
        if genes[i].strip() == '':
            continue
        for gene_s in genes[i].strip().split('*'):
            if gene_s.strip() == '':
                continue
            mgi_gene_disease.append([gene_s.strip(),diseases[i].strip().lower()])
    print('mgi_gene_disease: ' + str(len(mgi_gene_disease)))

    filter_wrong_disorders_fp = dir_path + 'label_by_disorders/filter_wrong_disorders.csv'
    filter_wrong_disorders = set()
    with codecs.open(filter_wrong_disorders_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().lower()
            filter_wrong_disorders.add(temp)
    print('filter_wrong_disorders: '+str(len(filter_wrong_disorders)))

    protein_symbol_disorder, various_labels, other_protein_symbol_disorder = generate_human_disorder_label_tree(
        dir_path)
    all_abbreviation_names = prepare_human_sperm_protein_name()

    protein_entries_labeled_by_disorder = set()
    disorders_by_labeled = set()
    break_flag = 0
    protein_out_file_path = dir_path + 'label_by_disorders/protein_entries_by_disorder_labeled.csv'
    disease_out_file_path = dir_path + 'label_by_disorders/disease_names_by_disorder_labeled.csv'
    if (not os.path.exists(disease_out_file_path)) or (not os.path.exists(protein_out_file_path)):
        for gd in mgi_gene_disease:
            for keyword in disorder_keywords_clean:
                if gd[1].strip() in filter_wrong_disorders:
                    continue
                if keyword == gd[1].strip():
                    protein_entries_labeled_by_disorder.add(gd[0])
                    disorders_by_labeled.add(gd[1])
                    break
                if keyword in gd[1].strip().split(' '):
                    if keyword == 'sperm':
                        for disease_word in gd[1].strip().split(' '):
                            if keyword in disease_word:
                                protein_entries_labeled_by_disorder.add(gd[0])
                                disorders_by_labeled.add(gd[1])
                                break_flag = 1
                                break
                        if break_flag == 1:
                            break
                    else:
                        protein_entries_labeled_by_disorder.add(gd[0])
                        disorders_by_labeled.add(gd[1])
                        break
        print('mgi_protein_entries_labeled_by_disorder: ' + str(len(protein_entries_labeled_by_disorder)))

        for psd in protein_symbol_disorder:
            flag = 0
            for keyword in disorder_keywords_clean:
                diseases = various_labels[psd[1]]
                if psd[1] in filter_wrong_disorders:
                    continue
                if keyword in diseases:
                    protein_entries_labeled_by_disorder.add(psd[0])
                    disorders_by_labeled.add(psd[1])
                    break
                for disease in diseases:
                    disease = re.sub(',',' ',disease)
                    disease = re.sub('/', ' ', disease)
                    disease_words = disease.strip().split(' ')
                    disease_words_temp = []
                    for dw in disease_words:
                        disease_words_temp.append(dw.strip())
                    if keyword in disease_words_temp:
                        protein_entries_labeled_by_disorder.add(psd[0])
                        disorders_by_labeled.add(psd[1])
                        flag = 1
                        break
                    if keyword == 'sperm':
                        for disease_word in disease_words_temp:
                            if keyword in disease_word:
                                protein_entries_labeled_by_disorder.add(psd[0])
                                disorders_by_labeled.add(psd[1])
                                flag = 1
                                break
                if flag == 1:
                    break
        print('genecards_protein_entries_labeled_by_disorder: ' + str(len(protein_entries_labeled_by_disorder)))

        for opsd in other_protein_symbol_disorder:
            flag = 0
            for keyword in disorder_keywords_clean:
                if keyword == opsd[1].strip():
                    protein_entries_labeled_by_disorder.add(opsd[0])
                    disorders_by_labeled.add(opsd[1])
                    break
                if keyword in opsd[1].strip().split(' '):
                    protein_entries_labeled_by_disorder.add(opsd[0])
                    disorders_by_labeled.add(opsd[1])
                    break
                if keyword == 'sperm':
                    for disease_word in opsd[1].strip().split(' '):
                        if keyword in disease_word:
                            protein_entries_labeled_by_disorder.add(opsd[0])
                            disorders_by_labeled.add(opsd[1])
                            flag = 1
                            break
                if flag == 1:
                    break
        print('humsavar_clinvar_protein_entries_labeled_by_disorder: ' + str(len(protein_entries_labeled_by_disorder)))

        with open(disease_out_file_path, 'w') as output_file:
            for diease in disorders_by_labeled:
                output_file.write(diease.strip() + '\n')

        is_protein_names = set()
        with open(protein_out_file_path, 'w') as output_file:
            for protein_symbol in protein_entries_labeled_by_disorder:
                # print(protein_symbol)
                protein_symbol = protein_symbol.split('   ')[0].strip()
                flag = 0
                if (protein_symbol.strip() in all_abbreviation_names.keys()) and (protein_symbol.strip() not in is_protein_names):
                    output_file.write(protein_symbol.strip() + '\t' + '1' + '\n')
                    is_protein_names.add(protein_symbol.strip())
                    continue
                else:
                    for key in all_abbreviation_names.keys():
                        other_names = all_abbreviation_names[key]
                        if (protein_symbol.strip() in other_names) and (key.strip() not in is_protein_names):
                            is_protein_names.add(key.strip())
                            output_file.write(key.strip() + '\t' + '1' + '\n')
                            flag = 1
                            break
                if flag == 0:
                    output_file.write(protein_symbol.strip() + '\t' + '0' + '\n')

        print(len(is_protein_names))

    # 所有基因与疾病对应关系
    file_path_testis_overexpression = dir_path + 'training_set_extension_by_testis/label_by_disorder/collect_disorder_differential_expression_from_genecards_deseq2_v2.csv'
    with codecs.open(file_path_testis_overexpression, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if [temp[0].strip(), temp[1].strip().lower()] not in protein_symbol_disorder:
                protein_symbol_disorder.append([temp[0].strip(), temp[1].strip().lower()])

    several_gene_disease_db_gene_disorders = {'db_mgi':mgi_gene_disease,'db_humsavar_clinvar':other_protein_symbol_disorder,'db_genecards':protein_symbol_disorder}
    file_path_disease_gene = dir_path + 'label_related_dataset_final/several_gene_disease_db_summary.csv'

    if not os.path.exists(file_path_disease_gene):
        with open(file_path_disease_gene, 'w') as output_file:
            for key in several_gene_disease_db_gene_disorders.keys():
                db_gene_diseases = several_gene_disease_db_gene_disorders[key]
                for gene_diseas in db_gene_diseases:
                    output_file.write(key + '$' + gene_diseas[0].strip() + '$' + gene_diseas[1].strip() + '\n')

# Summary of Positive Examples
def merge_labeled_train_data_positive_example(dir_path):

    file_path_phenotypes = dir_path + 'label_by_phenotypes/protein_entries_by_phenotype_labeled.csv'
    file_path_disorders = dir_path + 'label_by_disorders/protein_entries_by_disorder_labeled.csv'
    gene_names_for_positive_example = set()
    file_pathes = [file_path_phenotypes, file_path_disorders]
    for file_path in file_pathes:
        with codecs.open(file_path, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('\t')
                gene_names_for_positive_example.add(temp[0].strip())
    # print('mgi_gene_names_for_positive_example: '+str(len(gene_names_for_positive_example)))

    return gene_names_for_positive_example

# Negative sample screening
def labeled_train_data_negative_example_by_mammalian_phenotype(dir_path):

    protein_entries_labeled = merge_labeled_train_data_positive_example(dir_path)

    # 匹配是否条件敲除
    file_path_gene_knockout = dir_path + 'raw_data/mgi_mouse_gene_knockout_protein_info.csv'

    gene_names_from_cko = set()
    total_samples = set()
    with codecs.open(file_path_gene_knockout, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            total_samples.add(temp[0].strip())
            results = temp[4].strip().split('!!')
            for re in results:
                if re[:2] == 'cn':
                    if temp[0].strip() not in protein_entries_labeled:
                        gene_names_from_cko.add(temp[0].strip())
                        break
    # print(len(gene_names_from_cko))

    gene_names_for_negative_example = (total_samples.difference(protein_entries_labeled)).difference(gene_names_from_cko)

    print('mgi_gene_names_for_negative_example_hs: '+str(len(gene_names_for_negative_example)))

    return gene_names_for_negative_example

# Due to the change of training set generation strategy,
# only features were extracted from human data,
# and the positive label of mouse phenotype was combined
# with the positive label of human disease
def merge_phenotype_disorder_for_train_data_label(dir_path):

    # 阳性例子
    gene_names_for_positive_example = merge_labeled_train_data_positive_example(dir_path)

    # 阴性例子
    gene_names_for_negative_example = labeled_train_data_negative_example_by_mammalian_phenotype(dir_path)

    gene_names_from_train_data_samples = set(list(gene_names_for_positive_example)+list(gene_names_for_negative_example))

    print(len(gene_names_from_train_data_samples))

    out_file_path = dir_path + 'gene_symbols_from_train_data_samples.csv'
    if not os.path.exists(out_file_path):
        with open(out_file_path, 'w') as output_file:
            for gene_name in gene_names_from_train_data_samples:
                if gene_name in gene_names_for_positive_example:
                    output_file.write(gene_name + ',' + '1' + '\n')
                else:
                    output_file.write(gene_name + ',' + '0' + '\n')
    return gene_names_from_train_data_samples

# Ensure that all samples in the training data set are in the full sample database
def temp_match_train_dataset_total_dataset(dir_path):

    train_dataset_fp = dir_path + 'gene_symbols_from_train_data_samples.csv'
    total_dataset_fp = dir_path + 'human_sperm_protein_name_all_abbreviation_aliases.csv'
    train_dataset_gene_names = set()
    total_dataset_gene_names = {}
    with codecs.open(train_dataset_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            train_dataset_gene_names.add(temp[0].strip())

    with codecs.open(total_dataset_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            total_dataset_gene_names[temp[0].strip()] = set(temp)

    suppliment_gene_names = set()
    for gene_name in train_dataset_gene_names:
        flag = 0
        if gene_name in total_dataset_gene_names.keys():
            continue
        else:
            for key in total_dataset_gene_names.keys():
                if gene_name in total_dataset_gene_names[key]:
                    flag = 1
                    break
        if flag == 0:
            print(gene_name)
            suppliment_gene_names.add(gene_name)
    print(len(suppliment_gene_names))

################################################################
## 2. Training set labeling based on testes highly expressed protein sets
##
################################################################

# MGI phenotypic data that meets the LABEL condition,
# and those that are already included in the training dataset are removed
def mapping_total_potential_label_dataset_from_mgi(dir_path):

    mpheno_terms_fp = dir_path + 'training_set_extension_by_testis/label_by_phenotype/MPheno_OBO.ontology'
    mp_name_id = {}
    ont = pronto.Ontology(mpheno_terms_fp)
    for term in ont:
        mp_name_id[term.name.strip()] = term.id.strip()
    print(len(mp_name_id.keys()))

    phenotype_terms, _ = phenotype_terms_disorder_keywords(dir_path)

    phenotype_terms_ids_for_label = set()
    for phenotype_term in phenotype_terms:
        if phenotype_term not in mp_name_id.keys():
            print('miss this phenotype term: '+phenotype_term)
            continue
        phenotype_terms_ids_for_label.add(mp_name_id[phenotype_term])
    print(len(phenotype_terms_ids_for_label))

    mgi_pheno_geno_mp_fp = dir_path + 'training_set_extension_by_testis/label_by_phenotype/MGI_PhenoGenoMP.rpt'
    mgi_gene_phenoes = {}
    mgi_mouse_gene_symbols_have_phenoes = set()
    allele_symbols_ste_set = set()
    with codecs.open(mgi_pheno_geno_mp_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            allele_symbols_str = temp[1].strip()
            temp_gene_names = set()
            if allele_symbols_str not in allele_symbols_ste_set:
                allele_symbols = allele_symbols_str.split('|')
                for allele_symbol in allele_symbols:
                    if '<' not in allele_symbol:
                        continue
                    gene_name = allele_symbol.strip().split('<')[0].strip()
                    if '/' in gene_name:
                        ges = gene_name.split('/')
                        for ge in ges:
                            mgi_gene_phenoes[ge.strip()] = set()
                            temp_gene_names.add(ge.strip())
                            mgi_mouse_gene_symbols_have_phenoes.add(ge.strip())
                    else:
                        mgi_gene_phenoes[gene_name] = set()
                        temp_gene_names.add(gene_name)
                        mgi_mouse_gene_symbols_have_phenoes.add(gene_name)
                allele_symbols_ste_set.add(allele_symbols_str)
                for gene_name in temp_gene_names:
                    mgi_gene_phenoes[gene_name].add(temp[3].strip())
    print(len(mgi_gene_phenoes.keys()))

    mouse_gene_names_for_label = set()
    for gene_name in mgi_gene_phenoes.keys():
        if len(mgi_gene_phenoes[gene_name]) == 0:
            continue
        for pheno in mgi_gene_phenoes[gene_name]:
            if pheno in phenotype_terms_ids_for_label:
                mouse_gene_names_for_label.add(gene_name)
                break
    print('gene_names_for_label: '+ str(len(mouse_gene_names_for_label)))

    mgi_homolo_gene_fp = dir_path + 'training_set_extension_by_testis/label_by_phenotype/HMD_HumanPhenotype.rpt'
    mouse_human_homolo_gene = {}
    with codecs.open(mgi_homolo_gene_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) < 5:
                continue
            mouse_human_homolo_gene[temp[4].strip()] = temp[0].strip()
    print('The number of the mgi mouse gene: '+str(len(mouse_human_homolo_gene.keys())))

    # 有表型数据且和人类同源基因
    mgi_human_homologous_gene_symbols_have_phenoes = set()
    mgi_mouse_gene_symbols_have_phenoes_fp = dir_path + 'training_set_extension_by_testis/label_by_phenotype/mgi_mouse_human_homologous_gene_symbols_have_phenoes.csv'
    if not os.path.exists(mgi_mouse_gene_symbols_have_phenoes_fp):
        with open(mgi_mouse_gene_symbols_have_phenoes_fp, 'w') as output_file:
            for key in mouse_human_homolo_gene.keys():
                if key in mgi_mouse_gene_symbols_have_phenoes:
                    mgi_human_homologous_gene_symbols_have_phenoes.add(mouse_human_homolo_gene[key])
                    if key in mouse_gene_names_for_label:
                        output_file.write(key + ',' + mouse_human_homolo_gene[key] + ',' + '1' + '\n')
                    else:
                        output_file.write(key + ',' + mouse_human_homolo_gene[key] + ',' + '0' + '\n')
    print('mgi_human_homologous_gene_symbols_have_phenoes: ' + str(len(mgi_human_homologous_gene_symbols_have_phenoes)))

    gene_names_for_positive_example = merge_labeled_train_data_positive_example(dir_path)

    human_gene_names_for_label_by_mpheno = set()

    for gene_name in mouse_gene_names_for_label:
        # print(gene_name)
        if gene_name in mouse_human_homolo_gene.keys():
            if mouse_human_homolo_gene[gene_name] not in gene_names_for_positive_example:
                human_gene_names_for_label_by_mpheno.add(mouse_human_homolo_gene[gene_name])
                # print(mouse_human_homolo_gene[gene_name])
    print('The homolo gene of the mouse phenotypic: ' + str(len(human_gene_names_for_label_by_mpheno)))

    out_file_path = dir_path + 'training_set_extension_by_testis/label_by_phenotype/protein_entries_by_phenotype_labeled_total_mgi.csv'
    if not os.path.exists(out_file_path):
        with open(out_file_path, 'w') as output_file:
            for gene_name in human_gene_names_for_label_by_mpheno:
                output_file.write(gene_name + '\n')

    # 汇总所有采集的MGI 数据库上有表型的人类同源基因相关的数据集，包括下载采集的，以及直接网络爬虫采集的

    # 表型id与表型名称
    special_mp_names = {'abnormal joint capsule morphology':'MP:0000997','no suckling reflex':'MP:0001435'}
    mp_id_name = {}
    for mp_name in mp_name_id.keys():
        mp_id_name[mp_name_id[mp_name]] = mp_name
    print('mp_id_name: ' + str(len(mp_id_name)))
    human_mouse_homolo_gene = {}
    for mouse_gene in mouse_human_homolo_gene.keys():
        human_mouse_homolo_gene[mouse_human_homolo_gene[mouse_gene]] = mouse_gene
    print('human_mouse_homolo_gene: ' + str(len(human_mouse_homolo_gene)))

    mgi_mouse_gene_knockout_protein_info_fp = dir_path + 'mgi_mouse_gene_knockout_protein_info.csv'
    human_gene_phenoes_ids = {}
    gene_knockout_protein_names = set()
    with codecs.open(mgi_mouse_gene_knockout_protein_info_fp, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if len(temp) < 4:
                continue
            if temp[0].strip() not in gene_knockout_protein_names:
                human_gene_phenoes_ids[temp[0].strip()] = set()
                gene_knockout_protein_names.add(temp[0].strip())
            if temp[3].strip() not in mp_name_id.keys():
                if temp[3].strip() in special_mp_names.keys():
                    human_gene_phenoes_ids[temp[0].strip()].add(special_mp_names[temp[3].strip()])
                    continue
                print(temp[3].strip())
            mp_id = mp_name_id[temp[3].strip()]
            human_gene_phenoes_ids[temp[0].strip()].add(mp_id)
    print('gene_knockout_protein_names: '+str(len(gene_knockout_protein_names)))

    # 小鼠的基因名称与表型mgi_gene_phenoes
    # 小鼠的基因名称与对应人类的同源基因mouse_human_homolo_gene
    human_mouse_homolo_gene_mp_ids = {}
    for mouse_gene in mgi_gene_phenoes.keys():
        if mouse_gene not in mouse_human_homolo_gene.keys():
            continue
        human_mouse_homolo_gene_mp_ids[mouse_human_homolo_gene[mouse_gene] + '_' + mouse_gene] = mgi_gene_phenoes[mouse_gene]

    for human_gene in human_gene_phenoes_ids.keys():
        if human_gene not in human_mouse_homolo_gene.keys():
            continue
        if human_gene+ '_' + human_mouse_homolo_gene[human_gene] in human_mouse_homolo_gene_mp_ids.keys():
            mgi_mp_ids = set(list(human_gene_phenoes_ids[human_gene])+list(human_mouse_homolo_gene_mp_ids[human_gene+ '_' + human_mouse_homolo_gene[human_gene]]))
            human_mouse_homolo_gene_mp_ids[human_gene+ '_' + human_mouse_homolo_gene[human_gene]] = mgi_mp_ids
        else:
            human_mouse_homolo_gene_mp_ids[human_gene + '_' + human_mouse_homolo_gene[human_gene]] =human_gene_phenoes_ids[human_gene]
    print('human_mouse_homolo_gene_mp_ids: ' + str(len(human_mouse_homolo_gene_mp_ids.keys())))

    out_mgi_db_info_fp = dir_path + 'label_related_dataset_final/mgi_gene_phenoes_db_summary.csv'
    if not os.path.exists(out_mgi_db_info_fp):
        with open(out_mgi_db_info_fp, 'w') as output_file:
            for homolo_gene in human_mouse_homolo_gene_mp_ids.keys():
                if len(human_mouse_homolo_gene_mp_ids[homolo_gene]) == 0:
                    continue
                mp_ids = list(human_mouse_homolo_gene_mp_ids[homolo_gene])
                mp_names = []
                for mp_id in mp_ids:
                    mp_names.append(mp_id_name[mp_id])
                output_file.write(homolo_gene+ '$' + ';'.join(mp_ids) + '$' + ';'.join(mp_names) + '\n')

    return human_gene_names_for_label_by_mpheno

# According to the expression data of GTEX RNASEQ gene in Testis,
# the set of highly expressed genes was calculated.
# Data preparation was completed in this part and R program was completed
def preaper_testis_gene_expression_base_gtex(dir_path):

    gtex_raw_data_fp = 'D:/2_deep_profiling/vector_protein_exp/raw_data/GTEx/'
    expression_in_tissue_samples_fp = gtex_raw_data_fp + 'GTEx_v7_Annotations_SampleAttributesDS.txt'

    flag_terms = ['Blood	Whole Blood', 'Brain	Brain', 'Lung	Lung', 'Muscle	Muscle', 'Bone Marrow	Cells',
                  'Heart	Heart', 'Skin	Skin', 'Nerve	Nerve', 'Thyroid	Thyroid', 'Esophagus	Esophagus',
                  'Adipose Tissue', 'Blood Vessel', 'Blood	Cells', 'Pituitary	Pituitary', 'Testis	Testis',
                  'Pancreas	Pancreas', 'Prostate	Prostate', 'Skin	Cells']

    tissue_sample_ids = OrderedDict()
    sample_attribute_ids = set()
    with codecs.open(expression_in_tissue_samples_fp, "rb", "utf-8") as input_file:
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
    print(len(sample_attribute_ids))

    gene_official_symbol_names = prepare_human_protein_names()
    gene_official_symbols = set(gene_official_symbol_names.keys())

    gtex_expression_values_fp = gtex_raw_data_fp + 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz'
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
                print('sample_ids: ' + str(len(sample_ids)))
                continue
            temp = line.strip().split('	')
            if temp[1].strip() in gene_official_symbols:
                gene_name_values[temp[1].strip()] = [x.strip() for x in temp[2:]]
    print('gene_name_values: ' + str(len(gene_name_values.keys())))

    for tissue in tissue_sample_ids.keys():
        file_name = re.sub(' ','_','_'.join([x.strip() for x in tissue.split('-')]))
        file_name = re.sub('\(', '', file_name)
        file_name = re.sub('\)', '', file_name)
        samples_expression_from_gtex_fp = dir_path + 'training_set_extension_by_testis/testis_hyperexpression/'+file_name.lower()+'.csv'
        headers = []
        headers_index = []
        print(tissue+'_sample_ids: ' + str(len(tissue_sample_ids[tissue])))

        for sample_id in tissue_sample_ids[tissue]:
            if sample_id not in sample_ids:
                continue
            headers_index.append(sample_ids.index(sample_id))
            headers.append(sample_id)
        print(len(headers), len(headers_index))

        if len(headers) == 0:
            continue
        tissue_gene_name_values = {}
        for gene_name in gene_name_values.keys():
            temp_values = []
            for index_i in headers_index:
                temp_values.append(gene_name_values[gene_name][index_i])
            if len(temp_values) != len(headers_index):
                print(','.join(temp_values))
            tissue_gene_name_values[gene_name] = ','.join(temp_values)

        with open(samples_expression_from_gtex_fp, 'w') as output_file:
            output_file.write('gene_official_symbol'+','+','.join(headers)+'\n')
            for key in tissue_gene_name_values.keys():
                output_file.write(key+','+tissue_gene_name_values[key]+'\n')

# The list of all human genes, matching the list of genes
# with gene expression data on GTEX, and eliminating the list of human sperm specific protein genes
def temp_prepare_gene_list_for_collect_diorder_from_genecards(dir_path):

    tissue_expression_from_gtex_fp = dir_path + 'training_set_extension_by_testis/testis_hyperexpression/testis.csv'
    human_gene_from_gtex = set()
    with codecs.open(tissue_expression_from_gtex_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split(',')
            human_gene_from_gtex.add(temp[0].strip())
    print('human_gene_from_gtex: '+str(len(human_gene_from_gtex)))

    human_sperm_all_abbreviation_names = prepare_human_sperm_protein_name()
    human_sperm_genes = set()
    for key in human_sperm_all_abbreviation_names.keys():
        human_sperm_genes.add(key)
        for gene_name in human_sperm_all_abbreviation_names[key]:
            human_sperm_genes.add(gene_name)

    gene_names_from_train_data_samples = merge_phenotype_disorder_for_train_data_label(dir_path)
    human_gene_names_for_label_by_mpheno = mapping_total_potential_label_dataset_from_mgi(dir_path)

    need_collect_disorder_from_genecards = set()
    for gene_name in human_gene_from_gtex:
        if gene_name in human_sperm_genes:
            continue
        if gene_name in gene_names_from_train_data_samples:
            continue
        if gene_name in human_gene_names_for_label_by_mpheno:
            continue
        need_collect_disorder_from_genecards.add(gene_name)
    print('need_collect_disorder_from_genecards: '+str(len(need_collect_disorder_from_genecards)))

    out_file_path = dir_path + 'training_set_extension_by_testis/label_by_disorder/need_collect_disorder_from_genecards_all.csv'
    if not os.path.exists(out_file_path):
        with open(out_file_path, 'w') as output_file:
            for gene_name in need_collect_disorder_from_genecards:
                output_file.write(gene_name + '\n')
    return need_collect_disorder_from_genecards,human_gene_from_gtex,human_sperm_all_abbreviation_names

# A list of highly expressed genes obtained based on RNASEQ gene expression data
def prepare_gene_list_base_hight_expression_in_gtex(dir_path):
    files_path = dir_path + 'training_set_extension_by_testis/testis_hyperexpression/results_1211/'
    file_names = [f for f in listdir(files_path) if f.startswith('testis_')]
    need_collect_disorder_from_genecards,human_gene_from_gtex,human_sperm_all_abbreviation_names = temp_prepare_gene_list_for_collect_diorder_from_genecards(dir_path)
    need_collect_gene_name_list_by_high_expression = need_collect_disorder_from_genecards
    total_gene_name_list_by_high_expression = human_gene_from_gtex
    for file_name in file_names:
        temp_gene_name = set()
        with codecs.open(files_path+file_name, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 1, None):
                temp = line.strip().split(',')
                gene_name = re.sub('"','',temp[0]).strip()
                temp_gene_name.add(gene_name)
        total_gene_name_list_by_high_expression = total_gene_name_list_by_high_expression.intersection(temp_gene_name)
        need_collect_gene_name_list_by_high_expression = need_collect_gene_name_list_by_high_expression.intersection(temp_gene_name)
        if len(need_collect_gene_name_list_by_high_expression) == 0:
            break
    print('need_collect_gene_name_list_by_high_expression: '+str(len(need_collect_gene_name_list_by_high_expression)))
    print('total_gene_name_list_by_high_expression: ' + str(len(total_gene_name_list_by_high_expression)))

    out_file_path = dir_path + 'training_set_extension_by_testis/label_by_disorder/need_collect_disorder_from_genecards_deseq2_v2.csv'
    if (not os.path.exists(out_file_path)) and (len(need_collect_gene_name_list_by_high_expression) != 0):
        with open(out_file_path, 'w') as output_file:
            for gene_name in need_collect_gene_name_list_by_high_expression:
                output_file.write(gene_name + '\n')

    out_file_path_total = dir_path + 'training_set_extension_by_testis/total_gene_names_with_overexpressed_in_testis_by_deseq2_v2.csv'
    if (not os.path.exists(out_file_path_total)) and (len(total_gene_name_list_by_high_expression) != 0):
        with open(out_file_path_total, 'w') as output_file:
            for gene_name in total_gene_name_list_by_high_expression:
                if gene_name in human_sperm_all_abbreviation_names:
                    output_file.write(gene_name + ',' + 'hs' + '\n')
                else:
                    output_file.write(gene_name + ',' + 'ht' + '\n')

    return need_collect_gene_name_list_by_high_expression

# Collected disease data and differential expression data
# were processed and labeled according to disease keywords
def processing_disorder_differential_expression(dir_path):
    # fold change > 4, p-value < 0.05
    file_path = dir_path + 'training_set_extension_by_testis/label_by_disorder/collect_disorder_differential_expression_from_genecards_deseq2_v2.csv'
    gene_list_all = set()
    differential_expression_gene_list_by_gtex = set()
    various_labels = {}
    protein_symbol_disorder = []
    with codecs.open(file_path, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) < 3:
                continue
            if temp[1].strip() not in various_labels.keys():
                if temp[2].strip() == '':
                    various_labels[temp[1].strip().lower()] = [temp[1].strip().lower()]
                else:
                    various_labels[temp[1].strip().lower()] = [temp[1].strip().lower()] + [x.lower() for x in
                                                                                           temp[2].strip().split(';')]
            else:
                various_labels[temp[1].strip().lower()] += [x.lower() for x in temp[2].strip().split(';')]
            if [temp[0].strip(), temp[1].strip().lower()] not in protein_symbol_disorder:
                protein_symbol_disorder.append([temp[0].strip(), temp[1].strip().lower()])
            gene_list_all.add(temp[0].strip())
            if 'Testis' in temp[-1].strip():
                differential_expression_gene_list_by_gtex.add(temp[0].strip())
    print('differential_expression_gene_list_by_gtex: '+str(len(differential_expression_gene_list_by_gtex)))
    print('gene_list_all: ' + str(len(gene_list_all)))
    print(len(various_labels.keys()), len(protein_symbol_disorder))

    # 打标签的关键词
    _, disorder_keywords = phenotype_terms_disorder_keywords(dir_path)
    disorder_keywords_clean = set()
    for disorder_z in disorder_keywords:
        disorder_keywords_clean.add(disorder_z.lower())

    # 根据已经打过的标签统计的
    filter_wrong_disorders_fp = dir_path + 'label_by_disorders/filter_wrong_disorders.csv'
    filter_wrong_disorders = set()
    with codecs.open(filter_wrong_disorders_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().lower()
            filter_wrong_disorders.add(temp)
    print('filter_wrong_disorders: ' + str(len(filter_wrong_disorders)))

    # 打标签
    protein_entries_labeled_by_disorder = set()
    disorders_by_labeled = set()
    for psd in protein_symbol_disorder:
        flag = 0
        for keyword in disorder_keywords_clean:
            diseases = various_labels[psd[1]]
            if psd[1] in filter_wrong_disorders:
                continue
            if keyword in diseases:
                protein_entries_labeled_by_disorder.add(psd[0])
                disorders_by_labeled.add(psd[1])
                break
            for disease in diseases:
                disease = re.sub(',',' ',disease)
                disease = re.sub('/', ' ', disease)
                disease_words = disease.strip().split(' ')
                disease_words_temp = []
                for dw in disease_words:
                    disease_words_temp.append(dw.strip())
                if keyword in disease_words_temp:
                    protein_entries_labeled_by_disorder.add(psd[0])
                    disorders_by_labeled.add(psd[1])
                    flag = 1
                    break
                if keyword == 'sperm':
                    for disease_word in disease_words_temp:
                        if keyword in disease_word:
                            protein_entries_labeled_by_disorder.add(psd[0])
                            disorders_by_labeled.add(psd[1])
                            flag = 1
                            break
            if flag == 1:
                break
    print('genecards_protein_entries_labeled_by_disorder: ' + str(len(protein_entries_labeled_by_disorder)))

    disease_out_file_path = dir_path + 'training_set_extension_by_testis/label_by_disorder/disease_names_by_disorder_labeled.csv'
    if not os.path.exists(disease_out_file_path):
        with open(disease_out_file_path, 'w') as output_file:
            for diease in disorders_by_labeled:
                output_file.write(diease.strip() + '\n')

    out_file_path = dir_path + 'training_set_extension_by_testis/label_by_disorder/protein_entries_by_disorder_labeled_from_gtex.csv'
    if not os.path.exists(out_file_path):
        with open(out_file_path, 'w') as output_file:
            for gene_name in protein_entries_labeled_by_disorder:
                output_file.write(gene_name + '\n')
    return protein_entries_labeled_by_disorder

# Statistics the total amount of all kinds of data,
# as well as the number of all kinds of labels
def statistic_seed_gene_pool_and_train_dataset(dir_path):

    # uniprot 与 ncbi 数据库上人类蛋白质集合
    protein_official_symbol_names = prepare_human_protein_names()
    protein_official_and_alternative_symbols = set()
    for protein_official_symbol in protein_official_symbol_names.keys():
        protein_official_and_alternative_symbols.add(protein_official_symbol)
        for alternative_symbol in protein_official_symbol_names[protein_official_symbol]:
            protein_official_and_alternative_symbols.add(alternative_symbol)
    print('protein_official_and_alternative_symbols: ' + str(
        len(protein_official_and_alternative_symbols)))

    # 人类精子特定蛋白质/基因集合
    all_abbreviation_names = prepare_human_sperm_protein_name()
    human_sperm_gene_all_abbreviation_names = set()
    for key in all_abbreviation_names.keys():
        human_sperm_gene_all_abbreviation_names.add(key)
        for gene_symbol in all_abbreviation_names[key]:
            human_sperm_gene_all_abbreviation_names.add(gene_symbol)
    print('human_sperm_gene_all_abbreviation_names: ' + str(len(human_sperm_gene_all_abbreviation_names)))
    human_sperm_protein_and_human_protein_difference = set(all_abbreviation_names.keys()).difference(
        protein_official_and_alternative_symbols)

    print('human_sperm_protein_and_human_protein_difference: ' + str(
        len(human_sperm_protein_and_human_protein_difference)))
    print('-----------------------------------------')

    # 睾丸组织过表达基因列表（fold change > 4 and p-value < 0.05,GTEx 基因表达数据）
    human_testis_overexpresed_gene_symbols = dir_path + 'label_related_dataset_final/total_gene_names_with_overexpressed_in_testis_by_deseq2_v2.csv'
    testis_overexpresed_gene_symbols = set()
    testis_overexpresed_gene_symbols_hs = set()
    with codecs.open(human_testis_overexpresed_gene_symbols, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            testis_overexpresed_gene_symbols.add(temp[0].strip())
            if temp[0].strip() in human_sperm_gene_all_abbreviation_names:
                testis_overexpresed_gene_symbols_hs.add(temp[0].strip())
    print('testis_overexpresed_gene_symbols: ' + str(len(testis_overexpresed_gene_symbols)))
    print('testis_overexpresed_gene_symbols_hs: ' + str(len(testis_overexpresed_gene_symbols_hs)))

    # 只有属于人类蛋白质简称及其别称集合的睾丸高表达蛋白质集合才放入考虑范围
    testis_overexpresed_gene_symbols_in_human_proteins = testis_overexpresed_gene_symbols.intersection(
        protein_official_and_alternative_symbols)
    print('testis_overexpresed_gene_symbols_in_human_proteins: ' + str(
        len(testis_overexpresed_gene_symbols_in_human_proteins)))

    testis_overexpresed_gene_symbols_in_human_proteins_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/testis_overexpresed_gene_symbols_in_human_proteins.csv'
    if not os.path.exists(testis_overexpresed_gene_symbols_in_human_proteins_fp):
        with open(testis_overexpresed_gene_symbols_in_human_proteins_fp, 'w') as output_file:
            for gene_symbol in testis_overexpresed_gene_symbols_in_human_proteins:
                output_file.write(gene_symbol.strip() + '\n')

    print('-----------------------------------------')

    # 人类小鼠同源基因中有小鼠表型的集合
    human_mouse_homeotic_gene_symbols = dir_path + 'label_related_dataset_final/mgi_mouse_human_homologous_gene_symbols_have_phenoes.csv'
    mgi_homeotic_gene_symbols = set()
    mgi_homeotic_gene_symbols_hs = set()
    mgi_homeotic_gene_labeled_by_phenotype = set()
    with codecs.open(human_mouse_homeotic_gene_symbols, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            mgi_homeotic_gene_symbols.add(temp[1].strip())
            if temp[1].strip() in human_sperm_gene_all_abbreviation_names:
                mgi_homeotic_gene_symbols_hs.add(temp[1].strip())
            if temp[2].strip() == '1':
                mgi_homeotic_gene_labeled_by_phenotype.add(temp[1].strip())
    print('mgi_homeotic_gene_symbols: ' + str(len(mgi_homeotic_gene_symbols)))
    print('mgi_homeotic_gene_symbols_hs: ' + str(len(mgi_homeotic_gene_symbols_hs)))
    print('mgi_homeotic_gene_labeled_by_phenotype: ' + str(len(mgi_homeotic_gene_labeled_by_phenotype)))

    # 基于小鼠表型打上label的训练集，人类精子特定蛋白质结合
    file_path_phenotypes_hs = dir_path + 'label_related_dataset_final/protein_entries_by_phenotype_labeled_hs.csv'
    human_sperm_gene_labeled_by_mgi_phenotype = set()
    with codecs.open(file_path_phenotypes_hs, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('\t')
            human_sperm_gene_labeled_by_mgi_phenotype.add(temp[0].strip())
    print('human_sperm_gene_labeled_by_mgi_phenotype: ' + str(len(human_sperm_gene_labeled_by_mgi_phenotype)))

    print('-----------------------------------------')

    total_human_protein_labeled_by_mgi_phenotype = (mgi_homeotic_gene_labeled_by_phenotype.union(human_sperm_gene_labeled_by_mgi_phenotype)).intersection(protein_official_and_alternative_symbols)
    print('total_human_protein_labeled_by_mgi_phenotype: ' + str(len(total_human_protein_labeled_by_mgi_phenotype)))

    # 两类数据集上已经打上疾病label的训练数据
    file_path_disorders_hs = dir_path + 'label_related_dataset_final/protein_entries_by_disorder_labeled_hs.csv'
    human_sperm_gene_labeled_by_disorder = set()
    undetermined_gene_labeled_by_disorder = set()
    with codecs.open(file_path_disorders_hs, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('	')
            if temp[1].strip() == '1':
                human_sperm_gene_labeled_by_disorder.add(temp[0].strip())
            else:
                undetermined_gene_labeled_by_disorder.add(temp[0].strip())
    print('human_sperm_gene_labeled_by_disorder: ' + str(len(human_sperm_gene_labeled_by_disorder)))
    print('undetermined_gene_labeled_by_disorder: ' + str(len(undetermined_gene_labeled_by_disorder)))
    file_path_disorders_ht = dir_path + 'label_related_dataset_final/protein_entries_by_disorder_labeled_ht.csv'
    fragment_testis_overexpresed_gene_labeled_by_disorder = set()
    with codecs.open(file_path_disorders_ht, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            fragment_testis_overexpresed_gene_labeled_by_disorder.add(line.strip())
    print('fragment_testis_overexpresed_gene_labeled_by_disorder: ' + str(len(fragment_testis_overexpresed_gene_labeled_by_disorder)))
    gene_symbols_labeled_by_disease_keywords = (human_sperm_gene_labeled_by_disorder.
                                                union(undetermined_gene_labeled_by_disorder)).union(
        fragment_testis_overexpresed_gene_labeled_by_disorder)
    print('gene_symbols_labeled_by_disease_keywords: ' + str(len(gene_symbols_labeled_by_disease_keywords)))

    # 若干基因对应疾病数据库基因列表,包括基于关键词或者疾病名称扩展出来的训练样本
    file_path_disease_for_labeled = dir_path + 'label_related_dataset_final/disease_names_by_disorder_labeled.csv'
    disease_for_labeled = set()
    dieases_related_with_gene_for_label = set()
    with codecs.open(file_path_disease_for_labeled, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            disease_for_labeled.add(temp)
            dieases_related_with_gene_for_label.add(temp)
    print('disease_for_labeled: ' + str(len(disease_for_labeled)))
    # filter_wrong_disorders_fp = dir_path + 'label_related_dataset_final/filter_wrong_disorders.csv'
    filter_wrong_disorders_fp = dir_path + 'label_related_dataset_final/filter_wrong_disorders_20190117_add.csv'
    filter_wrong_disorders = set()
    with codecs.open(filter_wrong_disorders_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().lower()
            filter_wrong_disorders.add(temp)
    print('filter_wrong_disorders: ' + str(len(filter_wrong_disorders)))
    file_path_disease_gene = dir_path + 'label_related_dataset_final/several_gene_disease_db_summary.csv'
    gene_have_related_disease_symbols_list = set()
    gene_symbols_labeled_by_target_disease_name = set()
    with codecs.open(file_path_disease_gene, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            gene_abb = temp[1].strip().split('   ')[0].strip()
            gene_have_related_disease_symbols_list.add(gene_abb)
            if temp[2].strip() in filter_wrong_disorders:
                continue
            if temp[0].strip() != 'db_genecards':
                for stand_disease in disease_for_labeled:
                    similarity_ratio = Levenshtein.ratio(stand_disease, temp[2].strip())
                    if similarity_ratio > 0.9:
                        gene_symbols_labeled_by_target_disease_name.add(gene_abb)
                        dieases_related_with_gene_for_label.add(temp[2].strip())
                        # print(stand_disease)
                        # print(temp[2].strip())
                        break
            else:
                disease = temp[2].strip().split(';')[0].strip()
                for stand_disease in disease_for_labeled:
                    similarity_ratio = Levenshtein.ratio(stand_disease, disease)
                    if similarity_ratio > 0.9:
                        gene_symbols_labeled_by_target_disease_name.add(gene_abb)
                        dieases_related_with_gene_for_label.add(disease)
                        # print(stand_disease)
                        # print(disease)
                        break
    print('gene_symbols_labeled_by_target_disease_name: ' + str(len(gene_symbols_labeled_by_target_disease_name)))
    print('dieases_related_with_gene_for_label: ' + str(len(dieases_related_with_gene_for_label)))
    prepare_gene_symbols_for_train_model_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/dieases_related_with_gene_for_label.csv'
    if not os.path.exists(prepare_gene_symbols_for_train_model_fp):
        with open(prepare_gene_symbols_for_train_model_fp, 'w') as output_file:
            for disease in dieases_related_with_gene_for_label:
                output_file.write(disease + '\n')

    print('-----------------------------------------')

    print('gene_have_related_disease_symbols_list: ' + str(len(gene_have_related_disease_symbols_list)))
    total_gene_symbols_labeled_by_disease_name_and_keywords = (gene_symbols_labeled_by_target_disease_name.union(gene_symbols_labeled_by_disease_keywords)).intersection(protein_official_and_alternative_symbols)
    print('total_gene_symbols_labeled_by_disease_name_and_keywords: ' + str(len(total_gene_symbols_labeled_by_disease_name_and_keywords)))

    total_gene_symbols_labeled_by_disease_name_and_keywords_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/total_gene_symbols_labeled_by_disease_name_and_keywords.csv'
    if not os.path.exists(total_gene_symbols_labeled_by_disease_name_and_keywords_fp):
        with open(total_gene_symbols_labeled_by_disease_name_and_keywords_fp, 'w') as output_file:
            for disease in total_gene_symbols_labeled_by_disease_name_and_keywords:
                output_file.write(disease + '\n')

    # 若干基因与其对应的疾病数据库，为两类数据源打上label的基因统计，同时，包括遍历得到的不属于两类数据源的label
    testis_overexpresed_gene_with_disorder_label = total_gene_symbols_labeled_by_disease_name_and_keywords.intersection(testis_overexpresed_gene_symbols_in_human_proteins)
    human_sperm_gene_with_disorder_label = total_gene_symbols_labeled_by_disease_name_and_keywords.intersection(human_sperm_gene_all_abbreviation_names)
    other_gene_with_disorder_label = (total_gene_symbols_labeled_by_disease_name_and_keywords.difference(testis_overexpresed_gene_with_disorder_label.union(human_sperm_gene_with_disorder_label))).intersection(protein_official_and_alternative_symbols)
    print('testis_overexpresed_gene_with_disorder_label: ' + str(len(testis_overexpresed_gene_with_disorder_label)))
    print('human_sperm_gene_with_disorder_label: ' + str(len(human_sperm_gene_with_disorder_label)))
    print('other_gene_with_disorder_label: ' + str(len(other_gene_with_disorder_label)))

    print('-----------------------------------------')

    # 有表型数据的人类与小鼠同源基因列表
    human_homeotic_have_mgi_phenotype_gene_list_in_human_proteins = (mgi_homeotic_gene_symbols.union(
        human_sperm_gene_labeled_by_mgi_phenotype)).intersection(protein_official_and_alternative_symbols)
    total_gene_symbols_labeled_by_phenotype_name = (mgi_homeotic_gene_labeled_by_phenotype.union(
        human_sperm_gene_labeled_by_mgi_phenotype)).intersection(protein_official_and_alternative_symbols)
    print('human_homeotic_have_mgi_phenotype_gene_list_in_human_proteins: ' + str(len(human_homeotic_have_mgi_phenotype_gene_list_in_human_proteins)))
    print('total_gene_symbols_labeled_by_phenotype_name: ' + str(len(total_gene_symbols_labeled_by_phenotype_name)))

    human_homeotic_have_mgi_phenotype_gene_list_in_human_proteins_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/human_homeotic_have_mgi_phenotype_gene_list_in_human_proteins.csv'
    if not os.path.exists(human_homeotic_have_mgi_phenotype_gene_list_in_human_proteins_fp):
        with open(human_homeotic_have_mgi_phenotype_gene_list_in_human_proteins_fp, 'w') as output_file:
            for gene_symbol in human_homeotic_have_mgi_phenotype_gene_list_in_human_proteins:
                output_file.write(gene_symbol.strip() + '\n')

    # 基于小鼠同源基因表型，单独拥有的positive label,与精子特定蛋白质结合共同拥有的positive label
    testis_overexpresed_gene_with_phenotype_label = total_gene_symbols_labeled_by_phenotype_name.intersection(testis_overexpresed_gene_symbols_in_human_proteins)
    human_sperm_gene_with_phenotype_label = total_gene_symbols_labeled_by_phenotype_name.intersection(human_sperm_gene_all_abbreviation_names)
    other_gene_with_phenotype_label = (total_gene_symbols_labeled_by_phenotype_name.difference(testis_overexpresed_gene_with_phenotype_label.union(human_sperm_gene_with_phenotype_label))).intersection(protein_official_and_alternative_symbols)
    print('testis_overexpresed_gene_with_phenotype_label: ' + str(len(testis_overexpresed_gene_with_phenotype_label)))
    print('human_sperm_gene_with_phenotype_label: ' + str(len(human_sperm_gene_with_phenotype_label)))
    print('other_gene_with_phenotype_label: ' + str(len(other_gene_with_phenotype_label)))

    print('-----------------------------------------')

    total_train_dataset_with_positive_samples_by_disease_phenotype = total_gene_symbols_labeled_by_phenotype_name.union(total_gene_symbols_labeled_by_disease_name_and_keywords)
    print('total_train_dataset_with_positive_samples_by_disease_phenotype: ' + str(len(total_train_dataset_with_positive_samples_by_disease_phenotype)))

    total_train_dataset_with_positive_samples_by_disease_phenotype_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/total_train_dataset_with_positive_samples_by_disease_phenotype.csv'
    if not os.path.exists(total_train_dataset_with_positive_samples_by_disease_phenotype_fp):
        with open(total_train_dataset_with_positive_samples_by_disease_phenotype_fp, 'w') as output_file:
            for gene_symbol in total_train_dataset_with_positive_samples_by_disease_phenotype:
                output_file.write(gene_symbol.strip() + '\n')

    mgi_gene_names_for_negative_example_hs = labeled_train_data_negative_example_by_mammalian_phenotype(dir_path)
    total_train_dataset_with_negative_samples_by_phenotype = mgi_gene_names_for_negative_example_hs.difference(total_train_dataset_with_positive_samples_by_disease_phenotype)
    print('total_train_dataset_with_negative_samples_by_phenotype: ' + str(len(total_train_dataset_with_negative_samples_by_phenotype)))

    total_train_dataset_with_negative_samples_by_disease_phenotype_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/total_train_dataset_with_negative_samples_by_phenotype.csv'
    if not os.path.exists(total_train_dataset_with_negative_samples_by_disease_phenotype_fp):
        with open(total_train_dataset_with_negative_samples_by_disease_phenotype_fp, 'w') as output_file:
            for gene_symbol in total_train_dataset_with_negative_samples_by_phenotype:
                output_file.write(gene_symbol.strip() + '\n')

    # 需要优先准备的训练数据1075+2202，其他包括两类数据集剩余的部分
    prepare_gene_symbols_for_train_model_fp = dir_path + 'label_related_dataset_final/for_paper_supporting_document/prepare_gene_symbols_for_train_model.csv'
    prepare_gene_symbols_for_train_model = ((total_train_dataset_with_positive_samples_by_disease_phenotype.
        union(total_train_dataset_with_negative_samples_by_phenotype)).
        union(set(all_abbreviation_names.keys()))).union(testis_overexpresed_gene_symbols_in_human_proteins)
    print('prepare_gene_symbols_for_train_model: ' + str(
        len(prepare_gene_symbols_for_train_model)))
    if not os.path.exists(prepare_gene_symbols_for_train_model_fp):
        with open(prepare_gene_symbols_for_train_model_fp, 'w') as output_file:
            for gene_symbol in prepare_gene_symbols_for_train_model:
                output_file.write(gene_symbol.strip() + '\n')

################################################################
## 3. Determine the range of data collection and extraction
##
################################################################

def prepare_extract_feature_gene_list(dir_path):

    file_path_train_gene = dir_path + 'label_related_dataset_final/for_paper_supporting_document/prepare_gene_symbols_for_train_model.csv'
    train_gene_list = set()
    with codecs.open(file_path_train_gene, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            train_gene_list.add(temp)
    print('train_gene_list: ' + str(len(train_gene_list)))

    file_path_human_abb = dir_path + 'label_related_dataset_final/human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    human_gene_abb = {}
    human_genes = set()
    with codecs.open(file_path_human_abb, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split(',')
            human_gene_abb[temp[0].strip()] = temp[1:]
            human_genes.add(temp[0].strip())
    print('human_abb: ' + str(len(human_gene_abb.keys())))

    special_gene_symbols = {'Tex35':'TEX35','MESDC2':'MESD','C2orf47':'MAIP1'}
    file_path_train_abb = dir_path + 'label_related_dataset_final/for_paper_supporting_document/prepare_gene_symbols_abbreviation_for_train_model.csv'
    if not os.path.exists(file_path_train_abb):
        flag_genes = set()
        with open(file_path_train_abb, 'w') as output_file:
            for gene in train_gene_list:
                flag = 0
                if gene in special_gene_symbols.keys():
                    gene = special_gene_symbols[gene]
                if gene in flag_genes:
                    continue
                if gene in human_genes:
                    # if gene in flag_genes:
                    #     print('skip the gene: '+gene+'--------------------')
                    #     continue
                    output_file.write(gene + ',' + ','.join(human_gene_abb[gene]) + '\n')
                    flag_genes.add(gene)
                    flag = 1
                else:
                    for gn in human_genes:
                        if gene in human_gene_abb[gn]:
                            # if gn in flag_genes:
                            #     print('skip the gene: ' + gene+'++++++++++++++++++++++')
                            #     break
                            output_file.write(gene + ',' + ','.join(human_gene_abb[gn]) + '\n')
                            flag_genes.add(gene)
                            flag = 1
                            break
                if flag == 0:
                    print('fail to get the gene: ' + gene+'****************')
        print('The number of the final genes: ' + str(len(flag_genes)))

# To traverse and collect disease data of all human proteins
def prepare_human_protein_official_symbols_for_disorder_info(dir_path):

    total_human_protein_official_symbols_fp = 'D:/1_data_help/text_corpus/human_protein_official_symbol_ncbi_ids_uniprot_entries.csv'
    total_human_protein_official_symbols = set()
    identify_problems = ['HLA-DRB3$0301', 'HLA-DRB4$0101']
    with codecs.open(total_human_protein_official_symbols_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            for special_protein_name in identify_problems:
                if special_protein_name in line.strip():
                    total_human_protein_official_symbols.add(special_protein_name)
                    continue
            temp = line.strip().split('$')
            total_human_protein_official_symbols.add(temp[0].strip())

    human_sperm_protein_official_symbols_fp = 'C:/ggguo/1_data_help/protein_entities_standardized/human_sperm_protein_name_all_abbreviation_aliases.csv'
    human_sperm_protein_official_symbols = set()
    with codecs.open(human_sperm_protein_official_symbols_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            human_sperm_protein_official_symbols.add(temp[0].strip())

    human_protein_official_symbols_difference_by_sperm_protein_fp = dir_path + 'human_protein_official_symbols_difference_by_sperm_protein.csv'
    human_protein_official_symbols_difference_by_sperm_protein = total_human_protein_official_symbols.difference(human_sperm_protein_official_symbols)
    print('human_protein_official_symbols_difference_by_sperm_protein: ' + str(len(human_protein_official_symbols_difference_by_sperm_protein)))
    if not os.path.exists(human_protein_official_symbols_difference_by_sperm_protein_fp):
        with open(human_protein_official_symbols_difference_by_sperm_protein_fp, 'w') as output_file:
            for gene_symbol in human_protein_official_symbols_difference_by_sperm_protein:
                output_file.write(gene_symbol.strip() + '\n')

# On January 17, 2019, about 50 genes/proteins were removed due to further screening of selected disease names
def filter_wrong_selected_protein_symbols(dir_path):

    protein_offical_for_train_in_0115 = set()
    prepare_gene_symbols_abbreviation_for_train_model_0115_pf = dir_path + 'label_related_dataset_final/for_paper_supporting_document/prepare_gene_symbols_for_train_model_0115.csv'
    with codecs.open(prepare_gene_symbols_abbreviation_for_train_model_0115_pf, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            protein_offical_for_train_in_0115.add(temp[0].strip())

    protein_offical_for_train = set()
    prepare_gene_symbols_abbreviation_for_train_model_pf = dir_path + 'label_related_dataset_final/for_paper_supporting_document/prepare_gene_symbols_for_train_model.csv'
    with codecs.open(prepare_gene_symbols_abbreviation_for_train_model_pf, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            protein_offical_for_train.add(temp[0].strip())

    wrong_selected_protein_symbols = protein_offical_for_train_in_0115.difference(protein_offical_for_train)

    print(len(wrong_selected_protein_symbols))
    wrong_selected_protein_symbols_pf = dir_path + 'label_related_dataset_final/for_paper_supporting_document/wrong_selected_protein_symbols_for_train.csv'
    with open(wrong_selected_protein_symbols_pf, 'w') as output_file:
        for gene_symbol in wrong_selected_protein_symbols:
            output_file.write(gene_symbol.strip() + '\n')


if __name__ == '__main__':

    # prepare_human_protein_names()
    dir_path = ''
    phenotype_terms_disorder_keywords(dir_path)
    # labeled_train_data_positive_example_by_mammalian_phenotype(dir_path)
    # labeled_train_data_positive_example_by_human_disease(dir_path)
    # merge_phenotype_disorder_for_train_data_label(dir_path)
    # temp_match_train_dataset_total_dataset(dir_path)
    # mapping_total_potential_label_dataset_from_mgi(dir_path)
    # preaper_testis_gene_expression_base_gtex(dir_path)
    # temp_prepare_gene_list_for_collect_diorder_from_genecards(dir_path)
    # prepare_gene_list_base_hight_expression_in_gtex(dir_path)
    # processing_disorder_differential_expression(dir_path)
    # statistic_seed_gene_pool_and_train_dataset(dir_path)
    # prepare_extract_feature_gene_list(dir_path)
    # prepare_human_protein_official_symbols_for_disorder_info(dir_path)
    # filter_wrong_selected_protein_symbols(dir_path)



