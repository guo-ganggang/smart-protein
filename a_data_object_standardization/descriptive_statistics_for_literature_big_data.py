#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 25/10/2018 11:09 AM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : descriptive_statistics_for_literature_big_data.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import codecs
from itertools import islice
from os import listdir
import os
import pandas as pd
import datetime
from collections import OrderedDict
import xlrd
import re

import sys

reload(sys)
sys.setdefaultencoding('utf8')

################################################################
## Levenshtein Distance: It's a kind of edit distance.
##
################################################################
def levenshtein(str1, str2):
    m, n = len(str1) + 1, len(str2) + 1

    # 初始化矩阵
    matrix = [[0] * n for i in range(m)]
    matrix[0][0] = 0
    for i in range(1, m):
        matrix[i][0] = matrix[i - 1][0] + 1
        for j in range(1, n):
            matrix[0][j] = matrix[0][j - 1] + 1
    # 动态规划计算ld值
    for i in range(1, m):
        for j in range(1, n):
            if str1[i - 1] == str2[j - 1]:
                matrix[i][j] = matrix[i - 1][j - 1]
            else:
                matrix[i][j] = min(matrix[i - 1][j - 1], matrix[i - 1][j], matrix[i][j - 1]) + 1

    return matrix[m - 1][j - 1]

################################################################
## Descriptive statistics of literature data -- date,
# number of publications, number of references associated with a single protein,
# and quality of associated protein in a single literature
##
################################################################

# Get the name of all the files in the folder, parse out the ID,
# and the literature PMID and publication date
def read_parse_dataset(dir_path):
    id_publish_date = {}
    fileNameFeature = '.txt'
    fileNames = [f for f in listdir(dir_path) if f.endswith(fileNameFeature)]
    print(len(fileNames))
    for fileName in fileNames:
        temp_publish_date = []
        with codecs.open(dir_path+fileName, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                if len(temp) < 2:
                    continue
                if (not temp[0].strip().isdigit()) or (not temp[1].strip().isdigit()):
                    continue
                temp_publish_date.append(temp[:2])
        if len(temp_publish_date) == 0:
            continue
        id = fileName.strip().split('.')[0].strip().split('_')[-1].strip()
        id_publish_date[id] = temp_publish_date
    print(len(id_publish_date))

    # 返回 蛋白质/基因 ID，以及pmid和发表日期
    return id_publish_date

# Matches protein names and IDs from both databases
def id_protein_entry_match(dir_path):
    dir_path_ncbi = dir_path + 'ncbi/gene_symbols_id_human.csv'
    dir_path_ncbi_gene_entry_unpicking = dir_path + 'ncbi/protein_ncbi_entry_unpicking.csv'
    dir_path_uniprot = dir_path + 'uniprot/uniprot_gene_entry_clean_data.csv'
    dir_path_uniprot_gene_entry_unpicking = dir_path + 'uniprot/protein_uniprot_entry_unpicking.csv'


    id_protein_entry_ncbi = {}
    id_protein_entry_uniprot = {}
    with codecs.open(dir_path_ncbi, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            id_protein_entry_ncbi[temp[1].strip()] = temp[0].strip()

    with codecs.open(dir_path_ncbi_gene_entry_unpicking, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) <= 2:
                continue
            if temp[-1].strip() == 'null':
                continue
            if temp[-1].strip() not in id_protein_entry_ncbi.keys():
                id_protein_entry_ncbi[temp[-1].strip()] = temp[0].strip()
    print(len(id_protein_entry_ncbi.keys()))

    count = set()
    with codecs.open(dir_path_uniprot, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if temp[0].strip() == 'HUMAN':
                if temp[3].strip() in count:
                    print(temp[3].strip() + '-----------')
                    continue
                count.add(temp[3].strip())
                id_protein_entry_uniprot[temp[3].strip()] = temp[1].strip()

    with codecs.open(dir_path_uniprot_gene_entry_unpicking, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) <= 2:
                continue
            if temp[-1].strip() == 'null':
                continue
            if temp[-1].strip() not in id_protein_entry_uniprot.keys():
                id_protein_entry_uniprot[temp[-1].strip()] = temp[0].strip()
    print(len(id_protein_entry_uniprot.keys()))

    return id_protein_entry_ncbi, id_protein_entry_uniprot

# The number of literatures in the two databases and
# the literature PMID were counted by date. In addition,
# the publication date collection of statistical literatures by protein was included
def statistic_literature_by_date(dir_path,by_date='m'):
    id_pmid_publish_date = read_parse_dataset(dir_path)
    pmid_publish_dates = {}
    count_num = []
    id_publish_dates = {}
    for id in id_pmid_publish_date.keys():
        id_publish_dates[id] = []
        for pmid_date in id_pmid_publish_date[id]:
            year = pmid_date[1][:4]
            month = pmid_date[1][4:6]
            day = pmid_date[1][6:8]
            pmid_publish_dates[pmid_date[0]]=year+ '-' + month + '-' +day
            count_num.append(pmid_date[0])
            if by_date == 'm':
                if (year + '-' + month) not in id_publish_dates[id]:
                    id_publish_dates[id].append(year + '-' + month)
            else:
                if year not in id_publish_dates[id]:
                    id_publish_dates[id].append(year)
    print(len(pmid_publish_dates.keys()))
    print(len(count_num))

    # 获取所有文献最早日期与最晚日期
    df_date = pd.to_datetime(pmid_publish_dates.values(),format="%Y-%m-%d")
    min_date = datetime.datetime.now().date()
    max_date = (datetime.datetime.now() + datetime.timedelta(days=-1000)).date()
    print(min_date,max_date)
    for publish_date in df_date:
        # print(publish_date.date())
        days_differ_min = (publish_date.date() - min_date).days
        days_differ_max = (publish_date.date() - max_date).days
        if days_differ_min < 0:
            min_date = publish_date.date()
        if days_differ_max > 0:
            max_date = publish_date.date()
    print(min_date,max_date)

    p_date_count = OrderedDict()
    p_date_pmids = OrderedDict()

    if by_date == 'm':
        switch = 0
        for p_year in range(int(min_date.year),int(max_date.year)+1):
            for p_month in range(1,13):
                p_year_month = str(p_year) + '-' + str(p_month)
                if p_year_month == (str(min_date.year) + '-' + str(min_date.month)):
                    switch = 1
                if switch != 1:
                    continue
                if p_year_month == (str(max_date.year) + '-' + str(max_date.month)):
                    p_date_count[p_year_month] = 0
                    break
                else:
                    p_date_count[p_year_month] = 0
        print(len(p_date_count))

        for year_month_day in pmid_publish_dates.values():
            year_month = '-'.join([str(int(x)) for x in year_month_day.split('-')[:2]])
            p_date_count[year_month] += 1
    else:
        for p_year in range(int(min_date.year), int(max_date.year) + 1):
            p_date_count[str(p_year)] = 0
            p_date_pmids[str(p_year)] = []
        for key in pmid_publish_dates.keys():
            year_month_day = pmid_publish_dates[key]
            p_date_count[year_month_day.split('-')[0].strip()] += 1
            if key not in p_date_pmids[year_month_day.split('-')[0].strip()]:
                p_date_pmids[year_month_day.split('-')[0].strip()].append(key)

    return p_date_count,p_date_pmids,id_publish_dates

# The literature quantity and literature PMID of the two databases
# were counted according to protein, and related data were combined.
# In addition, protein collection was counted according to literature
def statistic_literature_by_protein(dir_path):

    dir_path_ncbi = dir_path + 'ncbi/parse_literatures_human/'
    dir_path_uniprot = dir_path + 'uniprot/parse_literatures_human/'

    protein_name_pmids_total = {}
    total_protein_name_pmids_file_path = dir_path + 'descriptive_statistics/protein_name_pmids_total.csv'
    if not os.path.exists(total_protein_name_pmids_file_path):
        # 按蛋白质合并两个数据库

        id_pmid_publish_date_ncbi = read_parse_dataset(dir_path_ncbi)
        id_pmid_publish_date_uniprot = read_parse_dataset(dir_path_uniprot)

        # 字典类型,key转换成蛋白质名称，value 为文献pmid
        id_protein_entry_ncbi, id_protein_entry_uniprot = id_protein_entry_match(dir_path)
        protein_name_pmids_ncbi = {}
        protein_name_pmids_uniprot = {}
        for id in id_pmid_publish_date_ncbi.keys():
            protein_name_pmids_ncbi[id_protein_entry_ncbi[id]] = []
            for pmid_date in id_pmid_publish_date_ncbi[id]:
                protein_name_pmids_ncbi[id_protein_entry_ncbi[id]].append(pmid_date[0])
        for id in id_pmid_publish_date_uniprot.keys():
            protein_name_pmids_uniprot[id_protein_entry_uniprot[id]] = []
            for pmid_date in id_pmid_publish_date_uniprot[id]:
                protein_name_pmids_uniprot[id_protein_entry_uniprot[id]].append(pmid_date[0])

        # 两个数据库文献合并成新的字典类型
        for protein_name in protein_name_pmids_ncbi.keys():
            protein_name_pmids_total[protein_name] = []
            for pmid in protein_name_pmids_ncbi[protein_name]:
                if pmid not in protein_name_pmids_total[protein_name]:
                    protein_name_pmids_total[protein_name].append(pmid)
            if protein_name in protein_name_pmids_uniprot.keys():
                for pmid in protein_name_pmids_uniprot[protein_name]:
                    if pmid not in protein_name_pmids_total[protein_name]:
                        protein_name_pmids_total[protein_name].append(pmid)
        for protein_name in protein_name_pmids_uniprot.keys():
            if protein_name not in protein_name_pmids_total.keys():
                protein_name_pmids_total[protein_name] = []
                for pmid in protein_name_pmids_uniprot[protein_name]:
                    if pmid not in protein_name_pmids_total[protein_name]:
                        protein_name_pmids_total[protein_name].append(pmid)
        print(len(protein_name_pmids_total))

        with open(total_protein_name_pmids_file_path, 'w') as output_file:
            for protein_name in protein_name_pmids_total.keys():
                output_file.write(protein_name + ',' + ','.join(protein_name_pmids_total[protein_name]) + '\n')
    else:
        with codecs.open(total_protein_name_pmids_file_path, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                protein_name_pmids_total[temp[0].strip()]=temp[1:]
        print(len(protein_name_pmids_total))

    total_pmid_date_file_path = dir_path + 'descriptive_statistics/pmid_publish_dates_total.csv'
    pmid_publish_dates_total = {}
    if not os.path.exists(total_pmid_date_file_path):
        # 获取全部pmid与发布日期
        id_pmid_publish_date_ncbi = read_parse_dataset(dir_path_ncbi)
        id_pmid_publish_date_uniprot = read_parse_dataset(dir_path_uniprot)
        pmids = set()
        for id in id_pmid_publish_date_ncbi.keys():
            for pmid_date in id_pmid_publish_date_ncbi[id]:
                if pmid_date[0] not in pmids:
                    year = pmid_date[1][:4]
                    month = pmid_date[1][4:6]
                    day = pmid_date[1][6:8]
                    pmid_publish_dates_total[pmid_date[0]]=year+ '-' + month + '-' +day
                    pmids.add(pmid_date[0])
        for id in id_pmid_publish_date_uniprot.keys():
            for pmid_date in id_pmid_publish_date_uniprot[id]:
                if pmid_date[0] not in pmids:
                    year = pmid_date[1][:4]
                    month = pmid_date[1][4:6]
                    day = pmid_date[1][6:8]
                    pmid_publish_dates_total[pmid_date[0]]=year+ '-' + month + '-' +day
                    pmids.add(pmid_date[0])
        print(len(pmid_publish_dates_total.keys()))

        with open(total_pmid_date_file_path, 'w') as output_file:
            for pmid in pmid_publish_dates_total.keys():
                output_file.write(pmid + ',' +pmid_publish_dates_total[pmid] + '\n')
    else:
        with codecs.open(total_pmid_date_file_path, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                pmid_publish_dates_total[temp[0].strip()]=temp[1].strip()
        print(len(pmid_publish_dates_total))

    return protein_name_pmids_total,pmid_publish_dates_total

# Perform statistics and prepare the data in the appropriate format for visual presentation
def exec_by_date_statistic_for_visualization(dir_path):

    '''
        按 月/年 统计NCBI与Uniprot上 6871 个蛋白质相关文献（人类、小鼠）
    '''

    dir_path_ncbi = dir_path + 'ncbi/parse_literatures_human/'
    dir_path_uniprot = dir_path + 'uniprot/parse_literatures_human/'
    p_date_count_ncbi, p_date_pmids_ncbi, id_publish_dates_ncbi = statistic_literature_by_date(dir_path_ncbi, 'y')
    p_date_count_uniprot, p_date_pmids_uniprot, id_publish_dates_uniprot = statistic_literature_by_date(
        dir_path_uniprot, 'y')

    # 双气泡图可视化方式，呈现两个文献数据库上历年文献数量，以及某年交集文献量占两个库当年并集文献总量的比率
    print(p_date_count_ncbi.keys())
    print(p_date_count_ncbi.values())
    print(p_date_count_uniprot.keys())
    print(p_date_count_uniprot.values())

    common_pmids_num = []
    for year in range(1966, 2019):
        intersection_pmids = list(set(p_date_pmids_ncbi[str(year)]).intersection(set(p_date_pmids_uniprot[str(year)])))
        common_pmids_num.append(round(len(intersection_pmids) / (float(p_date_count_ncbi[str(year)])+float(p_date_count_uniprot[str(year)])),4))
    print(common_pmids_num)

    # 双柱状图可视化方式，呈现两个文献数据库上历年文献数量，以及某年交集文献量分别占两个库当年总量的比率
    data_ncbi = []
    for year in range(1966,2019):
        temp = []
        intersection_pmids = list(set(p_date_pmids_ncbi[str(year)]).intersection(set(p_date_pmids_uniprot[str(year)])))
        temp.append(str(year))
        temp.append(round(len(intersection_pmids) / float(p_date_count_ncbi[str(year)]),4))
        temp.append(p_date_count_ncbi[str(year)])
        temp.append(str(year))
        temp.append('ncbi')
        data_ncbi.append(temp)
    print(p_date_count_ncbi.keys())

    data_uniprot = []
    for year in range(1966,2019):
        temp = []
        intersection_pmids = list(set(p_date_pmids_ncbi[str(year)]).intersection(set(p_date_pmids_uniprot[str(year)])))
        temp.append(str(year))
        temp.append(round(len(intersection_pmids) / float(p_date_count_uniprot[str(year)]),4))
        temp.append(p_date_count_uniprot[str(year)])
        temp.append(str(year))
        temp.append('uniprot')
        data_uniprot.append(temp)
    print(p_date_count_uniprot.keys())

    total_data = []
    total_data.append(data_ncbi)
    total_data.append(data_uniprot)
    print(total_data)

    # 冰山折线图可视化方式，呈现两个文献数据库上历年研究的蛋白质数量，以及当年新增蛋白质数据占总量的比重
    id_protein_entry_ncbi, id_protein_entry_uniprot = id_protein_entry_match(dir_path)
    year_protein_entry = OrderedDict()
    for year in range(1965, 2019):
        year_protein_entry[str(year)]=[]
        for nid in id_publish_dates_ncbi.keys():
            if str(year) in id_publish_dates_ncbi[nid]:
                year_protein_entry[str(year)].append(id_protein_entry_ncbi[nid])
        if year == 1965:
            continue
        for uid in id_publish_dates_uniprot.keys():
            if str(year) in id_publish_dates_uniprot[uid]:
                if id_protein_entry_uniprot[uid] not in year_protein_entry[str(year)]:
                    year_protein_entry[str(year)].append(id_protein_entry_uniprot[uid])

    past_protein_entry = set()
    year_count_protein_entry = []
    # {
    #     year: "<=215USD",
    #     value: 813,
    #     ratio: 0.54
    # },
    for key in year_protein_entry.keys():
        count = 0
        temp_dict = {}
        for protein in year_protein_entry[key]:
            if protein not in past_protein_entry:
                count += 1
                past_protein_entry.add(protein)

        temp_dict['year'] = key
        temp_dict['value'] = len(year_protein_entry[key])
        temp_dict['ratio'] = round(float(count) / len(year_protein_entry[key]),4)
        year_count_protein_entry.append(temp_dict)
    print(year_count_protein_entry)

# Perform statistics and prepare the data in the appropriate format for visual presentation
def exec_by_protein_statistic_for_visualization(dir_path):

    '''
        按 蛋白质 统计NCBI与Uniprot上 6871 个蛋白质相关文献（人类、小鼠）
    '''

    protein_name_pmids_total, pmid_publish_dates_total = statistic_literature_by_protein(dir_path)

    # 筛选出文献量大于2000的蛋白质，并按年统计，定义为“蛋白质研究热度”
    if not os.path.exists(dir_path+'descriptive_statistics/show_protein_pmids_mt2000_years_count.txt'):
        protein_pmids_mt2000_years_count = OrderedDict()
        for protein_name in protein_name_pmids_total.keys():
            if len(protein_name_pmids_total[protein_name]) > 2000:
                print(protein_name)
                publish_date_count = OrderedDict()
                for year in range(1970,2019):
                    publish_date_count[str(year)] = 0
                for pmid in protein_name_pmids_total[protein_name]:
                    publish_date_count[pmid_publish_dates_total[pmid].strip().split('-')[0].strip()] += 1
                protein_pmids_mt2000_years_count[protein_name] = publish_date_count.values()

        years = []
        for year in range(1970,2019):
            years.append(str(year))
        print(years)
        print([str(x) for x in protein_pmids_mt2000_years_count.keys()])
        show_protein_pmids_mt2000_years_count = []
        for j in range(len(protein_pmids_mt2000_years_count.values()[0])):
            for i in range(len(protein_pmids_mt2000_years_count.keys())):
                value = protein_pmids_mt2000_years_count[protein_pmids_mt2000_years_count.keys()[i]][j]
                show_protein_pmids_mt2000_years_count.append([j,i,value])
        # print(show_protein_pmids_mt1000_years_count)

        # 可视化格式的输出内容太多，写入文件后复制
        file_object = open(dir_path+'descriptive_statistics/show_protein_pmids_mt2000_years_count.txt', 'w')
        file_object.writelines(str(show_protein_pmids_mt2000_years_count))
        file_object.close()

    # 计算文献关联的蛋白质的数量，并计算每个蛋白质文献集合中，基于关联蛋白质的不同数量计算的类型比率，定义为“蛋白质研究黏度”
    # 比如 A 蛋白质 关联t篇文献，每篇文献关联k个蛋白质
    pmid_protein_set = {}
    pmid_protein_count = {}
    if not os.path.exists(dir_path+'descriptive_statistics/pmid_protein_set.csv'):
        for pmid in pmid_publish_dates_total.keys():
            pmid_protein_set[pmid] = set()

        for protein_name in protein_name_pmids_total.keys():
            for pmid in protein_name_pmids_total[protein_name]:
                pmid_protein_set[pmid].add(protein_name)
        print(len(pmid_protein_set.keys()))

        with open(dir_path + 'descriptive_statistics/pmid_protein_set.csv', 'a') as output_file:
            for pmid in pmid_protein_set.keys():
                temp = [pmid] + list(pmid_protein_set[pmid])
                output_file.write('$'.join(temp) + '\n')
                pmid_protein_count[temp[0]] = len(temp[1:])
    else:
        with codecs.open(dir_path + 'descriptive_statistics/pmid_protein_set.csv', "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                pmid_protein_count[temp[0]] = len(temp[1:])
                pmid_protein_set[temp[0]] = set(temp[1:])
    print(len(pmid_protein_count.keys()))

    # pmid_protein_count_sort = sorted(pmid_protein_count.items(), key=lambda e: e[1], reverse=True)
    # # intersection_set = set()
    # for pmid_count in pmid_protein_count_sort:
    #     pmid = pmid_count[0]
    #     proteins = pmid_protein_set[pmid]
    #     if len(proteins) >= 2:
    #         print(len(proteins))
    #     # if len(intersection_set) == 0:
    #     #     intersection_set = proteins
    #     #     continue
    #     # intersection_set = intersection_set.intersection(proteins)
    #     # if len(intersection_set) == 0:
    #     #     break
    #     # print('intersection..........'+str(len(intersection_set)))
    # # print(len(intersection_set))

    proteins_associated_pmid_statistic_categories = ["Ratio of associated proteins\'number =1",
                                                     'Ratio of associated proteins\'number =2',
                                                     'Ratio of associated proteins\'number >=3 and <10',
                                                     'Ratio of associated proteins\'number >=10 and <100',
                                                     'Ratio of associated proteins\'number >100']
    category_protein_ratio = OrderedDict()
    for c in range(len(proteins_associated_pmid_statistic_categories)):
        category_protein_ratio[proteins_associated_pmid_statistic_categories[c]] = []
    protein_names_pmids_count_mt_2000 = []
    for protein_name in protein_name_pmids_total.keys():
        if len(protein_name_pmids_total[protein_name]) > 2000:
            protein_names_pmids_count_mt_2000.append(protein_name)
            temp_dict = OrderedDict()
            for c in range(len(proteins_associated_pmid_statistic_categories)):
                temp_dict[proteins_associated_pmid_statistic_categories[c]] = 0
            for pmid in protein_name_pmids_total[protein_name]:
                if pmid_protein_count[pmid] == 1:
                    temp_dict[proteins_associated_pmid_statistic_categories[0]] += 1
                elif pmid_protein_count[pmid] == 2:
                    temp_dict[proteins_associated_pmid_statistic_categories[1]] += 1
                elif (pmid_protein_count[pmid] >= 3) and (pmid_protein_count[pmid] < 10):
                    temp_dict[proteins_associated_pmid_statistic_categories[2]] += 1
                elif (pmid_protein_count[pmid] >= 10) and (pmid_protein_count[pmid] < 100):
                    temp_dict[proteins_associated_pmid_statistic_categories[3]] += 1
                else:
                    temp_dict[proteins_associated_pmid_statistic_categories[4]] += 1
            print(temp_dict.values())
            for key in temp_dict.keys():
                category_protein_ratio[key].append(round(temp_dict[key] / float(len(protein_name_pmids_total[protein_name])),6))
    print([str(x) for x in protein_names_pmids_count_mt_2000])
    for category in category_protein_ratio.keys():
        print(category_protein_ratio[category])


################################################################
## Descriptive statistics of literature data -- date,
# number of references associated with a single protein,
# total amount of impact factors
##
################################################################

def journal_name_match(dir_path,colnameindex = 1):

    pubmed_journal_total = set()
    journal_total_file_path = dir_path + 'descriptive_statistics/pubmed_journal_total.csv'
    # journal_total_file_path = dir_path + 'descriptive_statistics/matched_journal_name_left.csv'
    with codecs.open(journal_total_file_path, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            pubmed_journal_total.add(temp)
            # temp = line.strip().split('$$$')
            # pubmed_journal_total.add(temp[3].strip())
    print(len(pubmed_journal_total))

    scopus_journal_name_citescore_file_path = dir_path + 'descriptive_statistics/CiteScore_Metrics_2011-2017_Download_25May2018.xlsx'

    try:
        data = xlrd.open_workbook(scopus_journal_name_citescore_file_path)
    except Exception, e:
        print str(e)

    table_sheets = ['2017 All','2016 All','2015 All','2014 All','2013 All','2012 All','2011 All']
    case_sensitives = [0,1]
    pubmed_journal_total_matched = set()
    matched_journal_name_last_time = {}
    for case_sensitive in case_sensitives:
        if case_sensitive == 0:
            critical_value = 0.85
        else:
            critical_value = 0.88
        for sheet in table_sheets:
            table = data.sheet_by_name(sheet)  # 通过名称获取
            # ncols = table.ncols
            # print(ncols)
            header_name = [u'Title',u'CiteScore']#'Scopus Sub-Subject Area','Title'
            colnames = table.row_values(colnameindex)
            print(colnames)
            titles = table.col_values(colnames.index(header_name[0]))
            print(len(titles))
            sci_mago_journal_rank = table.col_values(colnames.index(header_name[1]))
            print(len(sci_mago_journal_rank))
            titles_sci_mago_journal_rank = {}
            unique_titles_set = set()
            if len(titles) == len(sci_mago_journal_rank):
                for t in range(1,len(titles)):
                    title = titles[t].strip()
                    if title in unique_titles_set:
                        continue
                    value = str(sci_mago_journal_rank[t]).strip()
                    titles_sci_mago_journal_rank[title] = value
                    unique_titles_set.add(title)
            print(len(titles_sci_mago_journal_rank.keys()))

            matched_journal_name = {}
            for pubmed_journal_name in pubmed_journal_total:
                if pubmed_journal_name in pubmed_journal_total_matched:
                    continue
                max_levenshtein_ratio = 0.0
                max_levenshtein_str = ''
                clean_pubmed_journal_name_et = pubmed_journal_name.strip().split(' = ')
                clean_pubmed_journal_name_et.append(pubmed_journal_name.strip()) # 分割前的也包括在内
                for c in range(len(clean_pubmed_journal_name_et)):
                    clean_pubmed_journal_name = clean_pubmed_journal_name_et[c]
                    if c != (len(clean_pubmed_journal_name_et)-1):
                        clean_pubmed_journal_name = clean_pubmed_journal_name.strip().split(' : ')[0].strip()
                        clean_pubmed_journal_name = re.sub('The ', '', clean_pubmed_journal_name)
                        clean_pubmed_journal_name = re.sub(' & ', ' and ', clean_pubmed_journal_name)
                        clean_pubmed_journal_name = clean_pubmed_journal_name.strip().split('(')[0].strip()
                    max_levenshtein_ratio_i = 0.0
                    max_levenshtein_str_i = ''
                    for scopus_journal_name in titles_sci_mago_journal_rank.keys():
                        if case_sensitive == 0:
                            similarity_ratio = levenshtein(scopus_journal_name, clean_pubmed_journal_name)
                        else:
                            similarity_ratio = levenshtein(scopus_journal_name.lower(), clean_pubmed_journal_name.lower())
                        if similarity_ratio > max_levenshtein_ratio_i:
                            max_levenshtein_ratio_i = similarity_ratio
                            max_levenshtein_str_i = scopus_journal_name
                        if max_levenshtein_ratio_i == 1.0:
                            break
                    if max_levenshtein_ratio_i > max_levenshtein_ratio:
                        max_levenshtein_ratio = max_levenshtein_ratio_i
                        max_levenshtein_str = max_levenshtein_str_i
                    if max_levenshtein_ratio == 1.0:
                        break
                matched_journal_name[pubmed_journal_name] = [str(max_levenshtein_ratio),
                                                                 titles_sci_mago_journal_rank[max_levenshtein_str],
                                                                 max_levenshtein_str]
            print(len(matched_journal_name))

            matched_journal_name_fp = dir_path + 'descriptive_statistics/matched_journal_name_completed.csv'  # mte70lt85,lt70,mte85,lt85
            with open(matched_journal_name_fp, 'a') as output_file:
                for pubmed_journal_name in matched_journal_name.keys():
                    if float(matched_journal_name[pubmed_journal_name][0]) >= critical_value:
                        output_file.write('$$$'.join(matched_journal_name[pubmed_journal_name]) + '$$$' + pubmed_journal_name + '\n')
                        pubmed_journal_total_matched.add(pubmed_journal_name)
            if (case_sensitive == case_sensitives[-1]) and (sheet == table_sheets[-1]):
                matched_journal_name_last_time = matched_journal_name
    print(len(pubmed_journal_total_matched))

    matched_journal_name_left_fp = dir_path + 'descriptive_statistics/matched_journal_name_left.csv'  # mte70lt85,lt70,mte85,lt85
    with open(matched_journal_name_left_fp, 'w') as output_file:
        for pubmed_journal_name in matched_journal_name_last_time.keys():
            if pubmed_journal_name not in pubmed_journal_total_matched:
                output_file.write('$$$'.join(matched_journal_name_last_time[pubmed_journal_name]) + '$$$' + pubmed_journal_name + '\n')

def read_parse_dataset_for_influence(dir_path):
    id_pmid_date_country_journal = {}
    fileNameFeature = '.txt'
    fileNames = [f for f in listdir(dir_path) if f.endswith(fileNameFeature)]
    print(len(fileNames))
    for fileName in fileNames:
        temp_country_journal = []
        with codecs.open(dir_path+fileName, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                if len(temp) < 4:
                    continue
                if (not temp[0].strip().isdigit()) or (not temp[1].strip().isdigit()):
                    continue
                temp_country_journal.append(temp)
        if len(temp_country_journal) == 0:
            continue
        id = fileName.strip().split('.')[0].strip().split('_')[-1].strip()
        id_pmid_date_country_journal[id] = temp_country_journal
    print(len(id_pmid_date_country_journal))

    # 返回 蛋白质/基因 ID，以及pmid和发表日期
    return id_pmid_date_country_journal

def statistic_literature_by_protein_for_influence(dir_path):

    total_pmid_influence_file_path = dir_path + 'descriptive_statistics/pmid_influence_total.csv'
    pmid_influence_total = {}
    journal_name_set = set()
    country_set = set()

    if not os.path.exists(total_pmid_influence_file_path):
        # 获取全部pmid与发布日期
        dir_path_ncbi_ds = dir_path + 'descriptive_statistics/parse_literatures_human_for_ncbi_ds/'
        dir_path_uniprot_ds = dir_path + 'descriptive_statistics/parse_literatures_human_for_uniprot_ds/'
        id_pmid_date_country_journal_ncbi_ds = read_parse_dataset_for_influence(dir_path_ncbi_ds)
        id_pmid_date_country_journal_uniprot_ds = read_parse_dataset_for_influence(dir_path_uniprot_ds)
        pmids = set()
        for id in id_pmid_date_country_journal_ncbi_ds.keys():
            for pmid_influence in id_pmid_date_country_journal_ncbi_ds[id]:
                if pmid_influence[0].strip() not in pmids:
                    pmid_influence_total[pmid_influence[0].strip()] = pmid_influence[2:]
                    pmids.add(pmid_influence[0].strip())
        for id in id_pmid_date_country_journal_uniprot_ds.keys():
            for pmid_influence in id_pmid_date_country_journal_uniprot_ds[id]:
                if pmid_influence[0].strip() not in pmids:
                    pmid_influence_total[pmid_influence[0].strip()] = pmid_influence[2:]
                    pmids.add(pmid_influence[0].strip())
        print(len(pmid_influence_total.keys()))

        with open(total_pmid_influence_file_path, 'w') as output_file:
            for pmid in pmid_influence_total.keys():
                country_set.add(pmid_influence_total[pmid][0].strip())
                journal_name_set.add(pmid_influence_total[pmid][1].strip())
                output_file.write(pmid + ',' + ','.join(pmid_influence_total[pmid]) + '\n')
    else:
        with codecs.open(total_pmid_influence_file_path, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                if len(temp) < 3:
                    continue
                country_set.add(temp[1].strip())
                journal_name_set.add(temp[2].strip())
                pmid_influence_total[temp[0].strip()] = temp[1:]
        print(len(pmid_influence_total))

    journal_total_file_path = dir_path + 'descriptive_statistics/pubmed_journal_total.csv'
    with open(journal_total_file_path, 'w') as output_file:
        for journal in journal_name_set:
            output_file.write(journal.strip() + '\n')

    protein_name_pmids_total, pmid_publish_dates_total = statistic_literature_by_protein(dir_path)

    print(len(journal_name_set))
    print(len(country_set))

def exec_by_influence_statistic_for_visualization(dir_path):

    '''
        Statistics of 6871 proteins-related literatures on NCBI and Uniprot by protein (human and mouse)
    '''

    statistic_literature_by_protein_for_influence(dir_path)


if __name__ == '__main__':

    dir_path = 'data_help/'
    # exec_by_date_statistic_for_visualization(dir_path)
    # exec_by_protein_statistic_for_visualization(dir_path)
    # exec_by_influence_statistic_for_visualization(dir_path)

    journal_name_match(dir_path)