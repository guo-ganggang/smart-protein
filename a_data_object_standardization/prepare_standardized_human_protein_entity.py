#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 23/11/2018 5:10 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : prepare_standardized_human_protein_entity.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import gzip
import shutil
import xlrd
import codecs,re,os
from itertools import islice

'''
    The specific set of proteins in human sperm is named by alias
'''

# compress
def un_gz(file_name):
    """ungz zip file"""
    f_name = file_name.replace(".gz", "")
    with gzip.open(file_name, 'rb') as f_in:
        with open(f_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return f_name

# Preparation of Human Sperm Specific Protein Collection abbreviations and their acronyms from three databases, genecards,uniprot,ncbi
def prepare_human_sperm_protein_name():

    human_sperm_protein_all_abbreviation_names = {}

    dir_path = 'protein_entities_standardized/'
    human_sperm_protein_aliases_names_fp = dir_path + 'human_sperm_protein_name_all_abbreviation_aliases.csv'
    if (not os.path.exists(human_sperm_protein_aliases_names_fp)) or (
            os.path.getsize(human_sperm_protein_aliases_names_fp) == 0):
        human_sperm_protein_names_aliases = {}
        human_sperm_protein_names = set()
        human_sperm_protein_symbols_fp = dir_path + 'gene_symbols_mgi_links.csv'
        special_gene_names = {'LOC101928108':'LY6L', 'LOC388849':'CCDC188', 'Tex35':'TEX35'}
        skip_gene_names = ['ZNRD1ASP']
        with codecs.open(human_sperm_protein_symbols_fp, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                if temp[0].strip() in skip_gene_names:
                    continue
                if temp[0].strip() in special_gene_names.keys():
                    human_sperm_protein_names.add(special_gene_names[temp[0].strip()])
                    continue
                human_sperm_protein_names.add(temp[0].strip())
        print('human_sperm_protein_names: ' + str(len(human_sperm_protein_names)))

        uniprot_gene_names = {}
        uniprot_gene_symbols_fp = dir_path + 'uniprot_human_gene_names.tab.gz'
        with gzip.open(uniprot_gene_symbols_fp, 'rb', 'utf-8') as gzipped_file:
            for line in gzipped_file:
                if 'Entry' in line:
                    continue
                temp = line.strip().split('	')
                if len(temp) != 3:
                    continue
                if ';' in temp[2].strip():
                    temp_names = re.sub(';', ' ', temp[2].strip())
                else:
                    temp_names = temp[2].strip()
                names = []
                for name in temp_names.split(' '):
                    if name.strip() != '':
                        names.append(name.strip())
                if len(names) == 1:
                    uniprot_gene_names[names[0]] = names
                else:
                    uniprot_gene_names[names[0]] = names[1:]
        print(len(uniprot_gene_names.keys()))

        ncbi_aliases_gene_names = {}
        ncbi_aliases_gene_symbols_fp = dir_path + 'protein_entities_standardized_ncbi.csv'
        with codecs.open(ncbi_aliases_gene_symbols_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 1, None):
                temp = line.strip().split('$')
                if len(temp) != 5:
                    continue
                names = []
                for name in temp[2].split(','):
                    if name.strip() != '':
                        names.append(name.strip())
                if len(names) == 0:
                    ncbi_aliases_gene_names[temp[0].strip()] = [temp[0].strip()]
                else:
                    ncbi_aliases_gene_names[temp[0].strip()] = names
        print(len(ncbi_aliases_gene_names.keys()))

        genecards_gene_names = {}
        genecards_gene_symbols_fp = dir_path + 'protein_entities_standardized_genecards.csv'
        with codecs.open(genecards_gene_symbols_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 1, None):
                temp = line.strip().split('$')
                if len(temp) < 2:
                    continue
                names = []
                for name in temp[1:]:
                    if name.strip() == '':
                        continue
                    if len(name.strip().split(' ')) == 1:
                        names.append(name.strip())
                if len(names) == 0:
                    genecards_gene_names[temp[0].strip()] = [temp[0].strip()]
                else:
                    genecards_gene_names[temp[0].strip()] = names
        print(len(genecards_gene_names.keys()))

        for protein_name in human_sperm_protein_names:
            if protein_name in uniprot_gene_names.keys():
                human_sperm_protein_names_aliases[protein_name] = uniprot_gene_names[protein_name]
            else:
                for key in uniprot_gene_names.keys():
                    if protein_name in uniprot_gene_names[key]:
                        human_sperm_protein_names_aliases[protein_name] = uniprot_gene_names[key]
                        break
            if protein_name not in human_sperm_protein_names_aliases.keys():
                if protein_name in ncbi_aliases_gene_names.keys():
                    human_sperm_protein_names_aliases[protein_name] = ncbi_aliases_gene_names[protein_name]
                else:
                    for key in ncbi_aliases_gene_names.keys():
                        if protein_name in ncbi_aliases_gene_names[key]:
                            human_sperm_protein_names_aliases[protein_name] = ncbi_aliases_gene_names[key]
                            break
            else:
                if protein_name in ncbi_aliases_gene_names.keys():
                    human_sperm_protein_names_aliases[protein_name] = ncbi_aliases_gene_names[protein_name]+human_sperm_protein_names_aliases[protein_name]

            if protein_name not in human_sperm_protein_names_aliases.keys():
                if protein_name in genecards_gene_names.keys():
                    human_sperm_protein_names_aliases[protein_name] = genecards_gene_names[protein_name]
                else:
                    for key in genecards_gene_names.keys():
                        if protein_name in genecards_gene_names[key]:
                            human_sperm_protein_names_aliases[protein_name] = genecards_gene_names[key]
                            break
            # else:
            #     if protein_name in genecards_gene_names.keys():
            #         human_sperm_protein_names_aliases[protein_name] = genecards_gene_names[protein_name]+human_sperm_protein_names_aliases[protein_name]

            if protein_name not in human_sperm_protein_names_aliases.keys():
                print(protein_name+'---------------------')
        print(len(human_sperm_protein_names_aliases.keys()))

        for key in human_sperm_protein_names_aliases.keys():
            if len(human_sperm_protein_names_aliases[key]) == 1:
                if human_sperm_protein_names_aliases[key][0] == key:
                    for pro in ncbi_aliases_gene_names.keys():
                        if key in ncbi_aliases_gene_names[pro]:
                            human_sperm_protein_names_aliases[key] = ncbi_aliases_gene_names[pro] + [pro]

        with open(human_sperm_protein_aliases_names_fp, 'w') as output_file:
            for protein_name in human_sperm_protein_names:
                if protein_name not in human_sperm_protein_names_aliases.keys():
                    print(protein_name+'+++++++++++++++++++++')
                    continue
                if len(human_sperm_protein_names_aliases[protein_name]) > 1:
                    while protein_name in human_sperm_protein_names_aliases[protein_name]:
                        human_sperm_protein_names_aliases[protein_name].remove(protein_name)
                aliases = human_sperm_protein_names_aliases[protein_name]
                if len(aliases) == 0:
                    aliases = [protein_name]
                human_sperm_protein_all_abbreviation_names[protein_name] = aliases
                output_file.write(protein_name + ',' + ','.join(list(set(aliases)))+'\n')
        print(len(human_sperm_protein_all_abbreviation_names.keys()))
    else:
        with codecs.open(human_sperm_protein_aliases_names_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                human_sperm_protein_all_abbreviation_names[temp[0]] = temp[1:]
        print(len(human_sperm_protein_all_abbreviation_names.keys()))

    return human_sperm_protein_all_abbreviation_names

# Count the number of Gene Entry for each gene in Uniprot database
def prepare_human_sperm_gene_uniprot_entry():

    dir_path = 'protein_entities_standardized/'

    human_sperm_protein_all_abbreviation_names = prepare_human_sperm_protein_name()

    # 处理从uniprot 上下载的人类蛋白质entry 与protein names

    uniprot_primary_name_gene_names = {}
    gene_primary_name_uniprot_entries = {}

    gz_fp = dir_path + 'uniprot_human_entry_protein_names.xlsx.gz'
    uniprot_entry_name_fp = un_gz(gz_fp)
    try:
        data = xlrd.open_workbook(uniprot_entry_name_fp)
    except Exception, e:
        print str(e)

    table = data.sheet_by_name('Sheet0')
    # colnames = table.row_values(0)
    # print(colnames)

    entries = table.col_values(0)
    gene_names = table.col_values(1)
    protein_names = table.col_values(2)

    print(len(entries), len(gene_names), len(protein_names))
    print(entries[1], gene_names[1], protein_names[1])

    for i in range(1, len(entries)):
        if (entries[i].strip() == '') or (gene_names[i].strip() == ''):
            continue
        if ';' in gene_names[i].strip():
            temp_names = re.sub(';', ' ', gene_names[i].strip())
            individual_gene_names = [x.strip() for x in temp_names.split(' ')]
        else:
            individual_gene_names = [x.strip() for x in gene_names[i].strip().split(' ')]

        primary_gene_name = individual_gene_names[0].strip()
        if primary_gene_name not in gene_primary_name_uniprot_entries.keys():
            gene_primary_name_uniprot_entries[primary_gene_name] = set()
            uniprot_primary_name_gene_names[primary_gene_name] = []
        gene_primary_name_uniprot_entries[primary_gene_name].add(entries[i].strip())
        uniprot_primary_name_gene_names[primary_gene_name] += individual_gene_names[1:]

    print(len(gene_primary_name_uniprot_entries.keys()))

    human_sperm_protein_uniprot_entry = {}
    filter_protein_names = set()
    for protein in human_sperm_protein_all_abbreviation_names.keys():
        flag = 0
        if protein in gene_primary_name_uniprot_entries.keys():
            human_sperm_protein_uniprot_entry[protein] = ','.join(list(gene_primary_name_uniprot_entries[protein]))
            continue
        else:
            for other_gene_abb in human_sperm_protein_all_abbreviation_names[protein]:
                if other_gene_abb in gene_primary_name_uniprot_entries.keys():
                    human_sperm_protein_uniprot_entry[protein] = ','.join(list(gene_primary_name_uniprot_entries[other_gene_abb]))
                    flag = 1
                    break
        if flag == 0:
            filter_protein_names.add(protein)

    for f_protein_name in filter_protein_names:
        for protein_name in uniprot_primary_name_gene_names.keys():
            if f_protein_name in uniprot_primary_name_gene_names[protein_name]:
                human_sperm_protein_uniprot_entry[f_protein_name] = ','.join(list(gene_primary_name_uniprot_entries[protein_name]))
                break

    print(len(human_sperm_protein_uniprot_entry.keys()))

    out_file_path = dir_path + 'human_sperm_protein_uniprot_entry_gene_names.csv'
    if not os.path.exists(out_file_path):
        with open(out_file_path, 'w') as output_file:
            for key in human_sperm_protein_uniprot_entry.keys():
                output_file.write(key.strip() + ',' + human_sperm_protein_uniprot_entry[key] + '\n')

# Prepared human protein abbreviation collections and their abbreviations, from two databases, unprot,ncbi
def prepare_human_protein_names():

    raw_data_fp = 'protein_entities_standardized/'
    out_file_path = raw_data_fp + 'human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    protein_official_symbol_names = {}
    if not os.path.exists(out_file_path):
        ncbi_gene_official_symbols_fp = raw_data_fp + 'protein_entities_standardized_ncbi.csv'
        ncbi_gene_info_fp = 'text_corpus/ncbi/gene_info'
        uniprot_gene_official_symbols_fp = raw_data_fp + 'uniprot_human_gene_names.tab.gz'
        ncbi_human_gene_types = {}
        with codecs.open(ncbi_gene_info_fp, "rb", "utf-8") as input_file:
            for line in input_file:
                if '#' == line.strip()[0]:
                    continue
                temp = line.strip().split('	')
                if len(temp) < 10:
                    continue
                if temp[0] != '9606':
                    continue
                ncbi_human_gene_types[temp[2].strip()] = temp[9].strip()
        print('ncbi_human_gene_types: '+str(len(ncbi_human_gene_types.keys())))

        ncbi_gene_official_symbols = {}
        uniprot_gene_official_symbols = {}
        with codecs.open(ncbi_gene_official_symbols_fp, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 1, None):
                temp = line.strip().split('$')
                if temp[2].strip() == '':
                    ncbi_gene_official_symbols[temp[0].strip()] = [temp[0].strip()]
                else:
                    names = [temp[0].strip()]+temp[2].strip().split(',')
                    ncbi_gene_official_symbols[temp[0].strip()] = [x.strip() for x in names]
        print(len(ncbi_gene_official_symbols.keys()))

        with gzip.open(uniprot_gene_official_symbols_fp, 'rb', 'utf-8') as gzipped_file:
            middle_variate_gene_symbols = set()
            for line in gzipped_file:
                if 'Entry' in line:
                    continue
                temp = line.strip().split('	')
                if len(temp) != 3:
                    continue
                if temp[2].strip() == '':
                    continue
                # if ';' in temp[2].strip():
                #     temp_names = re.sub(';', ' ', temp[2].strip())
                # else:
                #     temp_names = temp[2].strip()
                names = []
                gene_symbol_str = re.sub(' ', ';', temp[2].strip())
                # gene_symbol_str = eval(repr(gene_symbol_str).replace('/', ';'))
                gene_names = gene_symbol_str.split(';')
                for name in gene_names:
                    if name.strip() != '':
                        names.append(name.strip())
                if len(names) == 0:
                    continue
                if names[0] not in middle_variate_gene_symbols:
                    uniprot_gene_official_symbols[names[0]] = names
                    middle_variate_gene_symbols.add(names[0])
                else:
                    uniprot_gene_official_symbols[names[0]] += names
        print(len(uniprot_gene_official_symbols.keys()))

        # uniprot_gene_official_symbols_difference = set(uniprot_gene_official_symbols.keys()).difference(
        #     ncbi_gene_official_symbols.keys())
        #
        # ncbi_gene_official_symbols_difference = set(ncbi_gene_official_symbols.keys()).difference(
        #     uniprot_gene_official_symbols.keys())
        #
        # gene_official_symbols_intersection = set(uniprot_gene_official_symbols.keys()).intersection(
        #     ncbi_gene_official_symbols.keys())
        #
        # print('gene_official_symbols_intersection: '+str(len(gene_official_symbols_intersection)))
        #
        # for gene_official_symbol in gene_official_symbols_intersection:
        #     gene_official_symbol_names[gene_official_symbol] = set(ncbi_gene_official_symbols[gene_official_symbol] + uniprot_gene_official_symbols[gene_official_symbol])
        #
        # for gene_official_symbol in uniprot_gene_official_symbols_difference:
        #     flag = 0
        #     for gene_name in gene_official_symbols_intersection:
        #         if gene_official_symbol in gene_official_symbol_names[gene_name]:
        #             flag = 1
        #             break
        #     if flag == 0:
        #         gene_official_symbol_names[gene_official_symbol] = set(uniprot_gene_official_symbols[gene_official_symbol])
        #
        # for gene_official_symbol in ncbi_gene_official_symbols_difference:
        #     flag = 0
        #     for gene_name in gene_official_symbols_intersection:
        #         if gene_official_symbol in gene_official_symbol_names[gene_name]:
        #             flag = 1
        #             break
        #     if flag == 0:
        #         gene_official_symbol_names[gene_official_symbol] = set(ncbi_gene_official_symbols[gene_official_symbol])

        gene_official_symbols_union = set(uniprot_gene_official_symbols.keys()).union(
            ncbi_gene_official_symbols.keys())
        for gene_symbol in gene_official_symbols_union:
            if gene_symbol in uniprot_gene_official_symbols.keys():
                protein_official_symbol_names[gene_symbol] = uniprot_gene_official_symbols[gene_symbol]
                if gene_symbol in ncbi_gene_official_symbols.keys():
                    protein_official_symbol_names[gene_symbol] += ncbi_gene_official_symbols[gene_symbol]
            else:
                if gene_symbol in ncbi_human_gene_types.keys():
                    if ncbi_human_gene_types[gene_symbol] == 'protein-coding':
                        protein_official_symbol_names[gene_symbol] = ncbi_gene_official_symbols[gene_symbol]
        with open(out_file_path, 'w') as output_file:
            for gene_name in protein_official_symbol_names.keys():
                clean_gene_names = list(set(protein_official_symbol_names[gene_name]))
                output_file.write(gene_name + ',' + ','.join(clean_gene_names) + '\n')
    else:
        with codecs.open(out_file_path, 'rb', 'utf-8') as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                protein_official_symbol_names[temp[0]] = set(temp[1:])
    print(len(protein_official_symbol_names.keys()))

    return protein_official_symbol_names


if __name__ == '__main__':
    prepare_human_sperm_protein_name()
    # prepare_human_sperm_gene_uniprot_entry()
