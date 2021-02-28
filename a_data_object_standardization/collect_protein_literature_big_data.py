#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 21/10/2018 6:06 PM
# @Author  : ggguo
# @Site    : SYM PROJECT
# @File    : collect_protein_literature_big_data.py
# @Software: SYM application case

###############
###Intermediate process code, for user reference only
###############

from os import listdir
import os
import codecs
from itertools import islice
import re


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

# According to PMID get the title and abstract of the literature
def extract_protein_pmid_title_abstract(dir_path):

    extract_text_data_from_ncbi_dir_path = 'ncbi/parse_literatures_human/'
    extract_text_data_from_uniprot_dir_path = 'uniprot/parse_literatures_human/'
    file_type = '.txt'
    extract_text_data_dir_path = [extract_text_data_from_ncbi_dir_path,extract_text_data_from_uniprot_dir_path]
    extract_text_data_by_source = {}
    gene_ids_set = set()
    for etd_dir_path in extract_text_data_dir_path:
        file_names = [f for f in listdir(etd_dir_path) if f.endswith(file_type)]
        print(len(file_names))
        for file_name in file_names:
            key = re.sub('\.txt', '', file_name).strip().split('_')[-1].strip()
            if key in gene_ids_set:
                print(key)
                continue
            temp_country_journal = []
            with codecs.open(etd_dir_path + file_name, "rb", "utf-8") as input_file:
                for line in islice(input_file.readlines(), 0, None):
                    temp = line.strip().split('$')
                    if len(temp) < 4:
                        continue
                    if not temp[0].strip().isdigit():
                        continue
                    temp_country_journal.append([temp[0]]+temp[2:])
            gene_ids_set.add(key)
            extract_text_data_by_source[key] = temp_country_journal
        print(len(extract_text_data_by_source.keys()))

    dir_path_local = 'data_help/'
    id_protein_entry_ncbi, id_protein_entry_uniprot = id_protein_entry_match(dir_path_local)

    extract_text_data_by_protein = {}
    pmids_total_set = set()
    for id_gene in extract_text_data_by_source.keys():
        if id_gene in id_protein_entry_ncbi.keys():
            if id_protein_entry_ncbi[id_gene] in extract_text_data_by_protein.keys():
                for protein_name in extract_text_data_by_source[id_gene]:
                    if protein_name[0] in extract_text_data_by_protein[id_protein_entry_ncbi[id_gene]].keys():
                        continue
                    else:
                        pmids_total_set.add(protein_name[0])
                        extract_text_data_by_protein[id_protein_entry_ncbi[id_gene]][protein_name[0]] = '$'.join(protein_name[1:])
            else:
                extract_text_data_by_protein[id_protein_entry_ncbi[id_gene]] = {}
                for protein_name in extract_text_data_by_source[id_gene]:
                    pmids_total_set.add(protein_name[0])
                    extract_text_data_by_protein[id_protein_entry_ncbi[id_gene]][protein_name[0]] = '$'.join(protein_name[1:])
        elif id_gene in id_protein_entry_uniprot.keys():
            if id_protein_entry_uniprot[id_gene] in extract_text_data_by_protein.keys():
                for protein_name in extract_text_data_by_source[id_gene]:
                    if protein_name[0] in extract_text_data_by_protein[id_protein_entry_uniprot[id_gene]].keys():
                        continue
                    else:
                        pmids_total_set.add(protein_name[0])
                        extract_text_data_by_protein[id_protein_entry_uniprot[id_gene]][protein_name[0]] = '$'.join(protein_name[1:])
            else:
                extract_text_data_by_protein[id_protein_entry_uniprot[id_gene]] = {}
                for protein_name in extract_text_data_by_source[id_gene]:
                    pmids_total_set.add(protein_name[0])
                    extract_text_data_by_protein[id_protein_entry_uniprot[id_gene]][protein_name[0]] = '$'.join(protein_name[1:])
        else:
            print(id_gene)
            continue
    print(len(extract_text_data_by_protein.keys()),len(pmids_total_set))

    out_file_path_pta = dir_path + 'protein_entry_pmid_title_abstract.csv'
    if not os.path.exists(out_file_path_pta):
        with open(out_file_path_pta, 'a') as output_file:
            for protein_entry in extract_text_data_by_protein.keys():
                for pmid in extract_text_data_by_protein[protein_entry]:
                    output_file.write(protein_entry + '$' + pmid + '$' + extract_text_data_by_protein[protein_entry][pmid] + '\n')

    pmid_texts = {}
    protein_pmids = {}
    pmids = set()
    for protein_name in extract_text_data_by_protein:
        temp_pmids = []
        for pmid in extract_text_data_by_protein[protein_name]:
            if pmid in pmids:
                continue
            if extract_text_data_by_protein[protein_name][pmid].strip() == '$':
                continue
            temp_pmids.append(pmid)
            title_abstract = extract_text_data_by_protein[protein_name][pmid].strip()
            if '?' == title_abstract[-1]:
                pmid_texts[pmid] = title_abstract[:-1]
            else:
                pmid_texts[pmid] = title_abstract
            pmids.add(pmid)
        if len(temp_pmids) == 0:
            continue
        protein_pmids[protein_name] = temp_pmids

    return pmid_texts,protein_pmids

# Prepare big data for protein literature
def prepare_protein_literature_big_data(dir_path):

    protein_entities_dictionary = {}
    file_path_protein_dictionary = dir_path + 'protein_entity_standardized_dictionary.csv'
    with codecs.open(file_path_protein_dictionary, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            protein_entities_dictionary[temp[0].strip()] = temp[2:]
    print('The amount of protein_entities_dictionary: '+ str(len(protein_entities_dictionary)))

    pmid_texts, protein_pmids = extract_protein_pmid_title_abstract(dir_path)

    for protein_name in protein_pmids:
        sorted_protein_dictionary = sorted(protein_entities_dictionary[protein_name]+[protein_name],key = lambda i:len(i),reverse=True)
        sorted_protein_dictionary_filter = []
        for strr in sorted_protein_dictionary:
            if len(strr) >= 2: #len(protein_name)
                sorted_protein_dictionary_filter.append(strr)

        count = 0
        for pmid in protein_pmids[protein_name]:
            raw_text = re.sub('\$',' ',pmid_texts[pmid])
            if raw_text.strip() == '':
                continue
            for alternative_name in sorted_protein_dictionary_filter:
                if ' '+alternative_name+' ' in raw_text:
                    raw_text = raw_text.replace(' '+alternative_name+' ', ' ' + protein_name + ' ')
                    count += 1
                elif ' '+alternative_name+',' in raw_text:
                    raw_text = raw_text.replace(' '+alternative_name+',', ' ' + protein_name + ',')
                    count += 1
                elif ' ('+alternative_name+')' in raw_text:
                    raw_text = raw_text.replace(' ('+alternative_name+')', ' ' + protein_name)
                    count += 1
                elif ' [' + alternative_name + ']' in raw_text:
                    raw_text = raw_text.replace(' [' + alternative_name + ']', ' ' + protein_name)
                    count += 1
                elif ' '+alternative_name+'.' in raw_text:
                    raw_text = raw_text.replace(' '+alternative_name+'.',' '+protein_name+'.')
                    count += 1
                else:
                    continue
            pmid_texts[pmid] = raw_text
        print('The amount of the text is replaced! ' + str(count) + ' ' + protein_name)

    out_file_path_standardized = dir_path + 'protein_entry_pmid_title_abstract_standardized.csv'
    if not os.path.exists(out_file_path_standardized):
        with open(out_file_path_standardized, 'a') as output_file:
            for pmid in pmid_texts.keys():
                if pmid_texts[pmid] == '':
                    continue
                output_file.write(pmid + '$' + pmid_texts[pmid] + '\n')


if __name__ == '__main__':

    dir_path = 'protein_entities_standardized/'
    prepare_protein_literature_big_data(dir_path)