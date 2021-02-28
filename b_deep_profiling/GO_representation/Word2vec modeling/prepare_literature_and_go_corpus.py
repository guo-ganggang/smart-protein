#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 14/11/2018 3:49 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : prepare_literature_and_go_corpus.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import os,re
import gzip
import shutil
import time
import ebi_api_fetch
import json
import codecs
from os import listdir
import platform
from itertools import islice
import data_helper as dh

import sys

reload(sys)
sys.setdefaultencoding('utf-8')

species_dict = {'human': '9606', 'zebrafish': '7955', 'fruit_fly': '7227','rat': '10116', 'mouse': '10090'}

# UnZip File
def un_gz(file_name):
    """ungz zip file"""
    f_name = file_name.replace(".gz", "")
    with gzip.open(file_name, 'rb') as f_in:
        with open(f_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return f_name

# Collect GO definition information and
# comment text original information to prepare corpus
def fetch_protein_go_from_ebi(file_path,path_type):

    gz_go_corpus_fname = file_path + 'raw_data'+path_type+'goa_uniprot_all.gpa.gz'
    go_corpus_fname = file_path + 'raw_data'+path_type+'goa_uniprot_all.gpa'
    if not os.path.exists(go_corpus_fname):
        un_gz(gz_go_corpus_fname)

    file_name_feature = 'go_id_'
    completed_go_ids = set()
    out_file_dir = file_path + 'raw_data' + path_type
    file_names = [f for f in listdir(out_file_dir) if f.startswith(file_name_feature)]
    for file_name in file_names:
        with codecs.open(out_file_dir+file_name, 'rb') as input_file:
            for line in input_file:
                temp = line.strip().split(',')
                completed_go_ids.add(temp[0].strip())
    print(len(completed_go_ids))

    out_file_path_c = file_path + 'raw_data'+path_type+'go_id_name_definition_comment.csv'
    out_file_path_f = file_path + 'raw_data'+path_type+'go_id_for_failed_collect.csv'
    go_ids = set()
    with codecs.open(go_corpus_fname, 'rb') as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) < 4:
                continue
            go_id = temp[3].strip()
            if go_id in go_ids:
                continue
            # uniprot_id = temp[1].strip()
            print(go_id) #GO%3A0016765
            requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"+ go_id \
                         +"/ancestors"
            html = None
            for k in range(3):
                if html is None:
                    html = ebi_api_fetch.crawler_html_by_url(requestURL)
                else:
                    break
                if k > 0:
                    time.sleep(50)
            if html is None:
                with open(out_file_path_f, 'a') as output_file:
                    output_file.write(go_id + '\n')
                continue
            else:
                html_dict = json.loads(html)

                # print(html_dict.keys())
                # print(html_dict['results'])
                # print(uniprot_id)
                # results = html_dict['results'][0]
                # for result in results.keys():
                #     print(result)
                # break

                with open(out_file_path_c, 'a') as output_file:
                    output_file.write(go_id+','+str(html_dict['results']) + '\n')
                    go_ids.add(go_id)

# Extract GO definition information and comment text to prepare corpus
def parse_protein_go_definition_comment_corpus(file_path):
    out_file_path_c = file_path + 'raw_data' + path_type + 'go_id_name_definition_comment.csv'
    go_corpus_path = file_path + 'text_corpus' + path_type + 'corpus_go_id_name_definition.csv'
    with open(go_corpus_path, 'w') as output_file:
        with codecs.open(out_file_path_c, 'rb') as input_file:
            for line in input_file:
                temp = line.strip().split(',')
                dc_str = ','.join([str(x)for x in temp[1:]]).strip()
                # print(dc_str[1:-1])
                dc_json = eval(dc_str[1:-1])
                # print(dc_json['definition']['text'])

                if 'definition' in dc_json.keys():
                    name_synonyms = dc_json['name'].strip()
                    # definition_comment = dc_json['definition']['text'].strip()
                    # if 'comment' in dc_json.keys():
                    #     comment_text = re.sub('Note that ','',dc_json['comment'].strip())
                    #     definition_comment = definition_comment + dc_json['comment'].strip()
                    if 'synonyms' in dc_json.keys():
                        # print(temp[0].strip())
                        synonymses = []
                        for synonyms in dc_json['synonyms']:
                            synonymses.append(synonyms['name'].strip())
                        name_synonyms = name_synonyms + ';' + ';'.join(synonymses)
                    output_file.write(temp[0].strip() + '$' + name_synonyms + '$' + dc_json['definition']['text'].strip() + '\n')

# Literature on each human protein, as well as literature
# on human homologous mouse phenotypes, literature
# on human genes and diseases, literature on PFAM family, and literature on GO
def prepare_literature_corpus(dir_path):

    subdir_path_fname = 'data_help/text_corpus/extract_literatures/pubmed_id_date_title_abstract_journal.csv'
    pmid_title_abstract = {}
    pmids = set()
    with codecs.open(subdir_path_fname, 'rb', 'utf-8') as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if len(temp) < 6:
                continue
            title_abstrct = ''.join(temp[2:4])
            if (temp[0].strip() in pmids) or (title_abstrct == '??'):
                continue
            pmid_title_abstract[temp[0].strip()] = title_abstrct
            pmids.add(temp[0].strip())
    print('pmid_title_abstract: ' + str(len(pmid_title_abstract.keys())))

    write_out_fp = dir_path + 'text_corpus/corpus_literatures_title_abstract_by_pmid.csv'
    with open(write_out_fp, 'w') as output_file:
        for pmid in pmid_title_abstract.keys():
            output_file.write(pmid.strip() + '$' + pmid_title_abstract[pmid.strip()] + '\n')

# Training sets of homologous genes of human genes and other 4 species that can be combined
def prepare_protein_go_ids_for_train(dir_path):

    human_protein_symbol_ncbi_ids, human_protein_symbol_uniprot_entries, train_protein_symbols_aliases, _= dh.prepare_human_protein_symbols_ids()

    train_protein_symbol_ncbi_ids = {}
    train_ncbi_ids = set()
    train_protein_symbol_uniprot_entries = {}
    train_uniprot_entries= set()
    for protein_symbol in train_protein_symbols_aliases.keys():
        if protein_symbol in human_protein_symbol_ncbi_ids.keys():
            train_protein_symbol_ncbi_ids[protein_symbol] = human_protein_symbol_ncbi_ids[protein_symbol].split(';')
            train_ncbi_ids = train_ncbi_ids | set(human_protein_symbol_ncbi_ids[protein_symbol].split(';'))
        else:
            for protein in train_protein_symbols_aliases[protein_symbol]:
                if protein.strip() in train_protein_symbol_ncbi_ids.keys():
                    train_protein_symbol_ncbi_ids[protein_symbol] = human_protein_symbol_ncbi_ids[protein.strip()].split(';')
                    train_ncbi_ids = train_ncbi_ids | set(human_protein_symbol_ncbi_ids[protein.strip()].split(';'))
                    break
        if protein_symbol in human_protein_symbol_uniprot_entries.keys():
            train_protein_symbol_uniprot_entries[protein_symbol] = human_protein_symbol_uniprot_entries[protein_symbol].split(';')
            train_uniprot_entries = train_uniprot_entries | set(human_protein_symbol_uniprot_entries[protein_symbol].split(';'))
        else:
            for protein in train_protein_symbols_aliases[protein_symbol]:
                if protein.strip() in human_protein_symbol_uniprot_entries.keys():
                    train_protein_symbol_uniprot_entries[protein_symbol] = human_protein_symbol_uniprot_entries[protein.strip()].split(';')
                    train_uniprot_entries = train_uniprot_entries | set(human_protein_symbol_uniprot_entries[protein.strip()].split(';'))
                    break
    print('train_protein_symbol_ncbi_ids: ' + str(len(train_protein_symbol_ncbi_ids.keys())))
    print('train_protein_symbol_uniprot_entries: ' + str(len(train_protein_symbol_uniprot_entries.keys())))

    mapped_train_protein_symbols = set(train_protein_symbol_ncbi_ids.keys()).union(set(train_protein_symbol_uniprot_entries.keys()))
    print('mapped_train_protein_symbols: ' + str(len(mapped_train_protein_symbols)))

    species_ncbi_id_go_ids_train = {}
    human_ncbi_id_go_ids = dh.get_gene_pubmed_list_from_ncbi(dir_path,'human')
    species_ncbi_ids = set(human_ncbi_id_go_ids.keys())
    for ncbi_id in train_ncbi_ids:
        if ncbi_id in species_ncbi_ids:
            species_ncbi_id_go_ids_train[ncbi_id] = human_ncbi_id_go_ids[ncbi_id]
        else:
            species_ncbi_id_go_ids_train[ncbi_id] = []
    print('species_ncbi_id_go_ids_train: ' + str(len(species_ncbi_id_go_ids_train.keys())))

    species_uniprot_entry_go_ids_train = {}
    human_uniprot_entry_go_ids = dh.get_gene_pubmed_list_from_uniprot(dir_path,'human')
    species_uniprot_entries = set(human_uniprot_entry_go_ids.keys())
    for uniprot_entry in train_uniprot_entries:
        if uniprot_entry in species_uniprot_entries:
            species_uniprot_entry_go_ids_train[uniprot_entry] = human_uniprot_entry_go_ids[uniprot_entry]
        else:
            species_uniprot_entry_go_ids_train[uniprot_entry] = []
    print('species_uniprot_entry_go_ids_train: ' + str(len(species_uniprot_entry_go_ids_train.keys())))

    for key in species_dict.keys():
        if key == 'human':
            continue
        species_taxid = species_dict[key]
        homologous_gene_id_human_gene_id, homologous_gene_entry_human_gene_entry = dh.other_species_homologous_gene_names_ids(species_taxid)
        homologous_gene_human_gene_ids = set(homologous_gene_id_human_gene_id.values())
        homologous_gene_human_gene_entries = set(homologous_gene_entry_human_gene_entry.values())

        temp_human_gene_id_spcies_gene_ids = {}
        for animal_gene_id in homologous_gene_id_human_gene_id.keys():
            human_gene_id = homologous_gene_id_human_gene_id[animal_gene_id]
            if human_gene_id not in temp_human_gene_id_spcies_gene_ids.keys():
                temp_human_gene_id_spcies_gene_ids[human_gene_id] = set()
            temp_human_gene_id_spcies_gene_ids[human_gene_id].add(animal_gene_id)

        temp_human_gene_entry_spcies_gene_ids = {}
        for animal_gene_entry in homologous_gene_entry_human_gene_entry.keys():
            human_gene_entry = homologous_gene_entry_human_gene_entry[animal_gene_entry]
            if human_gene_entry not in temp_human_gene_entry_spcies_gene_ids.keys():
                temp_human_gene_entry_spcies_gene_ids[human_gene_entry] = set()
            temp_human_gene_entry_spcies_gene_ids[human_gene_entry].add(animal_gene_entry)

        species_ncbi_id_go_ids = dh.get_gene_pubmed_list_from_ncbi(dir_path,key)
        species_uniprot_entry_go_ids = dh.get_gene_pubmed_list_from_uniprot(dir_path,key)

        for ncbi_id in species_ncbi_id_go_ids_train.keys():
            if ncbi_id in homologous_gene_human_gene_ids:
                spcies_ncbi_ids = temp_human_gene_id_spcies_gene_ids[ncbi_id]
                for spcies_ncbi_id in spcies_ncbi_ids:
                    if spcies_ncbi_id in species_ncbi_id_go_ids.keys():
                        species_ncbi_id_go_ids_train[ncbi_id] += species_ncbi_id_go_ids[spcies_ncbi_id]
        print(key + ' species_ncbi_id_go_ids_train: ' + str(len(species_ncbi_id_go_ids_train.keys())))

        for uniprot_entry in species_uniprot_entry_go_ids_train.keys():
            if uniprot_entry in homologous_gene_human_gene_entries:
                spcies_uniprot_entries = temp_human_gene_entry_spcies_gene_ids[uniprot_entry]
                for spcies_uniprot_entry in spcies_uniprot_entries:
                    if spcies_uniprot_entry in species_uniprot_entry_go_ids.keys():
                        species_uniprot_entry_go_ids_train[uniprot_entry] += species_uniprot_entry_go_ids[spcies_uniprot_entry]
        print(key + ' species_uniprot_entry_go_ids_train: ' + str(len(species_uniprot_entry_go_ids_train.keys())))

    mapped_train_ncbi_ids = train_ncbi_ids.intersection(set(species_ncbi_id_go_ids_train.keys()))
    mapped_train_uniprot_entries = train_uniprot_entries.intersection(set(species_uniprot_entry_go_ids_train.keys()))
    print('mapped_train_ncbi_ids: ' + str(len(mapped_train_ncbi_ids)))
    print('mapped_train_uniprot_entries: ' + str(len(mapped_train_uniprot_entries)))

    train_protein_symbol_go_ids = {}
    temp_train_protein_symbols = set()
    for train_protein_symbol in mapped_train_protein_symbols:
        flag = 0
        if train_protein_symbol in train_protein_symbol_ncbi_ids.keys():
            for ncbi_id in train_protein_symbol_ncbi_ids[train_protein_symbol]:
                if ncbi_id.strip() == '':
                    continue
                if ncbi_id.strip() in mapped_train_ncbi_ids:
                    if train_protein_symbol not in temp_train_protein_symbols:
                        train_protein_symbol_go_ids[train_protein_symbol] = []
                        temp_train_protein_symbols.add(train_protein_symbol)
                    train_protein_symbol_go_ids[train_protein_symbol] += species_ncbi_id_go_ids_train[ncbi_id.strip()]
                    flag = 1
            if train_protein_symbol in train_protein_symbol_uniprot_entries.keys():
                for uniprot_entry in train_protein_symbol_uniprot_entries[train_protein_symbol]:
                    if uniprot_entry.strip() == '':
                        continue
                    if uniprot_entry.strip() in mapped_train_uniprot_entries:
                        if train_protein_symbol not in temp_train_protein_symbols:
                            train_protein_symbol_go_ids[train_protein_symbol] = []
                            temp_train_protein_symbols.add(train_protein_symbol)
                        train_protein_symbol_go_ids[train_protein_symbol] += species_uniprot_entry_go_ids_train[uniprot_entry.strip()]
                        flag = 1
        else:
            for uniprot_entry in train_protein_symbol_uniprot_entries[train_protein_symbol]:
                if uniprot_entry.strip() == '':
                    continue
                if uniprot_entry.strip() in mapped_train_uniprot_entries:
                    if train_protein_symbol not in temp_train_protein_symbols:
                        train_protein_symbol_go_ids[train_protein_symbol] = []
                        temp_train_protein_symbols.add(train_protein_symbol)
                    train_protein_symbol_go_ids[train_protein_symbol] += species_uniprot_entry_go_ids_train[
                        uniprot_entry.strip()]
                    flag = 1
        if flag == 0:
            print(train_protein_symbol)
    print('train_protein_symbol_go_ids: ' + str(len(train_protein_symbol_go_ids.keys())))

    write_out_fp = dir_path + 'text_corpus/train_protein_symbol_go_ids.csv'
    with open(write_out_fp, 'w') as output_file:
        for train_protein_symbol in train_protein_symbol_go_ids.keys():
            go_ids = list(set(train_protein_symbol_go_ids[train_protein_symbol]))
            if len(go_ids) == 0:
                continue
            output_file.write(train_protein_symbol.strip() + '$' + ';'.join(go_ids) + '\n')

# Prepare the GO TERM text message
def prepare_go_term_text(dir_path):

    go_term_ids_fname = dir_path + 'text_corpus/train_protein_symbol_go_ids.csv'
    train_go_ids = set()
    with codecs.open(go_term_ids_fname, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('$')
            if temp[1].strip() == '':
                continue
            go_ids = [x.strip() for x in temp[1].strip().split(';')]
            train_go_ids = train_go_ids | set(go_ids)
    print('train_go_ids: '+ str(len(train_go_ids)))

    go_id_term = {}
    ncbi_go_term_fname = dir_path + 'raw_data/gene2go'
    with codecs.open(ncbi_go_term_fname, "rb", "utf-8") as input_file:
        for line in input_file:
            if line.strip()[0] == '#':
                continue
            temp = line.strip().split('	')
            if len(temp) < 6:
                continue
            if (temp[2].strip() == '-') or (temp[5].strip() == '-'):
                continue
            go_id_term[temp[2].strip()] = temp[5].strip()
    print(len(go_id_term.keys()))

    interpro_go_term_fname = dir_path + 'raw_data/interpro2go.csv'
    with codecs.open(interpro_go_term_fname, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(' ; ')
            if len(temp) < 2:
                continue
            if '> GO:' not in temp[0].strip():
                continue
            if temp[1].strip() == '':
                continue
            go_term = temp[0].strip().split('> GO:')[1].strip()
            go_id_term[temp[1].strip()] = go_term
    print(len(go_id_term.keys()))

    total_go_id = set(go_id_term.keys())
    train_go_id_term = {}
    for go_id in train_go_ids:
        if go_id not in total_go_id:
            continue
        train_go_id_term[go_id] = go_id_term[go_id]
    print(len(train_go_id_term))

    print(len(train_go_ids.difference(set(train_go_id_term.keys()))))

    write_out_fp = dir_path + 'text_corpus/train_go_id_term.csv'
    with open(write_out_fp, 'w') as output_file:
        for go_id in train_go_id_term.keys():
            output_file.write(go_id.strip() + '$' + train_go_id_term[go_id.strip()] + '\n')

# Merger corpus
def merge_go_corpus(dir_path):
    corpus_fpath = dir_path + 'text_corpus/'
    file_names = [f for f in listdir(corpus_fpath) if f.endswith('csv')]
    write_out_fp = dir_path + 'corpus_literature_definition_for_go.csv'
    with open(write_out_fp, 'w') as output_file:
        for file_name in file_names:
            with codecs.open(corpus_fpath+file_name, "rb", "utf-8") as input_file:
                for line in input_file:
                    temp = line.strip().split('$')
                    if temp[-1].strip() == '':
                        continue
                    output_file.write(temp[-1].strip() + '\n')




if __name__ == '__main__':

    print(platform.system())
    file_path = ''
    path_type = '/'
    if platform.system() == 'Windows':
        file_path = 'vector_protein_go/'
        path_type = '/'

    # fetch_protein_go_from_ebi(file_path,path_type)
    # parse_protein_go_definition_comment_corpus(file_path)
    # prepare_literature_corpus(file_path)
    # prepare_protein_go_ids_for_train(file_path)
    prepare_go_term_text(file_path)
    # merge_go_corpus(file_path)