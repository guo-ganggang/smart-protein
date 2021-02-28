#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 19/12/2018 7:53 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : collect_ppi_score_network.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############


from itertools import islice
import os
import codecs
import pandas as pd
import numpy as np
import re
import data_helper as hhgtd
import gzip
import shutil
import xlrd
from os import listdir

species_dict = {'human': '9606', 'zebrafish': '7955', 'fruit_fly': '7227','rat': '10116', 'mouse': '10090'} #'test': '394'


'''
    1. Read in the PPI data and prepare node2vec network
'''

# UnZip File
def un_gz(file_name):
    """ungz zip file"""
    f_name = file_name.replace(".gz", "")
    with gzip.open(file_name, 'rb') as f_in:
        with open(f_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return f_name

# Read pairs of related protein relationships from large files
def prepare_ppi_network_for_node2vec(dir_path):

    out_fp = dir_path + '/species_protein_pair_score_update/protein_interaction_pair_combined_score'
    raw_data_fp = dir_path + '/raw_data/protein.links.detailed.v10.5.txt' #protein.links.v10.5.txt
    size = 5000000
    data = pd.read_table(raw_data_fp, sep = ' ',chunksize=size, iterator=True,dtype=str)
    for chunk in data:
        count = 0
        for species_taxid in species_dict.values():
            is_species_taxids = chunk['protein1'].astype(str).str.startswith(species_taxid,na=False)
            chunk['is_species_taxid'] = is_species_taxids
            if len(chunk[chunk["is_species_taxid"] == False]) == size:
                continue
            select_chunk = chunk[chunk["is_species_taxid"] == True]
            count += len(select_chunk['is_species_taxid'])
            del select_chunk['is_species_taxid']
            select_chunk.to_csv(out_fp, mode='a', header=False)
            if count == size:
                break
            # print(select_chunk)
            # temp_protein1 = select_chunk['protein1']
            # temp_protein2 = select_chunk['protein2']
            # temp_combined_score = select_chunk['combined_score']
            # print(temp_protein1)
            # print(temp_protein2)
            # print(temp_combined_score)
            # species_protein_stringid_pair_score = pd.concat([select_chunk,species_protein_stringid_pair_score],ignore_index=True)

    # 将大文件按照物种划分成小文件
    for species_name in species_dict.keys():
        write_out_fp = dir_path + '/species_protein_pair_score_update/protein_interaction_pair_combined_score_'+species_name+'.csv'
        if os.path.exists(write_out_fp):
            continue
        with open(write_out_fp, 'w') as output_file:
            with codecs.open(out_fp, "rb", "utf-8") as input_file:
                for line in input_file:
                    if line.strip()[0] == ',':
                        continue
                    temp = line.strip().split(',')
                    taxid = temp[1].strip().split('.')[0]
                    if taxid != species_dict[species_name]:
                        continue
                    output_file.write(line.strip() + '\n')

# Add the updated data to the extracted small file above.
# Since it is a duplicate file, there is no actual update
def suppliment_ppi_network_for_node2vec(dir_path):
    subdir_path = dir_path + '/species_protein_pair_score_update/'
    dir_file_names = [f for f in listdir(subdir_path) if f.endswith('txt')]
    for file_name in dir_file_names:
        species_taxid = file_name.strip().split('.')[0].strip()
        species_name = ''
        for name in species_dict.keys():
            if species_dict[name] == species_taxid:
                species_name = name
                break
        write_out_fp = dir_path + '/species_protein_pair_score/protein_interaction_pair_combined_score_' + species_name + '.csv'
        suppliment_data_fname = subdir_path + file_name
        size = 50000
        data = pd.read_table(suppliment_data_fname, sep=' ', chunksize=size, iterator=True, dtype=str)
        for chunk in data:
            is_species_taxids = chunk['protein1'].astype(str).str.startswith(species_taxid, na=False)
            chunk['is_species_taxid'] = is_species_taxids
            if len(chunk[chunk["is_species_taxid"] == False]) == size:
                continue
            select_chunk = chunk[chunk["is_species_taxid"] == True]
            del select_chunk['is_species_taxid']
            select_chunk.to_csv(write_out_fp, mode='a', header=False)

# Read in the protein relationship pair data of the selected species
def read_species_protein_stringid_pair_score(subdir_path,species_name,protein_ensembl_protein_identifiers):

    dir_file_names = [f for f in listdir(subdir_path) if f.endswith('csv')]
    protein_interaction_pair_score = {}
    for file_name in dir_file_names:
        if species_name not in file_name:
            continue
        with codecs.open(subdir_path+file_name, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                if line.strip()[0] == ',':
                    continue
                temp = line.strip().split(',')
                protein1 = temp[1].strip().split('.')[1].strip()
                protein2 = temp[2].strip().split('.')[1].strip()
                if (protein1 not in protein_ensembl_protein_identifiers) and (protein2 not in protein_ensembl_protein_identifiers):
                    continue
                protein_interaction_pair_score[protein1+','+protein2] = temp[3].strip()
    print(species_name + ' number: protein_interaction_pair_score: ' + str(
        len(protein_interaction_pair_score.keys())))

    return protein_interaction_pair_score

# Matches the relationship pairs of a given set of proteins
def mapping_required_human_protein_relationship_pair(dir_path):

    human_protein_official_symbol_ncbi_ids, human_protein_official_symbol_uniprot_entries, train_protein_official_symbols_aliases,_ = hhgtd.prepare_human_protein_symbols_ids()

    human_protein_ncbi_ids = set()
    ncbi_id_human_protein = {}
    for protein_official_symbol in human_protein_official_symbol_ncbi_ids.keys():
        ncbi_ids = human_protein_official_symbol_ncbi_ids[protein_official_symbol].strip().split(';')
        for ncbi_id in ncbi_ids:
            if ncbi_id.strip() == '':
                continue
            if ncbi_id.strip() not in human_protein_ncbi_ids:
                ncbi_id_human_protein[ncbi_id] = set()
                human_protein_ncbi_ids.add(ncbi_id)
            ncbi_id_human_protein[ncbi_id].add(protein_official_symbol)

    human_protein_uniprot_entries = set()
    uniprot_entry_human_protein = {}
    for protein_official_symbol in human_protein_official_symbol_uniprot_entries.keys():
        uniprot_entries = human_protein_official_symbol_uniprot_entries[protein_official_symbol].strip().split(';')
        for uniprot_entry in uniprot_entries:
            if uniprot_entry.strip() == '':
                continue
            if uniprot_entry.strip() not in human_protein_uniprot_entries:
                uniprot_entry_human_protein[uniprot_entry] = set()
                human_protein_uniprot_entries.add(uniprot_entry)
            uniprot_entry_human_protein[uniprot_entry].add(protein_official_symbol)

    mapping_uniprot_geneid_dir = dir_path + '/mapping_uniprot_vs_string_gene/'
    dir_file_names = [f for f in listdir(mapping_uniprot_geneid_dir) if f.endswith('gz')]
    print('dir_file_names: '+str(len(dir_file_names)))
    for tax_name in species_dict.keys():
        tax_id = species_dict[tax_name]

        # temp for del 为了临时目的
        # if tax_id != '9606':
        #     continue

        ncbi_homologous_gene_id_human_gene_id, uniprot_homologous_gene_entry_human_gene_entry = hhgtd.other_species_homologous_gene_names_ids(tax_id)
        if tax_id != '9606':
            protein_ncbi_ids = set(ncbi_homologous_gene_id_human_gene_id.keys())
            protein_uniprot_entries = set(uniprot_homologous_gene_entry_human_gene_entry.keys())
        else:
            protein_ncbi_ids = human_protein_ncbi_ids
            protein_uniprot_entries = human_protein_uniprot_entries

        # print('number: protein_ncbi_ids: ' + str(
        #     len(protein_ncbi_ids)))
        # print('number: protein_uniprot_entries: ' + str(
        #     len(protein_uniprot_entries)))

        # uniprot 与 string ID对应关系
        dir_file_name = ''
        for file_name in dir_file_names:
            if tax_name in file_name:
                dir_file_name = mapping_uniprot_geneid_dir + file_name
                break
        mapping_uniprot_geneid_xlsx = un_gz(dir_file_name)
        data = xlrd.open_workbook(mapping_uniprot_geneid_xlsx)
        table = data.sheet_by_name('Sheet0')
        uniprot_entries = table.col_values(0)[1:]
        string_protein_ids = table.col_values(2)[1:]
        ensembl_protein_identifier_uniprot_gene_entry = {}
        ensembl_protein_identifiers_uniprot = set()
        for i in range(len(uniprot_entries)):
            if string_protein_ids[i].strip() == '':
                continue
            if uniprot_entries[i].strip() not in protein_uniprot_entries:
                continue
            ensembl_protein_identifiers_line = string_protein_ids[i].strip().split(';')
            for ensembl_protein_identifier in ensembl_protein_identifiers_line:
                if ensembl_protein_identifier.strip() == '':
                    continue
                ensembl_protein_identifier = ensembl_protein_identifier.strip().split('.')[1].strip()
                if ensembl_protein_identifier.strip() not in ensembl_protein_identifiers_uniprot:
                    ensembl_protein_identifier_uniprot_gene_entry[ensembl_protein_identifier.strip()] = set()
                    ensembl_protein_identifiers_uniprot.add(ensembl_protein_identifier.strip())
                ensembl_protein_identifier_uniprot_gene_entry[ensembl_protein_identifier.strip()].add(uniprot_entries[i].strip())
        print(tax_name + ' number: ensembl_protein_identifier_uniprot_gene_entry: ' + str(
            len(ensembl_protein_identifier_uniprot_gene_entry.keys())))

        # NCBI 与 string ID对应关系
        ensembl_protein_identifier_ncbi_gene_id = {}
        ensembl_protein_identifiers_ncbi = set()
        gene_info_fp = dir_path + 'ftp_files/gene2ensembl'

        with codecs.open(gene_info_fp, "rb", "utf-8") as input_file:
            for line in input_file:
                if line.strip()[0] == '#':
                    continue
                temp = line.strip().split('	')
                if len(temp) != 7:
                    continue
                if temp[0].strip() != tax_id:
                    continue
                if temp[-1].strip() == '-':
                    continue
                if temp[1].strip() not in protein_ncbi_ids:
                    continue
                ensembl_protein_identifier = temp[-1].strip().split('.')[0].strip()
                if ensembl_protein_identifier not in ensembl_protein_identifiers_ncbi:
                    if ensembl_protein_identifier.strip() == '':
                        continue
                    ensembl_protein_identifier_ncbi_gene_id[ensembl_protein_identifier.strip()] = set()
                    ensembl_protein_identifiers_ncbi.add(ensembl_protein_identifier.strip())
                ensembl_protein_identifier_ncbi_gene_id[ensembl_protein_identifier.strip()].add(temp[1].strip())
        print(tax_name+' number: ensembl_protein_identifier_ncbi_gene_id: ' + str(len(ensembl_protein_identifier_ncbi_gene_id.keys())))

        # 蛋白名称，id与entry，ensembl id
        # human_protein_ensembl_protein_identifiers = ensembl_protein_identifiers_ncbi.union(ensembl_protein_identifiers_uniprot)
        # protein_stringid_pair_score_dir = dir_path + '/species_protein_pair_score_v1/'
        # write_out_fp = protein_stringid_pair_score_dir + 'protein_ensembl_identifier_homologous_gene_' + tax_name + '_symbol'
        # with open(write_out_fp, 'w') as output_file:
        #     for protein_ensembl_identifier in human_protein_ensembl_protein_identifiers:
        #         human_uniprot_entries = set()
        #         human_ncbi_ids = set()
        #         if protein_ensembl_identifier in ensembl_protein_identifier_uniprot_gene_entry.keys():
        #             for entry in ensembl_protein_identifier_uniprot_gene_entry[protein_ensembl_identifier]:
        #                 human_uniprot_entries.add(entry.strip())
        #             if protein_ensembl_identifier in ensembl_protein_identifier_ncbi_gene_id.keys():
        #                 for id in ensembl_protein_identifier_ncbi_gene_id[protein_ensembl_identifier]:
        #                     human_ncbi_ids.add(id)
        #
        #         else:
        #             if protein_ensembl_identifier in ensembl_protein_identifier_ncbi_gene_id.keys():
        #                 for id in ensembl_protein_identifier_ncbi_gene_id[protein_ensembl_identifier]:
        #                     human_ncbi_ids.add(id)
        #
        #         human_protein_names = set()
        #         for entry in human_uniprot_entries:
        #             human_protein_names = human_protein_names | uniprot_entry_human_protein[entry]
        #         for id in human_ncbi_ids:
        #             human_protein_names = human_protein_names | ncbi_id_human_protein[id]
        #         output_file.write(protein_ensembl_identifier + '$' + ';'.join(list(human_uniprot_entries)) + '$' + ';'.join(
        #                     list(human_ncbi_ids)) + '$' + ';'.join(list(human_protein_names)) + '\n')

        # 获取不同物种与人类之间的对应关系，非人类物种的ensemble identifier 与 其同源人类基因
        # if tax_id == '9606':
        #     human_protein_ensembl_protein_identifiers = ensembl_protein_identifiers_ncbi.union(ensembl_protein_identifiers_uniprot)
        #     protein_stringid_pair_score_dir = dir_path + '/species_protein_pair_score_v1/'
        #     write_out_fp = protein_stringid_pair_score_dir + 'protein_ensembl_identifier_homologous_gene_' + tax_name
        #     with open(write_out_fp, 'w') as output_file:
        #         for protein_ensembl_identifier in human_protein_ensembl_protein_identifiers:
        #             human_uniprot_entries = set()
        #             human_ncbi_ids = set()
        #             if protein_ensembl_identifier in ensembl_protein_identifier_uniprot_gene_entry.keys():
        #                 for entry in ensembl_protein_identifier_uniprot_gene_entry[protein_ensembl_identifier]:
        #                     human_uniprot_entries.add(entry.strip())
        #                 if protein_ensembl_identifier in ensembl_protein_identifier_ncbi_gene_id.keys():
        #                     for id in ensembl_protein_identifier_ncbi_gene_id[protein_ensembl_identifier]:
        #                         human_ncbi_ids.add(id)
        #
        #             else:
        #                 if protein_ensembl_identifier in ensembl_protein_identifier_ncbi_gene_id.keys():
        #                     for id in ensembl_protein_identifier_ncbi_gene_id[protein_ensembl_identifier]:
        #                         human_ncbi_ids.add(id)
        #             output_file.write(protein_ensembl_identifier + '$' + ';'.join(list(human_uniprot_entries)) + '$' + ';'.join(
        #                         list(human_ncbi_ids)) + '\n')
        # else:
        #     human_protein_ensembl_protein_identifiers = ensembl_protein_identifiers_ncbi.union(ensembl_protein_identifiers_uniprot)
        #     protein_stringid_pair_score_dir = dir_path + '/species_protein_pair_score_v1/'
        #     write_out_fp = protein_stringid_pair_score_dir + 'protein_ensembl_identifier_homologous_human_gene_' + tax_name
        #     with open(write_out_fp, 'w') as output_file:
        #         for protein_ensembl_identifier in human_protein_ensembl_protein_identifiers:
        #             human_uniprot_entries = set()
        #             human_ncbi_ids = set()
        #             if protein_ensembl_identifier in ensembl_protein_identifier_uniprot_gene_entry.keys():
        #                 for entry in ensembl_protein_identifier_uniprot_gene_entry[protein_ensembl_identifier]:
        #                     human_uniprot_entries.add(uniprot_homologous_gene_entry_human_gene_entry[entry.strip()])
        #                 if protein_ensembl_identifier in ensembl_protein_identifier_ncbi_gene_id.keys():
        #                     for id in ensembl_protein_identifier_ncbi_gene_id[protein_ensembl_identifier]:
        #                         human_ncbi_ids.add(ncbi_homologous_gene_id_human_gene_id[id])
        #
        #             else:
        #                 if protein_ensembl_identifier in ensembl_protein_identifier_ncbi_gene_id.keys():
        #                     for id in ensembl_protein_identifier_ncbi_gene_id[protein_ensembl_identifier]:
        #                         human_ncbi_ids.add(ncbi_homologous_gene_id_human_gene_id[id])
        #             output_file.write(protein_ensembl_identifier + '$' + ';'.join(list(human_uniprot_entries)) + '$' + ';'.join(
        #                         list(human_ncbi_ids)) + '\n')

        # 不同物种蛋白质关系对及其值
        # protein_ensembl_protein_identifiers = ensembl_protein_identifiers_ncbi.union(ensembl_protein_identifiers_uniprot)
        # print(tax_name + ' number: protein_ensembl_protein_identifiers: ' + str(
        #     len(protein_ensembl_protein_identifiers)))
        # protein_stringid_pair_score_dir = dir_path + '/species_protein_pair_score/'
        # protein_stringid_pair_score = read_species_protein_stringid_pair_score(protein_stringid_pair_score_dir, tax_name,protein_ensembl_protein_identifiers)
        # write_out_fp = protein_stringid_pair_score_dir + 'protein_stringid_pair_score_' + tax_name
        # with open(write_out_fp, 'w') as output_file:
        #     for protein_stringid_pair in protein_stringid_pair_score.keys():
        #         output_file.write(protein_stringid_pair + ',' + protein_stringid_pair_score[protein_stringid_pair] + '\n')

# Mapping human homologous relationship pairs into human protein relationship pairs
def mapping_homologous_gene_pair_for_human_protein_pair(dir_path):

    subdir_path = dir_path + '/species_protein_pair_score/'
    file_signs = ['protein_ensembl','protein_stringid']

    # 先准备人类基因id 与ensemble id的对应关系
    dir_file_names = [f for f in listdir(subdir_path) if f.startswith(file_signs[0])]
    print(len(dir_file_names))
    human_gene_entry_ensamble = {}
    human_gene_id_ensamble = {}
    for file_name in dir_file_names:
        if 'human' not in file_name:
            continue
        with codecs.open(subdir_path + file_name, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                if temp[1].strip() != '':
                    for entry in temp[1].strip().split(';'):
                        if entry == '':
                            continue
                        human_gene_entry_ensamble[entry] = temp[0].strip()
                if temp[2].strip() != '':
                    for id in temp[2].strip().split(';'):
                        if id == '':
                            continue
                        human_gene_id_ensamble[id] = temp[0].strip()
    print('human_gene_entry_ensamble: ' + str(len(human_gene_entry_ensamble.keys())))
    print('human_gene_id_ensamble: ' + str(len(human_gene_id_ensamble)))

    # 进行匹配
    for species_name in species_dict.keys():
        if species_name == 'human':
            continue
        print(species_name + '-----------------------------')
        dir_file_names_id = [f for f in listdir(subdir_path) if f.startswith(file_signs[0])]
        print('dir_file_names_id: ' + str(len(dir_file_names_id)))
        homolo_ensamble_id_human_ensamble_id = {}
        temp_homolo_ensamble_ids = set()
        for file_name in dir_file_names_id:
            if species_name not in file_name:
                continue
            with codecs.open(subdir_path + file_name, "rb", "utf-8") as input_file:
                for line in islice(input_file.readlines(), 0, None):
                    temp = line.strip().split('$')
                    if temp[1].strip() != '':
                        for entry in temp[1].strip().split(';'):
                            if entry == '':
                                continue
                            if entry in human_gene_entry_ensamble.keys():
                                if temp[0].strip() not in temp_homolo_ensamble_ids:
                                    homolo_ensamble_id_human_ensamble_id[temp[0].strip()] = set()
                                    temp_homolo_ensamble_ids.add(temp[0].strip())
                                homolo_ensamble_id_human_ensamble_id[temp[0].strip()].add(human_gene_entry_ensamble[entry])
                    if temp[2].strip() != '':
                        for id in temp[2].strip().split(';'):
                            if id == '':
                                continue
                            if id in human_gene_id_ensamble.keys():
                                if human_gene_id_ensamble[id] not in temp_homolo_ensamble_ids:
                                    homolo_ensamble_id_human_ensamble_id[temp[0].strip()] = set()
                                    temp_homolo_ensamble_ids.add(temp[0].strip())
                                homolo_ensamble_id_human_ensamble_id[temp[0].strip()].add(human_gene_id_ensamble[id])
        print('homolo_ensamble_id_human_ensamble_id: ' + str(len(homolo_ensamble_id_human_ensamble_id.keys())))
        print('temp_homolo_ensamble_ids: ' + str(len(temp_homolo_ensamble_ids)))

        dir_file_names_pair = [f for f in listdir(subdir_path) if f.startswith(file_signs[1])]
        print('dir_file_names_pair: ' + str(len(dir_file_names_pair)))
        human_protein_pair_extend_by_homolo_score = {}
        for file_name in dir_file_names_pair:
            if species_name not in file_name:
                continue
            with codecs.open(subdir_path + file_name, "rb", "utf-8") as input_file:
                for line in islice(input_file.readlines(), 0, None):
                    temp = line.strip().split(',')
                    if temp[0].strip() in temp_homolo_ensamble_ids:
                        for human_ensamble in homolo_ensamble_id_human_ensamble_id[temp[0].strip()]:
                            protein_pair = human_ensamble + ',' + temp[0].strip()
                            human_protein_pair_extend_by_homolo_score[protein_pair] = temp[2].strip()

                    if temp[1].strip() in temp_homolo_ensamble_ids:
                        for human_ensamble in homolo_ensamble_id_human_ensamble_id[temp[1].strip()]:
                            protein_pair = human_ensamble + ',' + temp[1].strip()
                            human_protein_pair_extend_by_homolo_score[protein_pair] = temp[2].strip()
        print('human_protein_pair_extend_by_homolo_score: ' + str(len(human_protein_pair_extend_by_homolo_score.keys())))
        print('---------------------------------')

        write_out_fp = subdir_path + 'extend_human_protein_pair_score_by_homolo_' + species_name
        with open(write_out_fp, 'w') as output_file:
            for extend_protein_stringid_pair in human_protein_pair_extend_by_homolo_score.keys():
                output_file.write(extend_protein_stringid_pair + ',' + human_protein_pair_extend_by_homolo_score[extend_protein_stringid_pair] + '\n')

# Confirm that the related proteins
# in the training set have the corresponding protein relationship pairs
def match_train_data_protein_pair_score(dir_path):
    _, _, train_protein_official_symbols_aliases,human_protein_official_symbol_names = hhgtd.prepare_human_protein_symbols_ids()

    subdir_path = dir_path + '/species_protein_pair_score/protein_ensembl_identifier_homologous_gene_human_symbol'
    human_protein_symbol_ensemble_ids = {}
    human_protein_symbols = set()
    with codecs.open(subdir_path, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if temp[-1].strip() == '':
                continue
            if ';' not in temp[-1].strip():
                if temp[0].strip() not in human_protein_symbols:
                    human_protein_symbol_ensemble_ids[temp[-1].strip()] = set()
                    human_protein_symbols.add(temp[-1].strip())
                human_protein_symbol_ensemble_ids[temp[-1].strip()].add(temp[0].strip())
            else:
                for protein_symbol in temp[-1].strip().split(';'):
                    if protein_symbol.strip() == '':
                        continue
                    if protein_symbol.strip() not in human_protein_symbols:
                        human_protein_symbol_ensemble_ids[protein_symbol.strip()] = set()
                        human_protein_symbols.add(protein_symbol.strip())
                    human_protein_symbol_ensemble_ids[protein_symbol.strip()].add(temp[0].strip())
    print('human_protein_symbol_ensemble_ids: ' + str(len(human_protein_symbol_ensemble_ids.keys())))
    miss_human_protein_symbols = set(human_protein_official_symbol_names.keys()).difference(human_protein_symbols)
    print('miss_human_protein_symbols: '+str(len(miss_human_protein_symbols)))

    train_protein_symbol_ensemble_ids = {}
    train_ensemble_ids = set()

    for train_protein in train_protein_official_symbols_aliases.keys():
        flag = 0
        if train_protein.strip() in human_protein_symbols:
            train_protein_symbol_ensemble_ids[train_protein.strip()] = human_protein_symbol_ensemble_ids[train_protein.strip()]
            train_ensemble_ids = train_ensemble_ids | human_protein_symbol_ensemble_ids[train_protein.strip()]
            flag = 1
        else:
            for protein in train_protein_official_symbols_aliases[train_protein]:
                if protein in human_protein_symbols:
                    train_protein_symbol_ensemble_ids[train_protein.strip()] = human_protein_symbol_ensemble_ids[
                        protein.strip()]
                    train_ensemble_ids = train_ensemble_ids | human_protein_symbol_ensemble_ids[protein.strip()]
                    flag = 1
            if flag == 1:
                continue
            for key in human_protein_official_symbol_names.keys():
                proteins = human_protein_official_symbol_names[key.strip()]
                if train_protein in proteins:
                    if key.strip() in human_protein_symbols:
                        train_protein_symbol_ensemble_ids[train_protein.strip()] = human_protein_symbol_ensemble_ids[key.strip()]
                        train_ensemble_ids = train_ensemble_ids | human_protein_symbol_ensemble_ids[key.strip()]
                        flag += 1
                        # if flag >= 2:
                        #     print(train_protein.strip(),key.strip())
        # if flag == 0:
        #     print(train_protein+'----------------------')

    print('train_protein_official_symbols_aliases: ' + str(len(train_protein_official_symbols_aliases.keys())))
    print('train_protein_symbol_ensemble_ids: ' + str(len(train_protein_symbol_ensemble_ids.keys())))
    print('miss_human_protein_symbols_for_train: '+str(len(set(train_protein_official_symbols_aliases.keys()).difference(set(train_protein_symbol_ensemble_ids.keys())))))

    write_out_fp = dir_path + '/train_protein_official_symbols_ensemble_ids.csv'
    if not os.path.exists(write_out_fp):
        with open(write_out_fp, 'w') as output_file:
            for protein_symbol in train_protein_official_symbols_aliases.keys():
                protein_name_aliases = ';'.join(list(train_protein_official_symbols_aliases[protein_symbol]))
                protein_ensemble_ids = ''
                if protein_symbol in train_protein_symbol_ensemble_ids.keys():
                    protein_ensemble_ids = ';'.join(train_protein_symbol_ensemble_ids[protein_symbol])
                output_file.write(protein_symbol + '$' + protein_name_aliases + '$' + protein_ensemble_ids + '\n')

    print('train_ensemble_ids: ' + str(len(train_ensemble_ids)))

    return train_ensemble_ids

# Integrating human protein-related pairs and
# homologous gene related pairs from four other species
def prepare_protein_pair_score(dir_path):

    train_ensemble_ids = match_train_data_protein_pair_score(dir_path)

    subdir_path = dir_path + '/species_protein_pair_score/'
    human_protein_ensemble_ids_fp = subdir_path + 'protein_ensembl_identifier_homologous_gene_human'
    human_protein_ensemble_ids = set()
    with codecs.open(human_protein_ensemble_ids_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            human_protein_ensemble_ids.add(temp[0].strip())
    print('human_protein_ensemble_ids: ' + str(len(human_protein_ensemble_ids)))

    human_protein_pairs_fp = subdir_path + 'protein_stringid_pair_score_human'
    human_protein_pair_score = {}
    human_protein_ensemble_ids_have_pair = set()
    with codecs.open(human_protein_pairs_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            # if int(temp[2]) < 300:
            #     continue
            if (temp[0].strip() in human_protein_ensemble_ids) or (temp[1].strip() in human_protein_ensemble_ids):
                human_protein_pair_score[','.join(temp[:2])] = temp[2].strip()
            human_protein_ensemble_ids_have_pair.add(temp[0].strip())
            human_protein_ensemble_ids_have_pair.add(temp[1].strip())
    print('human_protein_pair_score: '+ str(len(human_protein_pair_score.keys())))

    dir_file_names = [f for f in listdir(subdir_path) if f.startswith('extend_')]
    print('dir_file_names: ' + str(len(dir_file_names)))
    for file_name in dir_file_names:
        extend_protein_pairs_fp = subdir_path + file_name
        with codecs.open(extend_protein_pairs_fp, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                # if int(temp[2]) < 300:
                #     continue
                if (temp[0].strip() in human_protein_ensemble_ids) or (temp[1].strip() in human_protein_ensemble_ids):
                    human_protein_pair_score[','.join(temp[:2])] = temp[2].strip()
                human_protein_ensemble_ids_have_pair.add(temp[0].strip())
                # human_protein_ensemble_ids_have_pair.add(temp[1].strip())
    print('human_protein_pair_score: ' + str(len(human_protein_pair_score.keys())))
    print('human_protein_ensemble_ids_have_pair: ' + str(len(human_protein_ensemble_ids_have_pair)))

    write_out_fp = dir_path + '/human_protein_pair_score.csv'
    if not os.path.exists(write_out_fp):
        with open(write_out_fp, 'w') as output_file:
            for protein_pair in human_protein_pair_score.keys():
                output_file.write(protein_pair + ',' + human_protein_pair_score[protein_pair] + '\n')

    print(len(train_ensemble_ids.difference(human_protein_ensemble_ids_have_pair)))

    return human_protein_ensemble_ids

# Convert the node name to int and convert score to float
def prepare_ppi_pair_score_for_nodevec(dir_path):

    train_ensemble_ids = match_train_data_protein_pair_score(dir_path)

    combined_score = 200

    input_fp = dir_path + '/human_protein_pair_score.csv'
    ensemble_ids = set()
    max_v = 0
    min_v = 100
    with codecs.open(input_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            ensemble_ids.add(temp[0].strip())
            ensemble_ids.add(temp[1].strip())
            if int(temp[2].strip()) > max_v:
                max_v = int(temp[2].strip())
            if int(temp[2].strip()) < min_v:
                min_v = int(temp[2].strip())
    print('ensemble_ids: '+str(len(ensemble_ids)))

    ensemble_id_transfer_node_id = {}
    write_out_fp = dir_path + '/protein_ensemble_id_for_node_mapping.csv'
    if not os.path.exists(write_out_fp):
        node_int = 1
        for ensemble_id in ensemble_ids:
            ensemble_id_transfer_node_id[ensemble_id] = node_int
            node_int += 1
        with open(write_out_fp, 'w') as output_file:
            for key in ensemble_id_transfer_node_id.keys():
                output_file.write(key + ' ' + str(ensemble_id_transfer_node_id[key]) + '\n')
    else:
        with codecs.open(write_out_fp, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(' ')
                ensemble_id_transfer_node_id[temp[0].strip()] = temp[1].strip()
    print('ensemble_id_transfer_node_id: ' + str(len(ensemble_id_transfer_node_id.keys())))

    # 阶段1 过滤score
    write_out_fp = dir_path + '/human_protein_node_pair_by_ensemble_transfer_'+str(combined_score)+'.csv'
    human_protein_ensemble_ids_have_pair_filter = set()
    with open(write_out_fp, 'w') as output_file:
        with codecs.open(input_fp, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                # difference:100-1119-11 394 558,150-1119-11 394 558,200-1542-6 060 855
                if int(temp[2].strip()) < combined_score:
                    continue
                node1 = str(ensemble_id_transfer_node_id[temp[0].strip()])
                node2 = str(ensemble_id_transfer_node_id[temp[1].strip()])
                human_protein_ensemble_ids_have_pair_filter.add(temp[0].strip())
                human_protein_ensemble_ids_have_pair_filter.add(temp[1].strip())
                # score = round((float(temp[2].strip())-min_v) / float(max_v-min_v),4)
                output_file.write(node1 + ' ' + node2 + '\n') #+ ' ' + str(score)
    print('human_protein_ensemble_ids_have_pair_filter: ' + str(len(human_protein_ensemble_ids_have_pair_filter)))

    print('train_ensemble_ids_difference_human_protein_ensemble_ids_have_pair_filter: ' + str(len(train_ensemble_ids.difference(human_protein_ensemble_ids_have_pair_filter))))

    # 阶段2 过滤train ensemble id
    write_out_fp = dir_path + '/human_protein_node_pair_by_ensemble_transfer_train_'+str(combined_score)+'.csv'
    human_protein_ensemble_ids_have_pair_filter = set()
    with open(write_out_fp, 'w') as output_file:
        with codecs.open(input_fp, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                # difference:100-1119-11 394 558,150-1119-11 394 558,200-1542-6 060 855
                if int(temp[2].strip()) < combined_score:
                    continue
                if (temp[0].strip() in train_ensemble_ids) or (temp[1].strip() in train_ensemble_ids):
                    node1 = str(ensemble_id_transfer_node_id[temp[0].strip()])
                    node2 = str(ensemble_id_transfer_node_id[temp[1].strip()])
                    human_protein_ensemble_ids_have_pair_filter.add(temp[0].strip())
                    human_protein_ensemble_ids_have_pair_filter.add(temp[1].strip())
                    score = round((float(temp[2].strip())-min_v) / float(max_v-min_v),4)
                    output_file.write(node1 + ' ' + node2 + ' ' + str(score) + '\n')

    print('human_protein_ensemble_ids_have_pair_filter: ' + str(len(human_protein_ensemble_ids_have_pair_filter)))

    print('train_ensemble_ids_difference_human_protein_ensemble_ids_have_pair_filter: ' + str(
        len(train_ensemble_ids.difference(human_protein_ensemble_ids_have_pair_filter))))

# Output a vector of each protein name
def output_protein_ppi_vector_for_train_data(dir_path):

    ensemble_node_index = {}
    ensemble_index_fp = dir_path + '/protein_ensemble_id_for_node_mapping.csv'
    with codecs.open(ensemble_index_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(' ')
            ensemble_node_index[temp[1].strip()] = temp[0].strip()
    print('ensemble_node_index: '+str(len(ensemble_node_index.keys())))

    ensemble_vector = {}
    node_index_vector_fp = dir_path + '/human_protein_node_vector.csv'
    with codecs.open(node_index_vector_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(' ')
            if len(temp) < 129:
                continue
            ensemble_vector[ensemble_node_index[temp[0].strip()]] = np.array(map(float, temp[1:]), dtype=np.float32)
    print('ensemble_vector: ' + str(len(ensemble_vector.keys())))

    train_protein_symbol_vector = {}
    miss_train_protein_symbols = set()
    train_protein_symbol_vector_fp = dir_path + '/train_protein_official_symbols_ensemble_ids.csv'
    with codecs.open(train_protein_symbol_vector_fp, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            ensemble_ids = temp[2].strip().split(';')
            temp_vector = np.zeros((128))
            flag = 0
            for ensemble_id in ensemble_ids:
                if ensemble_id not in ensemble_vector.keys():
                    continue
                temp_vector += ensemble_vector[ensemble_id]
                flag += 1
            if flag == 0:
                # print(temp[0].strip())
                miss_train_protein_symbols.add(temp[0].strip())
                continue
            values = ','.join([str(x) for x in (temp_vector / float(len(ensemble_ids)))])
            # print(values)
            train_protein_symbol_vector[temp[0].strip()] = values
    print('train_protein_symbol_vector: ' + str(len(train_protein_symbol_vector.keys())))
    print('miss_train_protein_symbols: ' + str(len(miss_train_protein_symbols)))

    write_out_fp = dir_path + '/train_protein_symbol_ppi_vector.csv'
    if not os.path.exists(write_out_fp):
        with open(write_out_fp, 'w') as output_file:
            for train_protein_symbol in train_protein_symbol_vector.keys():
                output_file.write(train_protein_symbol + '$' + train_protein_symbol_vector[train_protein_symbol] + '\n')

    return miss_train_protein_symbols


'''
    2、Read in other PPI data and try the following two ideas
'''

# Since there are about 1000 protein exact data
# in the downloaded data training set, the data collected through API is read in
def read_collect_gene_pair_score_from_api(dir_path):
    miss_train_protein_symbols = output_protein_ppi_vector_for_train_data(dir_path)

    file_path = dir_path + 'crawler_ppi_score_by_api/genes_scores.csv'
    protein_symbols_for_train = set()
    train_ensemble_ids_real = set()

    count = 0
    with codecs.open(file_path, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) != 13:
                continue
            # if (temp[5].strip() == '') or (float(temp[5].strip()) <= 0.3):
            if temp[5].strip() == '':
                continue
            protein_symbols_for_train.add(temp[2].strip())
            protein_symbols_for_train.add(temp[3].strip())
            train_ensemble_ids_real.add(temp[0].strip())
            train_ensemble_ids_real.add(temp[1].strip())
            count += 1
    print('protein_symbols_for_train: ' + str(len(protein_symbols_for_train)))
    print('miss_train_protein_symbols_still: '+str(len(miss_train_protein_symbols.difference(protein_symbols_for_train))))

    _, _, train_protein_official_symbols_aliases,_ = hhgtd.prepare_human_protein_symbols_ids()
    print(len(set(train_protein_official_symbols_aliases.keys()).difference(protein_symbols_for_train)))

    train_ensemble_ids = match_train_data_protein_pair_score(dir_path)
    print(len(set(train_ensemble_ids).difference(train_ensemble_ids_real)))

    print(count)

# Since there are about 1000 gaps in both methods, consider using data from homologous families
def read_cog_between_orthologous_groups(dir_path):

    file_path = dir_path + '/COG.mappings.v10.5.txt'

    train_ensemble_ids = match_train_data_protein_pair_score(dir_path)

    orthologous_group_human_protein_ensembles = {}
    orthologous_groups = set()
    human_protein_ensembles = set()
    with codecs.open(file_path, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) < 4:
                continue
            species_id_protein_symbol = temp[0].strip().split('.')
            if len(species_id_protein_symbol) != 2:
                continue
            if species_id_protein_symbol[0].strip() != '9606':
                continue
            if temp[3].strip() not in orthologous_groups:
                orthologous_group_human_protein_ensembles[temp[3].strip()] = set()
                orthologous_groups.add(temp[3].strip())
            orthologous_group_human_protein_ensembles[temp[3].strip()].add(species_id_protein_symbol[1].strip())
            human_protein_ensembles.add(species_id_protein_symbol[1].strip())

    print('orthologous_groups: ' + str(len(orthologous_groups)))
    print('human_protein_ensembles: ' + str(len(human_protein_ensembles)))

    print('train_protein_official_symbols_mapped: ' + str(len(train_ensemble_ids.intersection(human_protein_ensembles))))

    write_out_fp = dir_path + 'orthologous_group_human_protein_symbols.csv'
    with open(write_out_fp, 'w') as output_file:
        for key in orthologous_group_human_protein_ensembles.keys():
            output_file.write(key + '$' + ';'.join(list(orthologous_group_human_protein_ensembles[key])) + '\n')


if __name__ == '__main__':

    dir_path = ''
    # prepare_ppi_network_for_node2vec(dir_path)
    # suppliment_ppi_network_for_node2vec(dir_path)
    # mapping_required_human_protein_relationship_pair(dir_path)
    # mapping_homologous_gene_pair_for_human_protein_pair(dir_path)
    # match_train_data_protein_pair_score(dir_path)
    # prepare_protein_pair_score(dir_path)
    # prepare_ppi_pair_score_for_nodevec(dir_path)
    output_protein_ppi_vector_for_train_data(dir_path)

    # read_collect_gene_pair_score_from_api(dir_path)
    # read_cog_between_orthologous_groups(dir_path)