#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/12/2018 6:54 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : prepare_homologous_gene_of_popular_species.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

import codecs
from os import listdir

'''
    1. Homologous GENE based on NCBI Gene ID, derived from NCBI
'''

def get_homologous_genes_from_ncbi(dir_path,species_taxid):

    species_taxids = ['9606','10090','7955','10116','7227']

    file_path = dir_path + 'NCBI_homologene.data'
    groupid_taxid_geneid_gene_symbol = {}
    group_ids = set()
    groupid_taxids = {}
    with codecs.open(file_path, "rb", "utf-8") as input_file:
        for line in input_file:
            temp = line.strip().split('	')
            if len(temp) < 5:
                continue
            if temp[1].strip() not in species_taxids:
                continue
            if temp[0].strip() not in group_ids:
                groupid_taxid_geneid_gene_symbol[temp[0].strip()] = []
                groupid_taxids[temp[0].strip()] = []
                group_ids.add(temp[0].strip())
            groupid_taxid_geneid_gene_symbol[temp[0].strip()].append([x.strip() for x in temp[1:4]])
            groupid_taxids[temp[0].strip()].append(temp[1].strip())

    # 44233 all,26400 popular species
    print('groupid_taxid_geneid_gene_symbol: ' + str(len(groupid_taxid_geneid_gene_symbol.keys())))

    human_gene_id_homologous_taxid_geneids = {}
    for group_id in group_ids:
        if ('9606' not in groupid_taxids[group_id]) or (len(groupid_taxids[group_id])<=1):
            continue
        temp_homologous_genes = []
        human_gene_symbol_id = ''
        for taxid_geneid_gene_symbol in groupid_taxid_geneid_gene_symbol[group_id]:
            if taxid_geneid_gene_symbol[0] == '9606':
                human_gene_symbol_id = '_'.join(taxid_geneid_gene_symbol[1:])
            else:
                temp_homologous_genes.append(taxid_geneid_gene_symbol[0:2])
        human_gene_id_homologous_taxid_geneids[human_gene_symbol_id] = temp_homologous_genes

    print('human_gene_id_homologous_taxid_geneids: ' + str(len(human_gene_id_homologous_taxid_geneids.keys())))

    human_gene_id_species_taxids = {}
    human_gene_ids = set()
    for human_gene_symbol_id in human_gene_id_homologous_taxid_geneids.keys():
        species_taxid_geneid_gene_symbols = human_gene_id_homologous_taxid_geneids[human_gene_symbol_id]
        for species_taxid_geneid_gene_symbol in species_taxid_geneid_gene_symbols:
            if species_taxid == species_taxid_geneid_gene_symbol[0].strip():
                human_gene_id = human_gene_symbol_id.strip().split('_')[0].strip()
                if human_gene_id not in human_gene_ids:
                    human_gene_id_species_taxids[human_gene_id] = set()
                    human_gene_ids.add(human_gene_id)
                human_gene_id_species_taxids[human_gene_id].add(species_taxid_geneid_gene_symbol[1].strip())

    print('human_gene_id_species_taxids: ' + str(len(human_gene_id_species_taxids.keys())))

    return human_gene_id_species_taxids

'''
    2. Homologous GENE based on Uniprot Gene Entry, from Inparanoid
'''

def get_homologous_genes_from_inparanoid(dir_path,species_taxid):

    file_names = [f for f in listdir(dir_path) if f.startswith('sqltable.')]
    species_tax_ids_names = {'7227':'D.melanogaster','7955':'D.rerio','10090':'M.musculus','10116':'R.norvegicus'}
    for file_name in file_names:
        if species_tax_ids_names[species_taxid] not in file_name:
            continue
        file_path = dir_path + file_name
        group_ids = set()
        human_species = []
        other_species = []
        human_species_other_species_pairs = []
        with codecs.open(file_path, "rb", "utf-8") as input_file:
            for line in input_file:
                temp = line.strip().split('	')
                if len(temp) < 5:
                    continue
                if temp[3].strip() != '1.000':
                    continue
                if temp[0].strip() not in group_ids:
                    if len(human_species) != 0:
                        if len(other_species) != 0:
                            for hs in human_species:
                                for os in other_species:
                                    human_species_other_species_pairs.append([hs,os])
                    group_ids.add(temp[0].strip())
                    human_species = []
                    other_species = []
                    if temp[2].strip() == 'H.sapiens':
                        human_species.append(temp[4].strip())
                    else:
                        other_species.append(temp[4].strip())
                else:
                    if temp[2].strip() == 'H.sapiens':
                        human_species.append(temp[4].strip())
                    else:
                        other_species.append(temp[4].strip())

        print('human_species_other_species_pairs: ' + str(len(human_species_other_species_pairs)))
        human_other_homologous_species_genes = {}
        human_gene_entries = set()
        for hos_pair in human_species_other_species_pairs:
            if hos_pair[0].strip() not in human_gene_entries:
                human_other_homologous_species_genes[hos_pair[0].strip()] = set()
                human_gene_entries.add(hos_pair[0].strip())
            human_other_homologous_species_genes[hos_pair[0].strip()].add(hos_pair[1].strip())

        print(file_name)
        print('human_other_homologous_species_genes: ' + str(len(human_other_homologous_species_genes.keys())))
        print('-----------------------------------')

        return human_other_homologous_species_genes

'''
    3. Write to a file and save it in a local folder
'''

def exec_species_homologous_gene_id_mapping(dir_path):

    species_taxids = ['9606', '10090', '7955', '10116', '7227']
    ncbi_homologous_gene_out_fp = dir_path + 'results/human_homologous_gene_by_ncbi_id.csv'
    uniprot_homologous_gene_out_fp = dir_path + 'results/human_homologous_gene_by_uniprot_entry.csv'

    for species_taxid in species_taxids[1:]:
        human_gene_ncbi_id_species_taxids = get_homologous_genes_from_ncbi(dir_path, species_taxid)
        human_gene_uniprot_entry_species_taxids = get_homologous_genes_from_inparanoid(dir_path, species_taxid)

        with open(ncbi_homologous_gene_out_fp, 'a') as output_file:
            for ncbi_gene_id in human_gene_ncbi_id_species_taxids.keys():
                homologous_gene_ids = list(human_gene_ncbi_id_species_taxids[ncbi_gene_id])
                output_file.write(species_taxid + '$' + ncbi_gene_id + '$' + ';'.join(homologous_gene_ids) + '\n')

        with open(uniprot_homologous_gene_out_fp, 'a') as output_file:
            for uniprot_gene_id in human_gene_uniprot_entry_species_taxids.keys():
                homologous_gene_ids = list(human_gene_uniprot_entry_species_taxids[uniprot_gene_id])
                output_file.write(species_taxid + '$' + uniprot_gene_id + '$' + ';'.join(homologous_gene_ids) + '\n')




if __name__ == '__main__':

    dir_path = 'text_corpus/homologene/'
    exec_species_homologous_gene_id_mapping(dir_path)
