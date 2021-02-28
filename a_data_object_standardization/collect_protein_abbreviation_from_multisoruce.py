#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 21/10/2018 6:06 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : collect_protein_abbreviation_from_multisoruce.py
# @Software: Mining from a specific set of proteins in human sperm

###############
###Intermediate process code, for user reference only
###############

from Bio import Entrez
from Bio import Medline
import urllib2
import socket
import ssl
from bioservices import UniProt
from os import listdir
import os
import codecs
from itertools import islice
import time
from selenium import webdriver
from selenium.common.exceptions import TimeoutException,WebDriverException
from bs4 import BeautifulSoup
from collections import OrderedDict
import re



import sys

reload(sys)
sys.setdefaultencoding('utf-8')

Entrez.email = "ggguo@smu.edu.sg"

################################################################
## Construct a symbol entry dictionary of proteins
##
################################################################

'''
    1. Collect alias for protein names on Genecards
'''

# Generate the initial URL
def generate_url_by_gene_symbol(data_path):

    # 读入基因名称
    filePath_f = data_path + 'gene_symbols_mgi_links_update.csv'
    gene_symbols = []
    with codecs.open(filePath_f, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            gene_symbols.append(temp[0].strip())
    print(len(gene_symbols))

    url_join_part = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='

    gene_symbols_urls = []
    for gene_symbol in gene_symbols:
        gene_symbols_urls.append(url_join_part + gene_symbol)
    print(len(gene_symbols_urls))

    return gene_symbols,gene_symbols_urls

# Click on the URL to get the HTML
def crawler_html_by_click_url(url,filePath):
    path = filePath + 'geckodriver-win64/geckodriver'
    browser = webdriver.Firefox(executable_path=path)
    browser.set_page_load_timeout(150)

    for i in range(6):
        try:
            browser.get(url)
            time.sleep(2)
            break
        except TimeoutException as msg:
            print u"find factor exception%s" % msg
            time.sleep(100)
        except WebDriverException as msg:
            print u"find factor exception%s" % msg
            time.sleep(100)

    html = browser.page_source
    browser.close()
    browser.quit()
    return html

# Parse data such as gene/ protein nicknames
def paser_html_for_aliases_info(html):
    gene_aliases = []
    try:
        soup = BeautifulSoup(html, 'lxml')
        ul = soup.find("ul", attrs={"class": 'list-unstyled list-spacious'})
        lis = ul.findAll("li")
        # print(len(lis))
        for li in lis:
            if li.find("span") is not None:
                span = li.find("span", attrs={"id": "aliasMainName"})
                alias = span.text.strip().split('                    ')[0]
                # print(alias.strip())
                gene_aliases.append(alias.strip())
            else:
                alias = li.text.strip().split('                    ')[0]
                # print(alias.strip())
                gene_aliases.append(alias.strip())
        # print(len(gene_aliases))

    except KeyError as e:
        print e
    except UnicodeEncodeError as e:
        print e
    except AttributeError  as e:
        print e
    except IOError  as e:
        print e
    except TypeError  as e:
        print e

    return gene_aliases

# Parse data such as gene/ protein nicknames
def paser_html_for_database_geneid(html,gene_entry,database_name='Entrez Gene'):
    database_geneids = []
    try:
        soup = BeautifulSoup(html, 'lxml')
        divs = soup.findAll("div", attrs={"class": 'gc-subsection'})
        print(len(divs))
        for div in divs:
            flag_str = 'External Ids for '+gene_entry.strip()+' Gene'
            if flag_str in div.text.strip():
                ul = div.find("ul", attrs={"class": 'list-inline'})
                lis = ul.findAll("li")
                database_geneid = 'null'
                for li in lis:
                    li_content = re.sub('\n','',li.text.strip())
                    database_geneids.append(li_content)
                    if database_name in li_content: #UniProtKB:
                        database_geneid = li_content.split(': ')[1].strip()
                database_geneids.append(database_geneid)
                print(database_geneid)
                break

    except KeyError as e:
        print e
    except UnicodeEncodeError as e:
        print e
    except AttributeError  as e:
        print e
    except IOError  as e:
        print e
    except TypeError  as e:
        print e

    return database_geneids

# Control the acquisition and analysis process
def exec_cralwer_GeneCards(filePath):

    '''
        爬虫采集genecards 上基因的别称信息
    '''
    outFilePath = filePath + 'protein_entities_standardized_genecards_update.csv'
    gene_symbols, gene_symbols_urls = generate_url_by_gene_symbol(filePath)

    complete_gene_symbols = []
    if os.path.exists(outFilePath):
        with codecs.open(outFilePath, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                complete_gene_symbols.append(temp[0].strip())
    print(len(complete_gene_symbols))

    for i in range(len(gene_symbols_urls)):
        if gene_symbols[i] in complete_gene_symbols:
            continue
        gc_url = gene_symbols_urls[i]
        print(gc_url)
        html = crawler_html_by_click_url(gc_url,filePath)
        gene_aliases = paser_html_for_aliases_info(html)
        with open(outFilePath, 'a') as output_file:
            if len(gene_aliases) != 0:
                output_file.write(gene_symbols[i] + '$' + '$'.join(gene_aliases) + '\n')
            else:
                output_file.write(gene_symbols[i] + '\n')

'''
    2. Collect alias for protein names on NCBI
'''

# Custom parsing text, due to NCBI API parsing bugs
def parse_txt_file(dir_path):
    read_txt_path = dir_path + 'protein_entities_standardized_ncbi.csv'
    key_struct = ['Official Symbol','Official Name','Other Aliases','Other Designations','ID'] #,'Chromosome','Annotation','MIM'
    temp_dict_txt = []
    with open(read_txt_path, 'a') as output_file:
        if os.path.getsize(read_txt_path) == 0:
            output_file.write('$'.join(key_struct) + '\n')
        with codecs.open(dir_path + '(Homo sapiens [Organism])AND HGNC.txt', "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip()
                if temp == '':
                    if len(temp_dict_txt) == 5 or len(temp_dict_txt) == 6:
                        if len(temp_dict_txt)== 5:
                            temp_dict_txt.append('')
                        if temp_dict_txt[1] != '':
                            output_file.write('$'.join(temp_dict_txt[1:]) + '\n')
                        temp_dict_txt = []
                    continue
                if len(temp_dict_txt) == 0:
                    eid = temp.split('. ')
                    if len(eid) != 2:
                        # print(temp)
                        temp_dict_txt.append('')
                    else:
                        temp_dict_txt.append(eid[1].strip())
                        continue
                if len(temp_dict_txt) == 1:
                    osy = temp.split('Official Symbol:')
                    if len(osy) != 2:
                        # print(temp)
                        temp_dict_txt.append('')
                    else:
                        temp_dict_txt.append(osy[1].strip())
                        continue
                if len(temp_dict_txt) == 2:
                    ona = temp.split('Official Name:')
                    if len(ona) != 2:
                        # print(temp)
                        temp_dict_txt.append('')
                    else:
                        temp_dict_txt.append(ona[1].strip())
                        continue
                if len(temp_dict_txt) == 3:
                    oal = temp.split('Other Aliases:')
                    if len(oal) != 2:
                        # print(temp)
                        temp_dict_txt.append('')
                    else:
                        temp_dict_txt.append(oal[1].strip())
                        continue
                if len(temp_dict_txt) == 4:
                    ode = temp.split('Other Designations:')
                    if len(ode) != 2:
                        # print(temp)
                        temp_dict_txt.append('')
                    else:
                        temp_dict_txt.append(ode[1].strip())
                        continue
                # if len(temp_dict_txt) == 5:
                #     cso = temp.split('Chromosome:')
                #     if len(cso) != 2:
                #         # print(temp)
                #         temp_dict_txt.append('')
                #     else:
                #         temp_dict_txt.append(cso[1].strip())
                #         continue
                # if len(temp_dict_txt) == 6:
                #     ata = temp.split('Annotation:')
                #     if len(ata) != 2:
                #         # print(temp)
                #         temp_dict_txt.append('')
                #     else:
                #         temp_dict_txt.append(ata[1].strip())
                #         continue
                # if len(temp_dict_txt) == 7:
                #     mim = temp.split('MIM:')
                #     if len(mim) != 2:
                #         # print(temp)
                #         temp_dict_txt.append('')
                #     else:
                #         temp_dict_txt.append(mim[1].strip())
                #         continue
                id = temp.split('ID:')
                if len(id) == 2 and id[0].strip() == '':
                    temp_dict_txt.append(id[1].strip())
                # print(temp)

# Obtain all human protein appellation data based on API
def search_fetch_entrez_gene(dir_path):

    keyword = '(Homo sapiens [Organism])AND HGNC'
    # 全局搜索
    try:
        handle_fetch = Entrez.esearch(db="Gene", term=keyword, usehistory="y")
    except urllib2.URLError as err:
        print err

    # 搜索
    search_results = Entrez.read(handle_fetch)
    count = int(search_results["Count"])
    print "Found %i results" % count
    batch_size = 1000
    outfilePath = dir_path + keyword + ".txt"
    out_handle = open(outfilePath, "w")
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        try:
            fetch_handle = Entrez.efetch(db="Gene",
                                         rettype="medline", retmode="text",
                                         retstart=start, retmax=batch_size,
                                         webenv=search_results["WebEnv"],
                                         query_key=search_results["QueryKey"])
        except urllib2.URLError as err:
            print err
            continue
        except socket.error:
            errno, errstr = sys.exc_info()[:2]
            if errno == socket.timeout:
                print "There was a timeout"
            else:
                print "There was some other socket error"
            continue
        print "Going to download record %i to %i" % (start + 1, end)
        try:
            data = fetch_handle.read()
        except ssl.SSLError as e:
            print e
            continue
        fetch_handle.close()
        out_handle.write(data)
    out_handle.close()

    parse_txt_file(dir_path)

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


# It was found that there were problems
# in the original collection of gene IDs on Uniprot,
# and the remaining and erroneous ones were supplemented
def crawler_ncbi_geneid_from_genecards(dir_path):


    # 获取出错以及需要补充采集的蛋白质名称列表
    total_gene_symbols_file_path = dir_path + 'genecards/' # genecards 上可以检索出结果的gene symbols
    gene_symbols, gene_symbols_urls = generate_url_by_gene_symbol(total_gene_symbols_file_path)
    protein_entry_id_ncbi,_ = id_protein_entry_match(dir_path)  # uniprot_gene_entry_clean_data.csv
    protein_entry_supplement = list(set(gene_symbols).difference(set(protein_entry_id_ncbi)))
    print(len(protein_entry_supplement))

    # 生成
    url_join_part = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='
    protein_entry_supplement_urls = []
    for gene_symbol in protein_entry_supplement:
        protein_entry_supplement_urls.append(url_join_part + gene_symbol)
        # print(url_join_part + gene_symbol)
    print(len(protein_entry_supplement_urls))

    file_path_u_u = dir_path + 'ncbi/protein_ncbi_entry_unpicking.csv'

    complete_gene_symbols = []
    if os.path.exists(file_path_u_u):
        with codecs.open(file_path_u_u, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                complete_gene_symbols.append(temp[0].strip())
    print(len(complete_gene_symbols))

    for i in range(len(protein_entry_supplement_urls)):
        if protein_entry_supplement[i] in complete_gene_symbols:
            continue
        gc_url = protein_entry_supplement_urls[i]
        # print(gc_url)
        # if gc_url != 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=PLS1':
        #     continue
        html = crawler_html_by_click_url(gc_url,dir_path + 'ncbi/')
        database_geneids = paser_html_for_database_geneid(html,protein_entry_supplement[i])

        with open(file_path_u_u, 'a') as output_file:
            output_file.write(protein_entry_supplement[i]+'$'+'$'.join(database_geneids)+'\n')


'''
    3. Collect nicknames for protein names on Uniprot
'''

# Retrieve the specified protein name data based on the API
def search_fetch_uniprot(dir_path):

    protein_uniprot_entry = []

    # # 获取人类特定蛋白质集合在uniprot 上的蛋白质标识
    # protein_uniprot_entry_dir_path = 'C:/ggguo/1_data_help/uniprot/parse_literatures_human/'
    # fileNameFeature = '.txt'
    # fileNames = [f for f in listdir(protein_uniprot_entry_dir_path) if f.endswith(fileNameFeature)]
    # print(len(fileNames))
    #
    # for fileName in fileNames:
    #     pue = fileName.strip().split('.')[0].strip().split('_')[-1].strip()
    #     protein_uniprot_entry.append(pue)
    # print(len(protein_uniprot_entry))

    with codecs.open(dir_path + 'uniprot/protein_uniprot_entry_unpicking.csv', "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split('$')
            if len(temp) <= 2:
                continue
            if temp[-1].strip() == 'null':
                continue
            protein_uniprot_entry.append(temp[-1].strip())
    print(len(protein_uniprot_entry))

    outFilePath = dir_path + 'uniprot/protein_entities_standardized_uniprot_unpicking_data.csv'
    entry_id = set()
    if os.path.exists(outFilePath):
        with codecs.open(outFilePath, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                entry_id.add(temp[0].strip())
        print(len(entry_id))

    for entry in protein_uniprot_entry:
        if entry in entry_id:
            continue
        u = UniProt() #verbose=False
        data = u.search(entry, frmt="tab",limit=1,columns="id,genes,protein names") #+'+AND+organism:9606'
        try:
            search_results = data.strip().split('\n')[1]
            write_results = search_results.strip().split('	')
            if len(write_results) != 3:
                print(search_results)
                continue
            with open(outFilePath, 'a') as output_file:
                output_file.write('$'.join(write_results) + '\n')
                time.sleep(1)
                print(entry)

        except AttributeError as e:
            print(e)
            continue


# It was found that there were problems
# in the original collection of gene IDs on Uniprot,
# and the remaining and erroneous ones were supplemented
def crawler_uniprot_geneid_from_genecards(dir_path):
    # # 将初始采集的uniprot上的geneID出错的脏数据清洗一下
    # _, id_protein_entry_uniprot = id_protein_entry_match(dir_path) # uniprot_gene_entry_clean_data.csv
    # filePath_u_cd = dir_path + 'protein_entities_standardized/protein_entities_standardized_uniprot_clean_data.csv'
    # filePath_u_dd = dir_path + 'protein_entities_standardized/protein_entities_standardized_uniprot_dirty_data.csv'
    # if not os.path.exists(filePath_u_cd):
    #     id_protein_entry_uniprot_new = set()
    #     with open(filePath_u_cd, 'a') as output_file:
    #         with codecs.open(filePath_u_dd, "rb", "utf-8") as input_file:
    #             for line in islice(input_file.readlines(), 1, None):
    #                 temp = line.strip().split('$')
    #                 if (temp[0].strip() in id_protein_entry_uniprot.keys()) and (temp[0].strip() not in id_protein_entry_uniprot_new):
    #                     recommended_synonyms = temp[1].strip().split(' ')
    #                     if recommended_synonyms[0].strip() == id_protein_entry_uniprot[temp[0].strip()]:
    #                         id_protein_entry_uniprot_new.add(temp[0].strip())
    #                         output_file.write(line.strip() + '\n')
    #                     else:
    #                         if len(recommended_synonyms) > 1:
    #                             for r in range(1,len(recommended_synonyms)):
    #                                 if recommended_synonyms[r].strip() == id_protein_entry_uniprot[temp[0].strip()]:
    #                                     id_protein_entry_uniprot_new.add(temp[0].strip())
    #                                     output_file.write(line.strip() + '\n')
    #                                     # print(temp[0].strip(),str(r), '----------')
    #                                     continue
    #                         # else:
    #                         #     print(temp[0].strip(), 'error----------')
    #                 # else:
    #                 #     print(temp[0].strip(), temp[1].strip().split(' ')[0].strip(),'++++++++++')

    # 获取出错以及需要补充采集的蛋白质名称列表
    total_gene_symbols_file_path = dir_path + 'genecards/' # genecards 上可以检索出结果的gene symbols
    gene_symbols, gene_symbols_urls = generate_url_by_gene_symbol(total_gene_symbols_file_path)
    _, id_protein_entry_uniprot = id_protein_entry_match(dir_path)  # uniprot_gene_entry_clean_data.csv
    protein_entry_supplement = list(set(gene_symbols).difference(set(id_protein_entry_uniprot.values())))
    print(len(protein_entry_supplement))
    uniprot_protein_id = set()
    file_path_u_cd = dir_path + 'protein_entities_standardized/uniprot/protein_entities_standardized_uniprot_clean_data.csv'
    with codecs.open(file_path_u_cd, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split('$')
            uniprot_protein_id.add(temp[0].strip())
    print(len(uniprot_protein_id))
    print(len(list(set(id_protein_entry_uniprot.keys()).difference(uniprot_protein_id))))


    print(len(list(set(gene_symbols).difference(uniprot_protein_id))))

    for id in list(set(id_protein_entry_uniprot.keys()).difference(uniprot_protein_id)):
        if id_protein_entry_uniprot[id] not in protein_entry_supplement:
            protein_entry_supplement.append(id_protein_entry_uniprot[id])
    print(len(protein_entry_supplement))

    # 生成
    url_join_part = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='
    protein_entry_supplement_urls = []
    for gene_symbol in protein_entry_supplement:
        protein_entry_supplement_urls.append(url_join_part + gene_symbol)
        # print(url_join_part + gene_symbol)
    print(len(protein_entry_supplement_urls))

    file_path_u_u = dir_path + 'protein_entities_standardized/uniprot/protein_uniprot_entry_unpicking.csv'

    complete_gene_symbols = []
    if os.path.exists(file_path_u_u):
        with codecs.open(file_path_u_u, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                complete_gene_symbols.append(temp[0].strip())
    print(len(complete_gene_symbols))

    for i in range(len(protein_entry_supplement_urls)):
        if protein_entry_supplement[i] in complete_gene_symbols:
            continue
        gc_url = protein_entry_supplement_urls[i]
        # print(gc_url)
        # if gc_url != 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=PLS1':
        #     continue
        html = crawler_html_by_click_url(gc_url,dir_path + 'protein_entities_standardized/')
        database_geneid = paser_html_for_database_geneid(html,protein_entry_supplement[i])

        with open(file_path_u_u, 'a') as output_file:
            output_file.write(protein_entry_supplement[i]+'$'+'$'.join(database_geneid)+'\n')


'''
    4. The dictionary is constructed according 
    to the official protein entity name - all gene/protein nicknames
'''

# According to the Official Symbol as KEY,
# protein name and gene name as value,
# respectively read the data collected from the three databases
def entity_standardized_dictionary(data_path):

    # 读入待检索的关键词
    # inFilePath = data_path + 'protein_entities_standardized/gene_symbols.csv'
    # search_keywords = set()
    # with codecs.open(inFilePath, "rb", "utf-8") as input_file:
    #     for line in islice(input_file.readlines(), 0, None):
    #         search_keywords.add(line.strip())
    # print('The number of raw gene name '+str(len(search_keywords)))

    inFilePath = data_path + 'protein_entities_standardized/human_official_symbol_protein_names_all_abbreviation_aliases.csv'
    search_keywords = set()
    with codecs.open(inFilePath, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            search_keywords.add(temp[0].strip())
    print('The number of raw gene name '+str(len(search_keywords)))

    protein_entity_standardized_dictionary = {}

    # uniprot
    # 读入补充采集的uniprot protein id 蛋白质名称
    # dir_path_uniprot_gene_entry = data_path + 'uniprot/uniprot_gene_entry_clean_data.csv'
    # dir_path_uniprot_gene_entry_unpicking = data_path + 'uniprot/protein_uniprot_entry_unpicking.csv'
    # id_protein_entry_uniprot = {}

    # with codecs.open(dir_path_uniprot_gene_entry, "rb", "utf-8") as input_file:
    #     for line in islice(input_file.readlines(), 0, None):
    #         temp = line.strip().split(',')
    #         if temp[0].strip() == 'HUMAN':
    #             if temp[3].strip() in id_protein_entry_uniprot.keys():
    #                 print(temp[3].strip() + '-----------')
    #                 continue
    #             id_protein_entry_uniprot[temp[3].strip()] = temp[1].strip()
    # with codecs.open(dir_path_uniprot_gene_entry_unpicking, "rb", "utf-8") as input_file:
    #     for line in islice(input_file.readlines(), 0, None):
    #         temp = line.strip().split('$')
    #         if len(temp) <= 2:
    #             continue
    #         if temp[-1].strip() == 'null':
    #             continue
    #         if temp[-1].strip() in id_protein_entry_uniprot.keys():
    #             if temp[0].strip() != id_protein_entry_uniprot[temp[-1].strip()]:
    #                 print(temp[-1].strip() + '+++++++++++++')
    #                 id_protein_entry_uniprot[temp[-1].strip()] = id_protein_entry_uniprot[temp[-1].strip()]+','+temp[0].strip()
    #         else:
    #             id_protein_entry_uniprot[temp[-1].strip()] = temp[0].strip()
    # print(len(id_protein_entry_uniprot.keys()))

    # protein_entry_offical_entry_id_uniprot = {}

    uniprot_entity_dict = {}
    filePath_u = data_path + 'protein_entities_standardized/protein_entities_standardized_uniprot_human.csv'
    extract_parentheses = re.compile(r' [(](.*?)[)]', re.S)
    extract_bracket = re.compile(r' [[](.*?)[]]', re.S)

    repeat_uniprot_symbols = set()
    with codecs.open(filePath_u, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split('$')
            if len(temp) != 3:
                print(line.strip())
                continue
            # key = temp[0].strip()
            # if key in repeat_uid:
            #     # print(key)
            #     continue
            # repeat_uid.add(key)
            # if key not in id_protein_entry_uniprot.keys():

            gene_symbol_str = re.sub(' ', ';', temp[1].strip())
            # gene_symbol_str = eval(repr(gene_symbol_str).replace('/', ';'))
            gene_names = gene_symbol_str.split(';')
            if len(gene_names) == 0:
                continue
            official_gene_name = ''
            if gene_names[0].strip() not in repeat_uniprot_symbols:
                official_gene_name = gene_names[0].strip()
                repeat_uniprot_symbols.add(gene_names[0].strip())
            #
            # # flag = 0
            #
            # for gene_name in gene_names:
            #     if gene_name in search_keywords:
            #         official_gene_name = gene_name
            #         # flag = 1
            #         break
            # if flag == 0:
            #     print(key.strip()+'?????????????')
            #     continue
            # else:
            #     official_gene_name = id_protein_entry_uniprot[key]

            gene_name_values_str = '$'.join([x.strip() for x in gene_names]).strip()
            protein_name_values_str_b = re.sub('\(Fragment\)','',temp[2].strip()).strip()
            if ' (' not in protein_name_values_str_b:
                if ' [' not in protein_name_values_str_b:
                    protein_name_values_str_a = protein_name_values_str_b
                else:
                    header_section = protein_name_values_str_b.split(' [')[0].strip()
                    tail_sections = re.findall(extract_bracket, protein_name_values_str_b)
                    temp_tail_section = []
                    for tail_section in tail_sections:
                        if ':' not in tail_section:
                            temp_tail_section.append(tail_section)
                        else:
                            temp_tail_section.append(tail_section.strip().split(':')[1].strip())
                    protein_name_values_str_a = header_section + '$' + '$'.join(temp_tail_section)
            else:
                if ' [' not in protein_name_values_str_b:

                    header_section = protein_name_values_str_b.split(' (')[0].strip()
                    tail_sections = re.findall(extract_parentheses, protein_name_values_str_b)
                    temp_tail_sections = []
                    for tail_section in tail_sections:
                        if '(' in tail_section:
                            temp_tail_sections.append(tail_section+')')
                        else:
                            temp_tail_sections.append(tail_section)
                    protein_name_values_str_a = header_section + '$' + '$'.join(temp_tail_sections)
                else:
                    header_section_bracket = protein_name_values_str_b.split(' [')[0].strip()
                    header_section_parentheses = protein_name_values_str_b.split(' (')[0].strip()
                    header_section = ''
                    if len(header_section_bracket) > len(header_section_parentheses):
                        header_section = header_section_parentheses
                    elif len(header_section_bracket) < len(header_section_parentheses):
                        header_section = header_section_bracket
                    else:
                        print('Error!')
                    tail_sections_bracket = re.findall(extract_bracket, protein_name_values_str_b)
                    temp_tail_sections_bracket = []
                    for tail_section_bracket in tail_sections_bracket:
                        if ':' not in tail_section_bracket:
                            temp_tail_sections_bracket.append(tail_section_bracket)
                        else:
                            temp_tail_sections_bracket.append(tail_section_bracket.strip().split(':')[1].strip())

                    tail_sections_parentheses = re.findall(extract_parentheses, protein_name_values_str_b)
                    temp_tail_sections_parentheses = []
                    for tail_section_parentheses in tail_sections_parentheses:
                        if '(' in tail_section_parentheses:
                            temp_tail_sections_parentheses.append(tail_section_parentheses + ')')
                        else:
                            temp_tail_sections_parentheses.append(tail_section_parentheses)
                    protein_name_values_str_a = header_section+'$'+'$'.join(temp_tail_sections_bracket)+'$'+'$'.join(temp_tail_sections_parentheses)

            # uniprot_entity_dict[official_gene_name] = gene_name_values_str +'$'+ protein_name_values_str_a
            if official_gene_name == '':
                uniprot_entity_dict[gene_names[0].strip()] += gene_name_values_str.strip().split(
                    '$') + protein_name_values_str_a.strip().split('$')
            else:
                if ',' not in official_gene_name:
                    uniprot_entity_dict[official_gene_name] = gene_name_values_str.strip().split('$') + protein_name_values_str_a.strip().split('$')
                    # protein_entry_offical_entry_id_uniprot[official_gene_name] = temp[1].strip().split(' ')[
                    #                                                                  0].strip()# + ',' + key
                else:
                    for name in official_gene_name.strip().split(','):
                        uniprot_entity_dict[name] = gene_name_values_str.strip().split('$') + protein_name_values_str_a.strip().split('$')
                    # protein_entry_offical_entry_id_uniprot[name] = temp[1].strip().split(' ')[
                    #                                                                  0].strip()# + ',' + key

    print('uniprot  '+str(len(uniprot_entity_dict.keys())))

    # outFilePath_poi = data_path + 'protein_entities_standardized/uniprot/protein_entry_offical_entry_id_uniprot.csv'
    # with open(outFilePath_poi, 'w') as output_file:
    #     for protein_entry in protein_entry_offical_entry_id_uniprot.keys():
    #         output_file.write(protein_entry + ',' + protein_entry_offical_entry_id_uniprot[protein_entry] + '\n')

    # # genecards
    # genecards_entity_dict = {}
    # filePath_g = data_path + 'protein_entities_standardized/protein_entities_standardized_genecards.csv'
    #
    # with codecs.open(filePath_g, "rb", "utf-8") as input_file:
    #     for line in islice(input_file.readlines(), 1, None):
    #         temp = line.strip().split('$')
    #         # genecards_entity_dict[temp[0].strip()] = '$'.join([x.strip() for x in temp[1:]])
    #         genecards_entity_dict[temp[0].strip()] = [x.strip() for x in temp[1:]]
    # print('genecards  '+str(len(genecards_entity_dict.keys())))

    # ncbi
    ncbi_entity_dict = {}
    filePath_n = data_path + 'protein_entities_standardized/protein_entities_standardized_ncbi.csv'
    with codecs.open(filePath_n, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 1, None):
            temp = line.strip().split('$')
            if len(temp) != 5:
                print('Error '+line.strip())
                continue
            official_name = re.sub(' \[Homo sapiens \(human\)\]','',temp[1].strip())
            other_aliases = temp[2].strip().split(',')
            other_designations = temp[3].strip().split(';')
            various_names = [official_name] + other_aliases + other_designations
            # ncbi_entity_dict[temp[0].strip()] = '$'.join([x.strip() for x in various_names])
            ncbi_entity_dict[temp[0].strip()] = [x.strip() for x in various_names]
    print('ncbi  '+str(len(ncbi_entity_dict.keys())))

    # print('union  '+str(len(uniprot_genecards_union)))
    # print('intersection  '+str(len(list(set(uniprot_genecards_union).intersection(set(ncbi_entity_dict.keys()))))))
    #
    # for gene_name in list(search_keywords.difference(uniprot_genecards_union)):
    #     print(gene_name)

    # # uniprot 与 genecards 取交集
    # uniprot_genecards_union = list(set(uniprot_entity_dict.keys()).union(set(genecards_entity_dict.keys())))
    # uni_genec_protein_entity_standardized_dictionary = {}
    # for ugu in uniprot_genecards_union:
    #     if ugu in uniprot_entity_dict.keys():
    #         uni_genec_protein_entity_standardized_dictionary[ugu] = uniprot_entity_dict[ugu]
    #         if ugu in genecards_entity_dict.keys():
    #             for ge in genecards_entity_dict[ugu]:
    #                 if ge not in uniprot_entity_dict[ugu]:
    #                     uni_genec_protein_entity_standardized_dictionary[ugu].append(ge)
    #     else:
    #         uni_genec_protein_entity_standardized_dictionary[ugu] = genecards_entity_dict[ugu]
    #
    # print('uniprot and genecards union  '+str(len(uni_genec_protein_entity_standardized_dictionary.keys())))

    # uniprot 与ncbi 取并集,生成最终的词典
    ncbi_uniprot_genecards_union = set(uniprot_entity_dict.keys()).union(set(ncbi_entity_dict.keys()))

    for ugu in ncbi_uniprot_genecards_union:
        if ugu in uniprot_entity_dict.keys():
            protein_entity_standardized_dictionary[ugu] = uniprot_entity_dict[ugu]
            if ugu in ncbi_entity_dict.keys():
                protein_entity_standardized_dictionary[ugu] += ncbi_entity_dict[ugu]
        else:
            if ugu in search_keywords:
                protein_entity_standardized_dictionary[ugu] = ncbi_entity_dict[ugu]
    print('ncbi uniprot and genecards union  ' + str(len(protein_entity_standardized_dictionary.keys())))

    # doubt_protein_entiry_pairs = set()
    outFilePath = data_path + 'protein_entities_standardized/protein_entity_standardized_dictionary.csv'
    completed_standardized_dict_gene_symbol_list = set()
    if os.path.exists(outFilePath):
        with codecs.open(outFilePath, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                completed_standardized_dict_gene_symbol_list.add(temp[0])

    with open(outFilePath, 'a') as output_file:
        for key in protein_entity_standardized_dictionary.keys():
            if key in completed_standardized_dict_gene_symbol_list:
                continue
            temp_values = set(protein_entity_standardized_dictionary[key])
            # temp_values = set()
            # for entiry in protein_entity_standardized_dictionary[key]:
            #     if entiry.strip() != key.strip():
            #         temp_values.add(entiry.strip())
            # for entiry_n in temp_values:
            #     if entiry_n in protein_entity_standardized_dictionary.keys():
            #         if (('0' + ',' + entiry_n+','+key) in doubt_protein_entiry_pairs) or (('1' + ',' + entiry_n+','+key) in doubt_protein_entiry_pairs):
            #             doubt_protein_entiry_pairs.add('1' + ',' + key + ',' + entiry_n)
            #         else:
            #             doubt_protein_entiry_pairs.add('0' + ',' + key + ',' + entiry_n)
            output_file.write(key + '$' + str(len(temp_values)) + '$' + '$'.join(list(temp_values)) + '\n')

    # outFilePath_dpe = data_path + 'protein_entities_standardized/doubt_protein_entiry.csv'
    # with open(outFilePath_dpe, 'a') as output_file:
    #     for entiry_pair in doubt_protein_entiry_pairs:
    #         output_file.write(entiry_pair+'\n')

# Clean dictionary error names and duplicate entries
def clean_protein_name_dictionary(data_path):
    file_path_input_dictionary = data_path + 'protein_entities_standardized/protein_entity_standardized_dictionary.csv'
    file_path_output_dictionary = data_path + 'protein_entities_standardized/protein_entity_standardized_dictionary_clean.csv'
    with open(file_path_output_dictionary, 'w') as output_file:
        with codecs.open(file_path_input_dictionary, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 1, None):
                temp = line.strip().split('$')
                if len(temp) < 3:
                    continue
                if temp[0].strip() == '':
                    continue
                clean_protein_names = set()
                for gene_name in temp[2:]:
                    temp_gene_name = re.sub('Uncharacterized protein', '', gene_name.strip())
                    if temp_gene_name.strip() == '':
                        continue
                    clean_protein_names.add(temp_gene_name)
                if len(clean_protein_names) == 0:
                    continue
                if len(clean_protein_names) == 1:
                    if temp[0].strip() == list(clean_protein_names)[0].strip():
                        continue
                output_file.write(temp[0].strip() + '$' + str(len(clean_protein_names)) + '$' + '$'.join(list(clean_protein_names)) + '\n')


if __name__ == '__main__':

    # dir_path = 'protein_entities_standardized/'
    # search_fetch_entrez_gene(dir_path)
    # search_fetch_uniprot(dir_path)
    # exec_cralwer_GeneCards(dir_path)

    dir_path = ''
    entity_standardized_dictionary(dir_path)
    clean_protein_name_dictionary(dir_path)
    # crawler_uniprot_geneid_from_genecards(dir_path)
    # crawler_ncbi_geneid_from_genecards(dir_path)

    # dir_path = 'vector_protein_mesh/'
    # exec_protein_entity_standardized(dir_path)
