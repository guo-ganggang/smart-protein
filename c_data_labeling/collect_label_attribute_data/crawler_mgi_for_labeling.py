#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 4/10/2018 10:42 PM
# @Author  : ggguo
# @Site    : 
# @File    : crawler_mgi_for_labeling.py
# @Software: PyCharm

###############
###Intermediate process code, for user reference only
###############

from random import randint
import random
import time
from bs4 import BeautifulSoup
import codecs
from itertools import islice
import socket
import os
import requests
import urllib2
import re
from collections import OrderedDict

# Determine whether the folder exists,
# it will be automatically generated if it does not exist
def createIfNotExists(path):
    if os.path.exists(path):
        pass
        # print "path already exists!"
    else:
        os.mkdir(path)
        # print "dir created"


################################################################
## Tag data were collected according to gene domain
##
################################################################

# Generate the initial URL
def generate_url_by_mgi_link_module(data_path,module='phenotable'):

    # 读入蛋白名称及其对应链接
    filePath = data_path + 'gene_symbols_mgi_links.csv'
    gene_symbols_links = {}
    with codecs.open(filePath, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip().split(',')
            if len(temp) < 2:
                print(line.strip())
                continue
            if temp[1] == '0':
                continue
            gene_symbols_links[temp[0]] = []
            for url in temp[1:]:
                if module=='phenotable':
                    table_url = '/'.join(url.strip().split('/')[:-1] + [module] + [url.strip().split('/')[-1]])
                else:
                    table_url = url.strip()
                gene_symbols_links[temp[0]].append(table_url)

    print(len(gene_symbols_links))

    return gene_symbols_links

# get html
def crawler_html_by_url(url):
    # context = ssl._create_unverified_context()
    user_agents = [
        'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/69.0.3497.100 Safari/537.36',
        'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:50.0) Gecko/20100101 Firefox/50.0',
        'Mozilla/5.0 (Windows NT 10.0; WOW64; Trident/7.0; rv:11.0) like Gecko'
    ]
    user_agent = random.choice(user_agents)
    headers = {
        'Accept': 'text / html, application / xhtml + xml, application / xml;q = 0.9, image / webp, * / *;q = 0.8', \
        # 'Accept-Encoding':'gzip, deflate, sdch',\
        'Accept-Language': 'zh-CN, zh;q = 0.8',
        'Cache-Control': 'max-age = 0',
        'Upgrade-Insecure-Requests': '1',
        'User-Agent': user_agent
    }
    html = None
    try:

        request = requests.get(url, headers,timeout=10)
        html = request.text.decode('utf-8')

    except urllib2.URLError, e:
        if hasattr(e, "reason"):
            print u"连接股吧失败,错误原因", e.reason
    except socket.timeout as e:
        print type(e)
    except socket.error as e:
        print type(e)
    except Exception as e:
        print type(e)

    rand = randint(1, 3)
    time.sleep(rand)
    return html

# Parses out data such as phenotypes on MGI HTML
def paser_mgi_html_for_phenotypes_info(html):
    parse_results = {}
    try:
        soup = BeautifulSoup(html, 'lxml')

        # gene_name_mgiID = soup.title.text.strip()
        # print(gene_name_mgiID)

        # references = []
        # refer_table = soup.find("table", attrs={"id": 'referenceTable'})
        # refer_links = refer_table.findAll("a", attrs={"class": "MP"})
        # for link in refer_links:
        #     references.append(link["href"])

        # 类别
        pheno_table = soup.find("table", attrs={"id": "phenotable_id", "class": "phenotable"})
        pheno_table_trs = pheno_table.findAll("tr")
        if len(pheno_table_trs) > 2:

            table_header = pheno_table_trs[0].findAll("a")
            header_link = OrderedDict()
            cell_size = []
            if len(table_header) == 0:
                print('table header is null, plase check it, maybe error!=============')
            else:
                header_name = pheno_table_trs[0].find("th", attrs={"id": 'phenoSystemTH'}).text.strip()
                # print(header_name)
                th_header = pheno_table_trs[0].findAll("th")
                for h in range(1, len(th_header)):
                    header = th_header[h].find('a')
                    cell_size.append(int(float(th_header[h]['colspan'].strip())))
                    header_link[header.text.strip()] = [header['href'].strip(), th_header[h]['colspan'].strip()]

            # for key in header_link.keys():
            #     print(key,header_link[key])

            # 性别
            table_subheader = pheno_table_trs[1].findAll("img")
            # subheader_link = OrderedDict()
            subheader_link = []
            sex_cell_size = []
            if len(table_subheader) != 0:
                table_subheader_ths = pheno_table_trs[1].findAll("th")
                subheader_name = table_subheader_ths[0].findAll('span')[-1].text.strip()
                # print(subheader_name)
                for i in range(1, len(table_subheader_ths)):
                    th = table_subheader_ths[i]
                    img = th.find("img")
                    if img is None:
                        # subheader_link['null']=['null',th['colspan'].strip()
                        subheader_link.append('null')
                        sex_cell_size.append(int(float(th['colspan'].strip())))
                        continue
                    symbol = re.sub('\.svg', '', img['src'].strip().split('/')[-1])
                    symbol = re.sub('_', ' ', symbol)
                    # subheader_link[symbol.strip()]=[img['src'].strip(),th['colspan'].strip()]
                    subheader_link.append(symbol.strip())
                    sex_cell_size.append(int(float(th['colspan'].strip())))

            # for key in subheader_link.keys():
            #     print(key,subheader_link[key])

            # 来源
            source_row = []
            is_source_row = 2
            if pheno_table_trs[2]['id'].strip() == 'sourceRow':
                is_source_row = 3
                table_source_ths = pheno_table_trs[2].findAll("th")
                source_name = table_source_ths[0].findAll('span')[-1].text.strip()
                # print(source_name)
                for i in range(1, len(table_source_ths)):
                    th = table_source_ths[i]
                    div = th.find("div")
                    if div is None:
                        source_row.append('null')
                        continue
                    source_push = re.sub(' - ', '-', div.text.strip())
                    source_row.append(source_push)
            # print(source_row)

            # 指标
            id_class_symbols = OrderedDict()
            for j in range(is_source_row, len(pheno_table_trs)):
                tr = pheno_table_trs[j]
                if 'id' in tr.attrs:
                    id_class_symbols[tr['id'].strip().split('_')[0].strip()] = []
                    continue
                subid_value = []
                tds = tr.findAll('td')
                subid = tds[0].find('span').text.strip()
                subid_value.append(subid)
                for k in range(1, len(tds)):
                    if tds[k].find('a') is not None:
                        subid_value.append(tds[k].find('a')['href'])
                    else:
                        subid_value.append('null')
                id_class_symbols[id_class_symbols.keys()[-1]].append(subid_value)

            # 性别对照
            create_sex_index = []
            if len(subheader_link) != 0:
                count_class_index = len(id_class_symbols[id_class_symbols.keys()[0]][0][1:])
                # print(count_class_index,len(subheader_link))
                if count_class_index != len(subheader_link):
                    if len(sex_cell_size) == 1:
                        for k in range(count_class_index):
                            create_sex_index.append(subheader_link[0])
                    else:
                        if min(sex_cell_size) != max(sex_cell_size):
                            for j in range(len(sex_cell_size)):
                                size = sex_cell_size[j]
                                for k in range(size / min(sex_cell_size)):
                                    create_sex_index.append(subheader_link[j])
                        else:
                            for j in range(len(sex_cell_size)):
                                for k in range(count_class_index / len(sex_cell_size)):
                                    create_sex_index.append(subheader_link[j])
                else:
                    create_sex_index = subheader_link

            # 匹配
            if len(id_class_symbols.keys()) == 0:
                print('Parse error! There are no main category--------')
            else:
                for main_class in id_class_symbols.keys():
                    for sub_class in id_class_symbols[main_class]:
                        count_category = len(sub_class) - 1
                        if len(header_link.keys()) == count_category:
                            match_item = []
                            for i in range(len(header_link.keys())):
                                if sub_class[1:][i] == 'null':
                                    continue
                                match_item.append(header_link.keys()[i])
                                match_item.append(sub_class[1:][i])
                                if len(create_sex_index) != 0:
                                    if create_sex_index[i] != 'null':
                                        match_item.append(create_sex_index[i])
                                # if len(source_row) != 0:
                                #     match_item.append(source_row[i])
                            # print(main_class,sub_class[0],match_item)
                            parse_results[main_class + '$' + sub_class[0]] = match_item
                        else:
                            # 选出尺寸最小的值，用比最小的大的除以最小值
                            create_category = []
                            if len(cell_size) == 1:
                                for k in range(count_category):
                                    create_category.append(header_link.keys()[0])
                            else:
                                if min(cell_size) != max(cell_size):
                                    for j in range(len(cell_size)):
                                        size = cell_size[j]
                                        for k in range(size / min(cell_size)):
                                            create_category.append(header_link.keys()[j])
                                else:
                                    for j in range(len(cell_size)):
                                        for k in range(count_category / len(cell_size)):
                                            create_category.append(header_link.keys()[j])
                            if len(create_category) == count_category:
                                match_item = []
                                for i in range(len(create_category)):
                                    if sub_class[1:][i] == 'null':
                                        continue
                                    match_item.append(create_category[i])
                                    match_item.append(sub_class[1:][i])
                                    if len(create_sex_index) != 0:
                                        if create_sex_index[i] != 'null':
                                            match_item.append(create_sex_index[i])
                                    # if len(source_row) != 0:
                                    #     if source_row[i] != 'null':
                                    #         match_item.append(source_row[i])
                                # print(main_class, sub_class[0], match_item)
                                parse_results[main_class + '$' + sub_class[0]] = match_item
                            else:
                                print('Parse error! There are no enough main category*********')


    except KeyError as e:
        print e
    except UnicodeEncodeError as e:
        print e
    except AttributeError  as e:
        print e
    except IOError  as e:
        print e
    except TypeError as e:
        print e
    except requests.exceptions.ChunkedEncodingError as e:
        print e
    except requests.exceptions.ConnectionError as e:
        print e
    return parse_results

# Parses out data such as phenotypes on MGI HTML
def paser_mgi_html_for_mutation_description_info(html):
    parse_results = ''
    try:
        soup = BeautifulSoup(html, 'lxml')

        # 类别
        mutation_description = soup.find("span", attrs={"id": "mutationDescription"}).text.strip()
        # print(mutation_description)
        mutation_description = re.sub('\n','',mutation_description)
        mutation_description = re.sub('		      ', '', mutation_description)
        mutation_description = re.sub(':', '', mutation_description)
        parse_results = mutation_description.strip()
        # print(parse_results)

    except KeyError as e:
        print e
    except UnicodeEncodeError as e:
        print e
    except AttributeError  as e:
        print e
    except IOError  as e:
        print e
    except TypeError as e:
        print e
    except requests.exceptions.ChunkedEncodingError as e:
        print e
    except requests.exceptions.ConnectionError as e:
        print e
    return parse_results

# Control crawler execution
def exec_cralwer_mgi(file_path,module='phenotable'):

    '''
        爬MGI数据库link的信息
    '''

    gene_symbols_links = generate_url_by_mgi_link_module(file_path,module)
    if module=='phenotable':
        outFilePath_h = file_path + 'mgi_mouse_gene_knockout_protein_info.csv'
        outFilePath_f = file_path + 'fail_mgi_mouse_gene_knockout_protein_info.csv'
    else:
        outFilePath_h = file_path + 'mgi_mouse_gene_knockout_protein_info_1.csv'
        outFilePath_f = file_path + 'fail_mgi_mouse_gene_knockout_protein_info_1.csv'

    skip_protein_name = set()
    if os.path.isfile(outFilePath_h):
        with codecs.open(outFilePath_h, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                skip_protein_name.add(temp[0].strip())

    for key in gene_symbols_links.keys():
        if key in skip_protein_name:
            continue
        for mgi_url in gene_symbols_links[key]:
            print(mgi_url)
            mgi_id = mgi_url.split('/')[-1].strip()
            html = None
            for k in range(3):
                if html is None:
                    html = crawler_html_by_url(mgi_url)
                else:
                    break
                if k > 0:
                    time.sleep(50)
            if html is None:
                with open(outFilePath_f, 'a') as output_file:
                    output_file.write(mgi_id + '\n')
                continue
            if module=='phenotable':
                parse_results = paser_mgi_html_for_phenotypes_info(html)
                if len(parse_results) == 0:
                    with open(outFilePath_f, 'a') as output_file:
                        output_file.write(mgi_id + '\n')
                    continue
                with open(outFilePath_h, 'a') as output_file:
                    for class_key in parse_results.keys():
                        output_file.write(
                            key + '$' + mgi_id + '$' + class_key + '$' + '!!'.join(parse_results[class_key]) + '\n')
            else:
                parse_results = paser_mgi_html_for_mutation_description_info(html)
                if parse_results == '':
                    with open(outFilePath_f, 'a') as output_file:
                        output_file.write(mgi_id + '\n')
                    continue
                with open(outFilePath_h, 'a') as output_file:
                    output_file.write(key + '$' + mgi_id + '$' + parse_results + '\n')


if __name__ == '__main__':

    file_path = 'mgi/'
    # generate_url_by_gene_symbol(file_path)
    module = 'mutationDescription'
    exec_cralwer_mgi(file_path,module)