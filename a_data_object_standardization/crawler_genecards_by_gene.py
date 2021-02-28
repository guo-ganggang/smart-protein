#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 4/10/2018 10:42 PM
# @Author  : ggguo
# @Site    : 
# @File    : crawler_genecards_by_gene.py
# @Software: PyCharm

###############
###Intermediate process code, for user reference only
###############

from random import randint
import urllib
import random
import time
from bs4 import BeautifulSoup
import codecs
from itertools import islice
import socket
import os
import chardet
import httplib
import requests
import urllib2
import re
from selenium import webdriver
from selenium.common.exceptions import TimeoutException,WebDriverException

# Determine whether the folder exists,
# it will be automatically generated if it does not exist
def createIfNotExists(path):
    if os.path.exists(path):
        pass
        # print "path already exists!"
    else:
        os.mkdir(path)
        # print "dir created"


# Generate the initial URL
def generate_url_by_gene_symbol(data_path):

    # 读入基因名称
    filePath_f = data_path + 'gene_symbols_update.csv'
    gene_symbols = []
    with codecs.open(filePath_f, "rb", "utf-8") as input_file:
        for line in islice(input_file.readlines(), 0, None):
            temp = line.strip()
            gene_symbols.append(temp)
    print(len(gene_symbols))

    url_join_part = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='

    gene_symbols_urls = []
    for gene_symbol in gene_symbols:
        gene_symbols_urls.append(url_join_part + gene_symbol)
    print(len(gene_symbols_urls))

    return gene_symbols,gene_symbols_urls


# Parse out the Gene ID and other data
def paser_html_for_mgi_links(html):
    paser_results = []
    try:
        soup = BeautifulSoup(html, 'lxml')
        # 判断是否检索出错
        is_error = soup.findAll("div", attrs={"class": "error"})
        if len(is_error) != 0:
            paser_results.append('-1')
        else:
            animal_models_for_gene = soup.findAll("div", attrs={"id":"func_animal","class":"gc-subsection"})
            if len(animal_models_for_gene) == 0:
                paser_results.append('0')
            elif len(animal_models_for_gene) >= 2:
                paser_results.append('1')
            else:
                mgi_links = animal_models_for_gene[0].findAll("a", attrs={"class":"gc-ga-link ","target":"_blank"})
                for mgi_link in mgi_links:
                    paser_results.append(mgi_link["href"])

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
    # print(paser_results)
    return paser_results


# Click on the URL to get the HTML
def crawler_html_by_click_url(url):
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


# Parse out the Gene ID and other data
def paser_html_for_mgi_info(html):
    title_links = []
    try:
        soup = BeautifulSoup(html, 'lxml')
        gene_name_mgiID = soup.title.text.strip()
        title_links.append(gene_name_mgiID)
        table = soup.find("table", attrs={"id": 'referenceTable'})
        links = table.findAll("a", attrs={"class": "MP"})
        for link in links:
            # print(link)
            title_links.append(link["href"])

        # search_results_in_first_page = soup.findAll("tr", attrs={"class": 'rprt'})
        # sapiens = ['[Homo sapiens (human)]', '[Mus musculus (house mouse)]']
        # for search_result in search_results_in_first_page:
        #     temp_paser_results = []
        #     geneID = search_result.find("span", attrs={"class": 'gene-id'}).text.strip()
        #     geneID = re.sub('ID:','',geneID).strip()
        #     href = '/gene/' + geneID
        #     geneName_info = search_result.find("a", attrs={"href": href})
        #     # print(geneName_info)
        #     if geneName_info is None:
        #         continue
        #     geneName = geneName_info.text.strip()
        #     if geneName.upper() != gene:
        #         continue
        #     search_results = search_result.findAll("td")
        #     sapiens_judge = ''
        #     for result in search_results:
        #         if sapiens[0] in result.text.strip():
        #             sapiens_judge = sapiens[0]
        #             break
        #         elif sapiens[1] in result.text.strip():
        #             sapiens_judge = sapiens[1]
        #             break
        #         else:
        #             sapiens_judge = 'Null'
        #     if sapiens_judge == 'Null':
        #         continue
        #     temp_paser_results.append(gene) #geneName
        #     temp_paser_results.append(geneID)
        #     temp_paser_results.append(sapiens_judge)
        #     # print(','.join(temp_paser_results))
        #     paser_results.append(temp_paser_results)
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

    return title_links


# GET HTML
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

        request = requests.get(url, headers)
        html = request.text.decode('utf-8')

    except urllib2.URLError, e:
        if hasattr(e, "reason"):
            print u"连接股吧失败,错误原因", e.reason
    except socket.timeout as e:
        print type(e)
    except socket.error as e:
        print type(e)

    rand = randint(1, 3)
    time.sleep(rand)
    return html


def exec_cralwer_GeneCards_for_mgi(filePath):

    '''
        爬虫采集genecards 上基因的MGI数据库的信息link
    '''
    outFilePath_m = filePath + 'gene_symbols_mgi_links_update.csv'
    outFilePath_f = filePath + 'fail_to_collect_genecards.csv'
    gene_symbols, gene_symbols_urls = generate_url_by_gene_symbol(filePath)

    complete_gene_symbols_f = []
    if os.path.exists(outFilePath_f):
        with codecs.open(outFilePath_f, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                complete_gene_symbols_f.append(temp[0].strip())
    print(len(complete_gene_symbols_f))

    complete_gene_symbols_m = []
    if os.path.exists(outFilePath_m):
        with codecs.open(outFilePath_m, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split(',')
                complete_gene_symbols_m.append(temp[0].strip())
    print(len(complete_gene_symbols_m))

    complete_gene_symbols = list(set(complete_gene_symbols_m).union(set(complete_gene_symbols_f)))
    print(len(complete_gene_symbols))

    for i in range(len(gene_symbols_urls)):

        '''
          采集GENECARDS上面mgi数据库的链接信息
        '''
        if gene_symbols[i] in complete_gene_symbols:
            continue
        gc_url = gene_symbols_urls[i]
        html = crawler_html_by_click_url(gc_url)
        paser_results = []
        if html is not None:
            paser_results = paser_html_for_mgi_links(html)
        if len(paser_results) == 0:
            print(gene_symbols[i] + ' fail to paser!')
            with open(outFilePath_f, 'a') as output_file:
                output_file.write(gene_symbols[i] + ',' + '-1' + '\n')
            continue
        elif paser_results[0] == '-1':
            print(gene_symbols[i] + ' fail to paser!')
            with open(outFilePath_f, 'a') as output_file:
                output_file.write(gene_symbols[i] + ',' + '-1' + '\n')
            continue
        elif paser_results[0] == '0':
            with open(outFilePath_m, 'a') as output_file:
                output_file.write(gene_symbols[i] + ',' + '0' + '\n')
            continue
        else:
            with open(outFilePath_m, 'a') as output_file:
                output_file.write(gene_symbols[i] + ',' + ','.join(paser_results) + '\n')
        # '''
        #   采集mgi数据库的链接信息
        # '''
        # for mgi_url in paser_results:
        #     html = None
        #     for k in range(3):
        #         if html is None:
        #             html = crawler_html_by_url(mgi_url)
        #         else:
        #             break
        #         if k > 0:
        #             time.sleep(100)
        #     if html is None:
        #         print(gene_symbols_urls[i] + ' fail to get html!')
        #         with open(outFilePath_f, 'a') as output_file:
        #             output_file.write(gene_symbols[i] + '\n')
        #         continue
        #     paser_results = paser_html_for_mgi_info(html)
        #     if len(paser_results) == 0:
        #         print(mgi_url + ' fail to paser!')
        #         with open(outFilePath_f, 'a') as output_file:
        #             output_file.write(gene_symbols[i] + '\n')
        #         continue
            # with open(outFilePath_h, 'a') as output_file:
            #    output_file.write(','.join(gene_id_species[0:2]) + '\n')


if __name__ == '__main__':

    # filePath = ''
    # search_Entrez_fetch_literatures(filePath)
    # retrieve_string_interaction_network(filePath)
    # gene_name = 'A1BG'  #sys.argv[1]
    # print 'Gene:' + '\t' + gene_name
    # print 'Function:' + '\t' + getfunction(gene_name)

    filePath = 'genecards/'
    exec_cralwer_GeneCards_for_mgi(filePath)
