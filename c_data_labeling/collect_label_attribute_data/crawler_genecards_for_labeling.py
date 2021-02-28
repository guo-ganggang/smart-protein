#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 21/11/2018 10:10 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : crawler_genecards_for_labeling.py
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


'''
    Collect the disorder module on the GeneCards
'''

# Generate the initial URL
def generate_url_by_gene_symbol(data_path):

    # 读入基因名称
    filePath_f = data_path + 'human_protein_official_symbols_difference_by_sperm_protein.csv'
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

# Click on the URL to get the HTML
def crawler_html_by_click_url(url):

    path = 'geckodriver-win64/geckodriver'
    browser = webdriver.Firefox(executable_path=path)
    browser.set_page_load_timeout(150)

    for i in range(6):
        try:
            browser.get(url)
            time.sleep(2)
            break
        except TimeoutException as msg:
            print u"find factor exception%s" % msg
            time.sleep(50)
        except WebDriverException as msg:
            print u"find factor exception%s" % msg
            time.sleep(50)

    html = browser.page_source
    browser.close()
    browser.quit()
    return html

# Analyze data such as gene/ protein disease information
def paser_html_for_disorder(html):
    disorder_aliases = []
    try:
        soup = BeautifulSoup(html, 'lxml')
        disorder_section = soup.find("section", attrs={"id": 'diseases',"data-ga-label": 'Disorders',"data-section": 'Disorders'})
        trs = disorder_section.findAll("tr", attrs={"class": 'ng-scope',"data-ng-repeat": 'row in displayedRows'})
        for tr in trs:
            tds = tr.findAll("td")
            disorder = tds[0].a.text.strip()
            # print(disorder)
            lis = tds[1].findAll("li")
            aliases = []
            for i in range(len(lis)-2):
                aliases.append(lis[i].span.text.strip())
                # print(lis[i].span.text.strip())
            disorder_aliases.append([disorder,aliases])

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

    return disorder_aliases

# The mRNA differential expression information of gene/protein was analyzed
def paser_html_for_differential_expression(html):
    differential_expression = []
    try:
        soup = BeautifulSoup(html, 'lxml')
        disorder_section = soup.find("section", attrs={"id": 'expression',"data-ga-label": 'Expression',"data-section": 'Expression'})
        div = disorder_section.find("div", attrs={"id": 'differential-expression',"class": 'gc-subsection'})
        header = div.find("div", attrs={"class": 'gc-subsection-header'})
        if header.a.text.strip() == 'GTEx':
            inner = div.find("div", attrs={"class": 'gc-subsection-inner-wrap'})
            differential_expression.append(inner.span.text.strip())
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

    return differential_expression

# Control the acquisition and analysis process
def exec_cralwer(file_path):

    '''
        爬虫采集genecards 上基因的disorder信息
    '''
    out_file_path = file_path + 'collect_disorder_differential_expression_from_genecards_deseq2_v2.csv'
    gene_symbols, gene_symbols_urls = generate_url_by_gene_symbol(file_path)

    complete_gene_symbols = set()
    if os.path.exists(out_file_path):
        with codecs.open(out_file_path, "rb", "utf-8") as input_file:
            for line in islice(input_file.readlines(), 0, None):
                temp = line.strip().split('$')
                complete_gene_symbols.add(temp[0].strip())
    print(len(complete_gene_symbols))


    for i in range(len(gene_symbols_urls)):
        # if gene_symbols[i] != 'GOPC':
        #     continue
        if gene_symbols[i] in complete_gene_symbols:
            continue

        gc_url = gene_symbols_urls[i]
        print(gc_url)
        with open(out_file_path, 'a') as output_file:
            output_file.write(gene_symbols[i] + '$' + '\n')
        html = crawler_html_by_click_url(gc_url)
        if html is None:
            continue
        disorder_aliases = paser_html_for_disorder(html)
        differential_expression = paser_html_for_differential_expression(html)
        if len(disorder_aliases) == 0:
            continue
        with open(out_file_path, 'a') as output_file:
            for result in disorder_aliases:
                output_file.write(gene_symbols[i] + '$' + result[0] + '$' + ';'.join(result[1]) + '$' + '$'.join(differential_expression) + '\n')


if __name__ == '__main__':

    dir_path = ''
    exec_cralwer(dir_path)