#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 14/11/2018 11:41 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : ebi_api_fetch.py
# @Software: Mining from a specific set of proteins in human sperm

from random import randint
import random
import time
import socket
import requests
import urllib2

import sys
reload(sys)
sys.setdefaultencoding( "utf-8" )

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
            print "Connection failed, error cause: ", e.reason
    except socket.timeout as e:
        print type(e)
    except socket.error as e:
        print type(e)
    except Exception as e:
        print type(e)

    rand = randint(1,2)
    time.sleep(rand)
    return html