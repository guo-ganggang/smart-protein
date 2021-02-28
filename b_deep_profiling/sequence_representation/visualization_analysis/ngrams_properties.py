#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 4/11/2018 3:59 PM
# @Author  : ganggang & xufang
# @Site    : Reproductive Medicine
# @File    : ngrams_properties.py
# @Software: Mining from a specific set of proteins in human sperm

import numpy as np
import collections
from random import *

# mass(average mass of amino acid), Molecular Weight（分子量）,Hydropathy Index,
#  hydrophobicity,charge,and van der Waals volume（等电点PI）//polarity,volume,

# a H1, hydrophobicity; H2, hydrophicility; H3, hydrogen bond; V, volumes of side chains; P1, polarity;
#  P2, polarizability; SASA, solvent-accessible surface area; NCI, net charge index of side chains;
#  MASS, average mass of amino acid.

NGRAM_PROPERTIES = {
    'A':[0.62,-0.5,2,27.5,8.1,0.046,1.181,0.007187,71.0788],
    'C':[0.29,-1,2,44.6,5.5,0.128,1.461,-0.03661,103.1388],
    'D':[-0.9,3,4,40,13,0.105,1.587,-0.02382,115.0886],
    'E':[-0.74,3,4,62,12.3,0.151,1.862,0.006802,129.1155],
    'F':[1.19,-2.5,2,115.5,5.2,0.29,2.228,0.037552,147.1766],
    'G':[0.48,0,2,0,9,0,0.881,0.179052,57.0519],
    'H':[-0.4,-0.5,4,79,10.4,0.23,2.025,-0.01069,137.1411],
    'I':[1.38,-1.8,2,93.5,5.2,0.186,1.81,0.021631,113.1594],
    'K':[-1.5,3,2,100,11.3,0.219,2.258,0.017708,128.1741],
    'L':[1.06,-1.8,2,93.5,4.9,0.186,1.931,0.051672,113.1594],
    'M':[0.64,-1.3,2,94.1,5.7,0.221,2.034,0.002683,131.1986],
    'N':[-0.78,2,4,58.7,11.6,0.134,1.655,0.005392,114.1039],
    'P':[0.12,0,2,41.9,8,0.131,1.468,0.239531,97.1167],
    'Q':[-0.85,0.2,4,80.7,10.5,0.18,1.932,0.049211,128.1307],
    'R':[-2.53,3,4,105,10.5,0.18,1.932,0.049211,156.1875],
    'S':[-0.18,0.3,4,29.3,9.2,0.062,1.298,0.004627,87.0782],
    'T':[-0.05,-0.4,4,51.3,8.6,0.108,1.525,0.003352,101.1051],
    'V':[1.08,-1.5,2,71.5,5.9,0.14,1.645,0.057004,99.1326],
    'W':[0.81,-3.4,3,145.5,5.4,0.409,2.663,0.037977,186.2132],
    'Y':[0.26,-2.3,3,117.3,6.2,0.298,2.368,0.023599,163.1760]
}

# sum up 3grams property
"""
sum_properties = property['A'] + property['B'] + property['C']
"""
def calculate_property(label):
    split_to_char = list(label)
    sum_properties = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0.])
    for char in split_to_char:
        sum_properties += np.array(pick_key(char))
    return sum_properties/3.


"""
list = [
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [2, 3, 4, 5, 6, 7, 8, 9, 10]
       ]
"""
def make_property_list(labels):
    property_list = []
    for label in labels:
        property_list += [calculate_property(label)]
    return property_list

def pick_key(char):
    rand_dict = { 1 : 'N', 2 : 'D', 3 : 'E', 4 : 'Q', 5 : 'L', 6 : 'I'}
    U_stop_codons = [0., 0., 0., 0., 0., 0., 0., 0., 168.053]
    O_stop_codons = [0., 0., 0., 0., 0., 0., 0., 0., 255.313]
    try:
        return NGRAM_PROPERTIES[char]
    #return NGRAM_PROPERTIES[char]
    except:
        if char == 'B':
            return NGRAM_PROPERTIES[rand_dict[randint(1, 2)]]
        elif char == 'Z':
            return NGRAM_PROPERTIES[rand_dict[randint(3, 4)]]
        elif char == 'J':
            return NGRAM_PROPERTIES[rand_dict[randint(5, 6)]]
        elif char == 'X':
            return NGRAM_PROPERTIES[list(NGRAM_PROPERTIES.keys())[randint(0, 19)]]
        elif char == 'U':
            return U_stop_codons
        elif char == 'O':
            return O_stop_codons
