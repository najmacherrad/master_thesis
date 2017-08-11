# ELM
# /Memoire/ELM

# coding=utf-8
from numpy import *
import petl as etl
from re import *
import operator
import glob 
import pandas as pd
import re
import csv

import wget
import os
from socket import error as SocketError
import errno
#------------------------------------------------------------------------------------
# Download data from elm website API
#------------------------------------------------------------------------------------
###########
# WILD TYPE
file_fasta_wt = 'DIDAgenes_wt.txt'
fasta = open(file_fasta_wt,'r')
lines = fasta.readlines()
sequences={}
for i in range(0,len(lines)):
    if lines[i].startswith('>'):
        splitline = lines[i].split('|') 
        accessorID = splitline[1]
        sequences[accessorID] = ''
    else:
        sequences[accessorID] = sequences[accessorID] + lines[i].rstrip('\n').rstrip('*')

url = 'http://elm.eu.org/start_search/'
for key in sequences:
    sequence = sequences[key]
    newname = key + '.txt'
    try:
        os.system('python2.7 -m wget ' + url + sequence + ' -o ' + newname)
    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            raise
        pass

# For them with no results I don't know why (maybe too long sequence) -> retry with this
url = 'http://elm.eu.org/start_search/'
newname = 'Q8WXG9' + '.txt'
sequence = sequences['Q8WXG9']
os.system('python2.7 -m wget ' + url + sequence + ' -o ' + newname)

#############3
# MUTANT
file_fasta_wt = 'DIDAgenes_mut.txt'
fasta = open(file_fasta_wt,'r')
lines = fasta.readlines()
sequences={}
for i in range(0,len(lines)):
    if lines[i].startswith('>'):
        splitline = lines[i].split('|') 
        accessorID = splitline[1] + '_' + splitline[3].split('_')[1].rstrip()
        sequences[accessorID] = ''
    else:
        sequences[accessorID] = sequences[accessorID] + lines[i].rstrip('\n').rstrip('*')
        
url = 'http://elm.eu.org/start_search/'
for key in sequences:
    sequence = sequences[key]
    newname = key + '.txt'
    try:
        os.system('python2.7 -m wget ' + url + sequence + ' -o ' + newname)
    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            raise
        pass

# For them with no results I don't know why (maybe too long sequence) -> retry with this
url = 'http://elm.eu.org/start_search/'
newname = 'Q9H251_9' + '.txt'
sequence = sequences['Q9H251_9']
os.system('python2.7 -m wget ' + url + sequence + ' -o ' + newname)
