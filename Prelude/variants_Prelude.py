# PRELUDE
# /Memoire/PreludeFugue

# coding=utf-8
from numpy import *
import petl as etl
from re import *
import operator
import glob 
import pandas as pd
import re
import csv

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Group Prelude results in one csv file
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#-----------------
# WILD TYPE
#------------------
# /Memoire/PreludeFugue/ResPreFug-Najma/PreFug_wt
prelude_wt=open('prelude_wt.csv','w')
prelude_wt.write('Uniprot_ACC\taa_nb\taa\tprediction\n')
for fil in glob.glob('*.pre'): 
    name = str.split(fil,'_')[0]
    f=open(fil, 'r')
    lines = f.readlines()
    seq = lines[7].split()[1] #lines[7]: protein sequence
    try:
        seq_pred = lines[9].split()[1]  #lines[9]: 1st lowest energy confo
    except IndexError:
        print(name)
        seq_pred = '0' * len(seq)
    position = 1
    for i in range(0,len(seq)):
        if seq_pred[i] == '0':
            prelude_wt.write(name + '\t' + str(position) + '\t'+ seq[i]+ '\t' + 'NA' +'\n')
        else:
            prelude_wt.write(name + '\t' + str(position) + '\t'+ seq[i]+ '\t' + seq_pred[i]+'\n')
        position = position + 1
    f.close()
prelude_wt.close()

#------------
# MUTANT
#------------
# /Memoire/PreludeFugue/ResPreFug-Najma/PreFug_mut
prelude_mut=open('prelude_mut.csv','w')
prelude_mut.write('Uniprot_ACC\tvariant\taa_nb\taa\tprediction\n')
for fil in glob.glob('*.pre'): 
    name = str.split(fil,'_')[0]
    x = str.split(fil,'_')[2]
    variant = str.split(x,'.')[0]
    f=open(fil, 'r')
    lines = f.readlines()
    seq = lines[7].split()[1]    #lines[7]: protein sequence
    try:
        seq_pred = lines[9].split()[1]   #lines[9]: 1st lowest energy confo
    except IndexError:
        print(name)
        seq_pred = '0' * len(seq)
    position = 1
    for i in range(0,len(seq)):
        if seq_pred[i] == '0':
            prelude_mut.write(name + '\t' + variant + '\t' + str(position) + '\t'+ seq[i]+ '\t' + 'NA' +'\n')
        else:
            prelude_mut.write(name + '\t' + variant + '\t' + str(position) + '\t'+ seq[i]+ '\t' + seq_pred[i]+'\n') 
        position = position + 1
    f.close()
prelude_mut.close()



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Create a new table with only the variants by using the key
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# /Memoire/PreludeFugue

#-----------------
# WILD TYPE
#------------------
# Join prelude table and didagenes table to have the column 'gene_name'
a = pd.read_csv('prelude_wt.csv','\t')
b = pd.read_csv('didagenes.csv','\t')
a = a.fillna('NA')
merged = a.merge(b, on='Uniprot_ACC')
merged.to_csv('_prelude_wt2.csv', index=False)

# Create the key
c = pd.read_csv('_prelude_wt2.csv',',')
c['key'] = ''
c['key'] = c.Gene_name + c.aa + c.aa_nb.map(str)
c = c.fillna('NA')
c.to_csv('_preludekey_wt.csv')

# Create a new table with only the variants by using the key
N=open('_preludekey_wt.csv','r')      
V=open ('didavariantskey.csv','r')
R=open ('_preluderesults_wt.csv','w')
c1 = csv.reader(N)
c2 = csv.reader(V)
c3 = csv.writer(R)
prelude = list(c1) #123'035
variants = list(c2) #364
results_row = ['ID','Gene_name','Variant_effect','prediction','prediction_envt']
c3.writerow(results_row)
for i in range(1,len(prelude)): #on commence à 1 car 1ere ligne ce sont les titres
    for j in range(1,len(variants)):
        if variants[j][5]==prelude[i][6]:
            if variants[j][3]=='missense' or variants[j][3]=='silent' or variants[j][3]=='deletion':
                envt = prelude[i-4][4] + prelude[i-3][4]+ prelude[i-2][4] + prelude[i-1][4] + '-' + prelude[i+1][4] + prelude[i+2][4] + prelude[i+3][4] + prelude[i+4][4]
                results_row = [variants[j][1],variants[j][2],variants[j][3],prelude[i][4],envt]
                c3.writerow(results_row)
            elif variants[j][3]=='insertion':
                envt = prelude[i-3][4]+ prelude[i-2][4]+ prelude[i-1][4] + prelude[i][4] + '-' + prelude[i+1][4] + prelude[i+2][4] + prelude[i+3][4] + prelude[i+4][4]
                results_row = [variants[j][1],variants[j][2],variants[j][3],'NA',envt]
                c3.writerow(results_row)
N.close()
V.close()
R.close()
 

#------------
# MUTANT
#------------
# Join prelude table and _didagenes table to have the column 'gene_name'
a = pd.read_csv('prelude_mut.csv','\t')
b = pd.read_csv('didagenes.csv','\t')
a = a.fillna('NA')
merged = a.merge(b, on='Uniprot_ACC')
merged.to_csv('_prelude_mut2.csv', index=False)

# Create the key
c = pd.read_csv('_prelude_mut2.csv',',')
c['key'] = ''
c['key'] = c.Gene_name + c.aa_nb.map(str) + '_' + c.variant.map(str)
c = c.fillna('NA')
c.to_csv('_preludekey_mut.csv')

# Create a new table with only the variants by using the key
N=open('_preludekey_mut.csv','r')      
V=open ('didavariantskey_mut.csv','r')
R=open ('_preluderesults_mut.csv','w')
c1 = csv.reader(N)
c2 = csv.reader(V)
c3 = csv.writer(R)
prelude = list(c1) #317'389
variants = list(c2) #364
results_row = ['ID','Gene_name','Variant_effect','prediction','prediction_envt']
c3.writerow(results_row)
for i in range(1,len(prelude)):
    for j in range(1,len(variants)):
        if variants[j][5]==prelude[i][7]:
            if variants[j][3]=='missense' or variants[j][3]=='silent':
                envt = prelude[i-4][5] + prelude[i-3][5] + prelude[i-2][5] + prelude[i-1][5] + '-' + prelude[i+1][5] + prelude[i+2][5] + prelude[i+3][5] + prelude[i+4][5]
                results_row = [variants[j][1],variants[j][2],variants[j][3],prelude[i][5],envt]
                c3.writerow(results_row)
            elif variants[j][3]=='deletion':
                envt = prelude[i-4][5] + prelude[i-3][5] + prelude[i-2][5] + prelude[i-1][5] + '-' + prelude[i+1][5] + prelude[i+2][5] + prelude[i+3][5] + prelude[i][5]
                results_row = [variants[j][1],variants[j][2],variants[j][3],'NA',envt]
                c3.writerow(results_row)
            elif variants[j][3]=='insertion':
                pred = prelude[i+1][5] + prelude[i+2][5] #position 59 -> G  +  position 60 -> R
                envt = prelude[i][5] + prelude[i-3][5] + prelude[i-2][5] + prelude[i-1][5] + '-' + prelude[i+3][5] + prelude[i+4][5] + prelude[i+5][5] + prelude[i+6][5]
                results_row = [variants[j][1],variants[j][2],variants[j][3],pred,envt]
                c3.writerow(results_row)
 
N.close()
V.close()
R.close()
 
