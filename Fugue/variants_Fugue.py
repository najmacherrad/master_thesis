# FUGUE
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
fugue_wt=open('fugue_wt.csv','w')
fugue_wt.write('Uniprot_ACC\taa_nb\taa\tprediction\tweight\n') #120'718
for fil in glob.glob('*.fug'): 
    name = str.split(fil,'_')[0]
    f=open(fil, 'r')
    lines = f.readlines()
    seq,seq_pred,A,C,B,P,G,E,O = '','','','','','','','',''
    for k in range(13,len(lines),10):
        seq=seq + lines[k].split()[2] #lines[13]: protein sequence and then every 10 lines
        seq_pred = seq_pred + lines[k+1].split()[2] #lines[14]: prediction 
        A = A + lines[k+2].split()[3]
        C = C + lines[k+3].split()[3]
        B = B + lines[k+4].split()[3]
        P = P + lines[k+5].split()[3]
        G = G + lines[k+6].split()[3]
        E = E + lines[k+7].split()[3]
        O = O + lines[k+8].split()[3]
    weight=''
    position = 1
    for i in range(0,len(seq)):
        if A[i]!='0':
            weight = weight + A[i]
        elif C[i]!='0':
            weight = weight + C[i]
        elif B[i]!='0':
            weight = weight + B[i]
        elif P[i]!='0':
            weight = weight + P[i]
        elif G[i]!='0':
            weight = weight + G[i]      
        elif E[i]!='0':
            weight = weight + E[i]
        elif O[i]!='0':
            weight = weight + O[i]
        else:
            weight = weight + '0'
        fugue_wt.write(name + '\t' + str(position) + '\t'+ seq[i]+ '\t' + seq_pred[i]+'\t'+ weight[i] +'\n')
        position = position + 1
    f.close()
    
fugue_wt.close()

#------------
# MUTANT
#------------
# /Memoire/PreludeFugue/ResPreFug-Najma/PreFug_mut
fugue_mut=open('fugue_mut.csv','w')
fugue_mut.write('Uniprot_ACC\tvariant\taa_nb\taa\tprediction\tweight\n') #311'164
for fil in glob.glob('*.fug'): 
    name = str.split(fil,'_')[0]
    x = str.split(fil,'_')[2]
    variant = str.split(x,'.')[0]
    f=open(fil, 'r')
    lines = f.readlines()
    seq,seq_pred,A,C,B,P,G,E,O = '','','','','','','','',''
    for k in range(13,len(lines),10):
        if lines[k].startswith(' >')!=True:
            seq=seq + lines[k].split()[2] #lines[13]: protein sequence and then every 10 lines
            seq_pred = seq_pred + lines[k+1].split()[2] #lines[14]: prediction 
            A = A + lines[k+2].split()[3]
            C = C + lines[k+3].split()[3]
            B = B + lines[k+4].split()[3]
            P = P + lines[k+5].split()[3]
            G = G + lines[k+6].split()[3]
            E = E + lines[k+7].split()[3]
            O = O + lines[k+8].split()[3]
        else:
            break
    weight=''
    position = 1
    for i in range(0,len(seq)):
        if A[i]!='0':
            weight = weight + A[i]
        elif C[i]!='0':
            weight = weight + C[i]
        elif B[i]!='0':
            weight = weight + B[i]
        elif P[i]!='0':
            weight = weight + P[i]
        elif G[i]!='0':
            weight = weight + G[i]      
        elif E[i]!='0':
            weight = weight + E[i]
        elif O[i]!='0':
            weight = weight + O[i]
        else:
            weight = weight + '0'
        fugue_mut.write(name + '\t' + variant + '\t' + str(position) + '\t'+ seq[i]+ '\t' + seq_pred[i]+'\t'+ weight[i]  +'\n')
        position = position + 1
    f.close()
fugue_mut.close()


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Create a new table with only the variants by using the key
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# /Memoire/PreludeFugue

#-----------------
# WILD TYPE
#------------------
# Join fugue table and _didagenes table to have the column 'gene_name'
a = pd.read_csv('fugue_wt.csv','\t')
b = pd.read_csv('didagenes.csv','\t')
a = a.fillna('NA')
merged = a.merge(b, on='Uniprot_ACC')
merged.to_csv('_fugue_wt2.csv', index=False)

# Create the key
c = pd.read_csv('_fugue_wt2.csv',',')
c['key'] = ''
c['key'] = c.Gene_name + c.aa + c.aa_nb.map(str)
c = c.fillna('NA')
c.to_csv('_fuguekey_wt.csv')
                
# Create a new table with only the variants by using the key
N=open('_fuguekey_wt.csv','r')      
V=open ('didavariantskey.csv','r')
R=open ('_fugueresults_wt.csv','w')
c1 = csv.reader(N)
c2 = csv.reader(V)
c3 = csv.writer(R)
fugue = list(c1) #120'719
variants = list(c2) #364
results_row = ['ID','Gene_name','Variant_effect','prediction','weight','prediction_envt','weight_envt']
c3.writerow(results_row)
for i in range(1,len(fugue)): 
    for j in range(1,len(variants)):
        if variants[j][5]==fugue[i][7]:
            if variants[j][3]=='missense' or variants[j][3]=='silent' or variants[j][3]=='deletion':
                envt = fugue[i-4][4] + fugue[i-3][4]+ fugue[i-2][4] + fugue[i-1][4] + '-' + fugue[i+1][4] + fugue[i+2][4] + fugue[i+3][4] + fugue[i+4][4]
                weight_envt = fugue[i-4][5] + fugue[i-3][5]+ fugue[i-2][5] + fugue[i-1][5] + '-' + fugue[i+1][5] + fugue[i+2][5] + fugue[i+3][5] + fugue[i+4][5]
                results_row = [variants[j][1],variants[j][2],variants[j][3],fugue[i][4],fugue[i][5],envt,weight_envt]
                c3.writerow(results_row)
            elif variants[j][3]=='insertion': 
                envt = fugue[i-3][4]+ fugue[i-2][4]+ fugue[i-1][4] + fugue[i][4] + '-' + fugue[i+1][4] + fugue[i+2][4] + fugue[i+3][4] + fugue[i+4][4]
                weight_envt = fugue[i-3][5]+ fugue[i-2][5]+ fugue[i-1][5] + fugue[i][5] + '-' + fugue[i+1][5] + fugue[i+2][5] + fugue[i+3][5] + fugue[i+4][5]
                results_row = [variants[j][1],variants[j][2],variants[j][3],'NA','NA',envt,weight_envt]
                c3.writerow(results_row)
N.close()
V.close()
R.close()

#------------
# MUTANT
#------------
# Join fugue table and _didagenes table to have the column 'gene_name'
a = pd.read_csv('fugue_mut.csv','\t')
b = pd.read_csv('didagenes.csv','\t')
a = a.fillna('NA')
merged = a.merge(b, on='Uniprot_ACC')
merged.to_csv('_fugue_mut2.csv', index=False)

# Create the key
c = pd.read_csv('_fugue_mut2.csv',',')
c['key'] = ''
c['key'] = c.Gene_name + c.aa_nb.map(str) + '_' + c.variant.map(str)
c = c.fillna('NA')
c.to_csv('_fuguekey_mut.csv')

# Create a new table with only the variants by using the key
N=open('_fuguekey_mut.csv','r')      
V=open ('didavariantskey_mut.csv','r')
R=open ('_fugueresults_mut.csv','w')
c1 = csv.reader(N)
c2 = csv.reader(V)
c3 = csv.writer(R)
fugue = list(c1) #295'323
variants = list(c2) #365
results_row = ['ID','Gene_name','Variant_effect','prediction','weight','prediction_envt','weight_envt']
c3.writerow(results_row)
for i in range(1,len(fugue)): 
    for j in range(1,len(variants)):
        if variants[j][5]==fugue[i][8]: #KEY
            if variants[j][3]=='missense' or variants[j][3]=='silent':
                envt = fugue[i-4][5]+fugue[i-3][5]+fugue[i-2][5]+fugue[i-1][5]  + '-'  +fugue[i+1][5]+fugue[i+2][5]+fugue[i+3][5]+fugue[i+4][5]
                weight_envt = fugue[i-4][6]+fugue[i-3][6]+fugue[i-2][6]+fugue[i-1][6]  + '-'  +fugue[i+1][6]+fugue[i+2][6]+fugue[i+3][6]+fugue[i+4][6]
                results_row = [variants[j][1],variants[j][2],variants[j][3],fugue[i][5],fugue[i][6],envt,weight_envt]
                c3.writerow(results_row)
            elif variants[j][3]=='deletion':
                envt = fugue[i-4][5]+fugue[i-3][5]+fugue[i-2][5]+fugue[i-1][5]  + '-' +fugue[i][5]+fugue[i+1][5]+fugue[i+2][5]+fugue[i+3][5]
                weight_envt = fugue[i-4][6]+fugue[i-3][6]+fugue[i-2][6]+fugue[i-1][6]  + '-' +fugue[i][6]+fugue[i+1][6]+fugue[i+2][6]+fugue[i+3][6]
                results_row = [variants[j][1],variants[j][2],variants[j][3],'NA','NA',envt,weight_envt]
                c3.writerow(results_row)
            elif variants[j][3]=='insertion':
                pred = fugue[i+1][5] + fugue[i+2][5] #position 59 -> G  +  position 60 -> R
                we = (int(fugue[i+1][6]) + int(fugue[i+2][6]))/2
                envt = fugue[i-3][5]+fugue[i-2][5]+fugue[i-1][5]+fugue[i][5]+    '-' + fugue[i+3][5]+fugue[i+4][5]+fugue[i+5][5]+fugue[i+6][5]
                weight_envt = fugue[i-3][6]+fugue[i-2][6]+fugue[i-1][6]+fugue[i][6]+    '-' + fugue[i+3][6]+fugue[i+4][6]+fugue[i+5][6]+fugue[i+6][6]
                results_row = [variants[j][1],variants[j][2],variants[j][3],pred,we,envt,weight_envt]
                c3.writerow(results_row)
N.close()
V.close()
R.close()    
    
