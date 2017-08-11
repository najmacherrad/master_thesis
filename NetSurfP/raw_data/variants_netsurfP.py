#NetsurfP
# /Memoire/NetSurfP

# coding=utf-8
from numpy import *
import petl as etl
import re
import operator
import pandas as pd
import csv

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Modify netsurfp results files
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# If results obtained online : 
# - copy and past the results in a text file
# - delete the Netsurfp header and the 'explain output' at the end
# - convert the file into a text file
# - change the extension into .csv and when you open it use the space separator and fusion the separators
# - save the file

#------------
# WILD TYPE
#------------
file_netsurfp_wt = 'NetSurfP_DIDAgenes_wt.csv'

# Add name to column and extract Uniprot ACC
table1 = etl.fromcsv(file_netsurfp_wt, delimiter=' ')
table2 = etl.pushheader(table1, ['Class_assignment','aa','Sequence_name','aa_nb','RSA','ASA',
                                 'Z_score','Proba_a_helix', 'Proba_B_strand','Proba_Coil'])
table3 = etl.split(table2, 'Sequence_name', '_', ['sp', 'Uniprot_ACC','Uniprot_KB','HUMAN'])
table4 = etl.cutout(table3,'sp','HUMAN')
table5 = etl.cut(table4,'Uniprot_ACC','Uniprot_KB','aa','aa_nb','Class_assignment','RSA','ASA',
                 'Z_score','Proba_a_helix', 'Proba_B_strand','Proba_Coil')
etl.totsv(table5, '_netsurfp_wt.csv') #157'384 lignes

'''
# Create the key in didavariants table
d = pd.read_csv('_didavariants.csv','\t')
d = d.fillna('none')
for i in range(0,len(d['Protein_change'])):
    if (d['Protein_change'][i] != 'none'):
        objet = re.search( r'[A-Z]{1}[0-9]+', d['Protein_change'][i])
        d['Protein_change'][i]=objet.group()
d['key'] = ''
d['key'] = d.Gene_name + d.Protein_change
d.to_csv('_didavariantskey.csv',index=False) #364 lines
'''

#------------
# MUTANT
#------------
file_netsurfp_mut = 'NetSurfP_DIDAgenes_wt.csv'

# Add name to column and extract Uniprot ACC
table1 = etl.fromcsv(file_netsurfp_mut, delimiter='\t')
table2 = etl.pushheader(table1, ['Class_assignment','aa','Sequence_name','aa_nb','RSA','ASA',
                                 'Z_score','Proba_a_helix', 'Proba_B_strand','Proba_Coil'])
table3 = etl.split(table2, 'Sequence_name', '_', ['sp', 'Uniprot_ACC','Uniprot_KB','HUMAN','v','variant'])
table4 = etl.cutout(table3,'sp','HUMAN','v') 
table5 = etl.cut(table4,'Uniprot_ACC','Uniprot_KB','variant','aa','aa_nb','Class_assignment',
                 'RSA','ASA','Z_score','Proba_a_helix', 'Proba_B_strand','Proba_Coil')
etl.totsv(table5, '_netsurfp_mut.csv') #399'481 lignes

'''
# Create the key in didavariants table
d = pd.read_csv('_didavariants.csv','\t')
d = d.fillna('none')
for i in range(0,len(d['Protein_change'])):
    if (d['Protein_change'][i] != 'none'):
        objet = re.search( r'[A-Z]{1}[0-9]+', d['Protein_change'][i])
        d['Protein_change'][i]=objet.group()
for i in range(0,len(d['Protein_change'])):
    if (d['Protein_change'][i] != 'none'):
        objet = re.search( r'[0-9]+', d['Protein_change'][i])  
        d['Protein_change'][i]=objet.group()
d['key'] = ''
d['key'] = d.Gene_name + d.Protein_change.map(str) + '_' + d.ID.map(str)
d.to_csv('_didavariantskey_mut.csv') #364 lines
'''

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Create a new table with only the variants by using the key
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#------------
# WILD TYPE
#------------
# Join _netsurfp table and _didagenes table to have the column 'gene_name'
a = pd.read_csv('_netsurfp_wt.csv','\t')
b = pd.read_csv('_didagenes.csv','\t')
b = b.dropna(axis=1)
c = a.merge(b, on='Uniprot_ACC')
c['key'] = ''
c['key'] = c.Gene_name + c.aa + c.aa_nb.map(str)
c.to_csv('_netsurfpkey_wt.csv',index=False)

#PROBLEME !!!!!!!!!!!!!
N = open('_netsurfpkey_wt.csv','r')
V = open('_didavariantskey.csv','r')
R = open('_netsurfpresults_wt.csv','w')
c1 = csv.reader(N)
c2 = csv.reader(V)
c3 = csv.writer(R)
netsurfp = list(c1)
variants = list(c2)
results_row = ['ID','Gene_name','Variant_effect','Class_assignment','RSA','RSA_envt',
               'Z_score','Z_score_envt','ASA','Proba_a_helix','Proba_B_strand','Proba_Coil']
c3.writerow(results_row)
for i in range(1,len(netsurfp)): #[13] = key est à la 14e position
    for j in range(1,len(variants)): #[5] = key est à la 6e position
        if variants[j][5]==netsurfp[i][13]:
            if variants[j][3]=='deletion' or variants[j][3]=='missense' or variants[j][3]=='silent':
                rsaenvt = (float(netsurfp[i-4][6])+float(netsurfp[i-3][6])+float(netsurfp[i-2][6])+float(netsurfp[i-1][6]) + float(netsurfp[i+1][6])+float(netsurfp[i+2][6])+float(netsurfp[i+3][6])+float(netsurfp[i+4][6]))/8
                zscoreenvt = (float(netsurfp[i-4][8])+float(netsurfp[i-3][8])+float(netsurfp[i-2][8])+float(netsurfp[i-1][8]) + float(netsurfp[i+1][8])+float(netsurfp[i+2][8])+float(netsurfp[i+3][8])+float(netsurfp[i+4][8]))/8
                results_row = [variants[j][1],variants[j][2],variants[j][3],netsurfp[i][5],netsurfp[i][6],rsaenvt,netsurfp[i][8],zscoreenvt,netsurfp[i][7],netsurfp[i][9],netsurfp[i][10],netsurfp[i][11]]
                c3.writerow(results_row)
            elif variants[j][3]=='insertion':
                rsaenvt = (float(netsurfp[i-3][6])+float(netsurfp[i-2][6])+float(netsurfp[i-1][6])+float(netsurfp[i][6]) + float(netsurfp[i+1][6])+float(netsurfp[i+2][6])+float(netsurfp[i+3][6])+float(netsurfp[i+4][6]))/8
                zscoreenvt = (float(netsurfp[i-3][8])+float(netsurfp[i-2][8])+float(netsurfp[i-1][8])+float(netsurfp[i][8]) + float(netsurfp[i+1][8])+float(netsurfp[i+2][8])+float(netsurfp[i+3][8])+float(netsurfp[i+4][8]))/8
                results_row = [variants[j][1],variants[j][2],variants[j][3],'NA','NA',rsaenvt,'NA',zscoreenvt,'NA','NA','NA','NA']
                c3.writerow(results_row)

N.close()
V.close()
R.close()


#------------
# MUTANT
#------------
# Join _netsurfp table and _didagenes table to have the column 'gene_name'
a = pd.read_csv('_netsurfp_mut.csv','\t')
b = pd.read_csv('_didagenes.csv','\t')
b = b.dropna(axis=1)
c = a.merge(b, on='Uniprot_ACC')
c['key'] = ''
c['key'] = c.Gene_name + c.aa_nb.map(str) + '_' + c.variant.map(str)
c.to_csv('_netsurfpkey_mut.csv')


N = open('_netsurfpkey_mut.csv','r')
V = open('_didavariantskey_mut.csv','r')
R = open('_netsurfpresults_mut.csv','w')
c1 = csv.reader(N)
c2 = csv.reader(V)
c3 = csv.writer(R)
netsurfp = list(c1)
variants = list(c2)
results_row = ['ID','Gene_name','Variant_effect','Class_assignment','RSA','RSA_envt','Z_score','Z_score_envt','ASA','Proba_a_helix','Proba_B_strand','Proba_Coil']
c3.writerow(results_row)
for i in range(1,len(netsurfp)): #[14] = key est à la 15e position
    for j in range(1,len(variants)): #[5] = key est à la 6e position
        if variants[j][5]==netsurfp[i][14]:
            if variants[j][3]=='silent' or variants[j][3]=='missense':
                rsaenvt = (float(netsurfp[i-4][7])+float(netsurfp[i-3][7])+float(netsurfp[i-2][7])+float(netsurfp[i-1][7]) + float(netsurfp[i+1][7])+float(netsurfp[i+2][7])+float(netsurfp[i+3][7])+float(netsurfp[i+4][7]))/8
                zscoreenvt = (float(netsurfp[i-4][9])+float(netsurfp[i-3][9])+float(netsurfp[i-2][9])+float(netsurfp[i-1][9]) + float(netsurfp[i+1][9])+float(netsurfp[i+2][9])+float(netsurfp[i+3][9])+float(netsurfp[i+4][9]))/8
                results_row = [variants[j][1],variants[j][2],variants[j][3],netsurfp[i][6],netsurfp[i][7],rsaenvt,netsurfp[i][9],zscoreenvt,netsurfp[i][8],netsurfp[i][10],netsurfp[i][11],netsurfp[i][12]]
                c3.writerow(results_row)
            elif variants[j][3]=='insertion':
                classassi = netsurfp[i+1][6] + netsurfp[i+2][6]
                rsa = (float(netsurfp[i+1][7]) + float(netsurfp[i+2][7]))/2
                asa = (float(netsurfp[i+1][8]) + float(netsurfp[i+2][8]))/2
                zscore = (float(netsurfp[i+1][9]) + float(netsurfp[i+2][9]))/2
                pa = (float(netsurfp[i+1][10]) + float(netsurfp[i+2][10]))/2
                pb = (float(netsurfp[i+1][11]) + float(netsurfp[i+2][11]))/2
                pc = (float(netsurfp[i+1][12]) + float(netsurfp[i+2][12]))/2
                rsaenvt = (float(netsurfp[i-3][7])+float(netsurfp[i-2][7])+float(netsurfp[i-1][7])+float(netsurfp[i][7]) + float(netsurfp[i+3][7])+float(netsurfp[i+4][7])+float(netsurfp[i+5][7])+float(netsurfp[i+6][7]))/8
                zscoreenvt = (float(netsurfp[i-3][9])+float(netsurfp[i-2][9])+float(netsurfp[i-1][9])+float(netsurfp[i][9]) + float(netsurfp[i+3][9])+float(netsurfp[i+4][9])+float(netsurfp[i+5][9])+float(netsurfp[i+6][9]))/8
                results_row = [variants[j][1],variants[j][2],variants[j][3],classassi,rsa,rsaenvt,
                               zscore,zscoreenvt,asa,pa,pb,pc]
                c3.writerow(results_row)
            elif variants[j][3]=='deletion':
                i=i+1
                rsaenvt = (float(netsurfp[i-4][7])+float(netsurfp[i-3][7])+float(netsurfp[i-2][7])+float(netsurfp[i-1][7]) + float(netsurfp[i][7])+float(netsurfp[i+1][7])+float(netsurfp[i+2][7])+float(netsurfp[i+3][7]))/8
                zscoreenvt = (float(netsurfp[i-4][9])+float(netsurfp[i-3][9])+float(netsurfp[i-2][9])+float(netsurfp[i-1][9]) + float(netsurfp[i][9])+float(netsurfp[i+1][9])+float(netsurfp[i+2][9])+float(netsurfp[i+3][9]))/8
                results_row = [variants[j][1],variants[j][2],variants[j][3],'NA','NA',rsaenvt,'NA',zscoreenvt,'NA','NA','NA','NA']
                c3.writerow(results_row)


N.close()
V.close()
R.close()







