# Pop Hot Music

# coding=utf-8
from numpy import *
import petl as etl
from re import *
import operator
import glob 
import pandas as pd
import re
import csv


dico_aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
#-----------------
# POPMusic
#-----------------
#recuperer tous les fichier .pop
#Begin at line 14

popmusic=open('popmusic.csv','w')
popmusic.write('PDB_ID,Type,Chain,aa_nb,aa_wt,aa_mut,Secondary_structure,Solvent_accessibility,DDG\n')
#PopHotMusic_pdb
for fil in glob.glob('PopHotMusic_pdb/*.pop'):
    name = str.split(fil,'.')[0].split('/')[1]
    f=open(fil, 'r')
    lines = f.readlines()
    for i in range(14,len(lines)):
        a = lines[i][:1] + ' ' + lines[i][1:]
        l=a.split()
        l[2]= dico_aa[l[2]]
        l[3]= dico_aa[l[3]]
        newline = ",".join(l)
        popmusic.write(name + ',PDB_code,' + newline + '\n')
    f.close()
#PopHotMusic_experimental OK
for fil in glob.glob('PopHotMusic_experimental/*.pop'):
    name = str.split(fil,'.')[0].split('/')[1].replace('pdb','')
    f=open(fil, 'r')
    lines = f.readlines()
    for i in range(14,len(lines)):
        l = lines[i].split()
        l[2]= dico_aa[l[2]]
        l[3]= dico_aa[l[3]]
        newline = ",".join(l)
        popmusic.write(name + ',Experimental,' + newline + '\n')
    f.close()
#PopHotMusic_modbase OK
for fil in glob.glob('PopHotMusic_modbase/*.pop'):
    name = str.split(fil,'.')[0].split('/')[1]
    f=open(fil, 'r')
    lines = f.readlines()
    for i in range(14,len(lines)):
        l = lines[i].split()
        l[1]= dico_aa[l[1]]
        l[2]= dico_aa[l[2]]
        newline = ",".join(l)
        popmusic.write(name + ',Modbase_model,NA,' + newline + '\n')
    f.close()
popmusic.close()

# Create a new table with only the variants
pop=open('popmusic.csv','r')
struct_variants=open('structural_DIDAvariants_hotpop.csv','r')
R=open ('popmusic_results.csv','w')
c1 = csv.reader(pop)
c2 = csv.reader(struct_variants,delimiter='\t')
c3 = csv.writer(R)
P = list(c1)
results_row = ['VariantID','UniprotID','PDBID','Type','Chain','Secondary_structure','Solvent_accessibility','DDG']
c3.writerow(results_row)
iterc2 = iter(c2)
next(iterc2)
for variants in iterc2:
    #PDB_code
    if len(variants[3])==4:
        PDB = variants[3].lower()
    #Experimental,Theoritical_model_1FV4 and Modbase
    else:
        PDB = variants[3].split('_')[-1]
    template = int(variants[5].split('-')[0])
    target = int(variants[6].split('-')[0])
    shift = target - template
    aa_wt = (variants[2].split('(')[1].split(')')[0])[0]
    aa_mut = (variants[2].split('(')[1].split(')')[0])[-1:]
    objet = re.search(r'[0-9]+',variants[2])
    positionPROT = int(objet.group())
    positionPDB = int(objet.group()) - shift   #position - shift
    #FOR PDB_code and Experimental types
    if variants[3].split('_')[0]!='Modbase':
        for row in P:
            if row[0] == PDB and row[3]==str(positionPDB) and row[4]==aa_wt and row[5]==aa_mut:
                results_row= [variants[0],variants[1],row[0],row[1],row[2],row[6],row[7],row[8]]
                c3.writerow(results_row)
    else:
        if variants[1] == 'P14222':
            key = 'P14222_6variants'
        elif variants[1] == 'Q07837':
            key = 'Q07837_variant56_59'
        else:
            key = variants[1] + 'variant' + variants[0]
        for row in P:
            if row[0]==key and row[3]==str(positionPROT) and row[4]==aa_wt and row[5]==aa_mut:
                results_row= [variants[0],variants[1],PDB,row[1],row[2],row[6],row[7],row[8]]
                c3.writerow(results_row)
pop.close()
struct_variants.close()
R.close()


#-----------------
# HOTMusic
#-----------------
#recuperer tous les fichier .hot
#Begin at line 17
hotmusic=open('hotmusic.csv','w')
hotmusic.write('PDB_ID,Type,Chain,aa_nb,aa_wt,aa_mut,Secondary_structure,Solvent_accessibility,DTm\n')
#PopHotMusic_pdb
for fil in glob.glob('PopHotMusic_pdb/*.hot'):
    name = str.split(fil,'.')[0].split('/')[1]
    f=open(fil, 'r')
    lines = f.readlines()
    for i in range(17,len(lines)):
        a = lines[i][:1] + ' ' + lines[i][1:]
        l=a.split()
        l[2]= dico_aa[l[2]]
        l[3]= dico_aa[l[3]]
        newline = ",".join(l)
        hotmusic.write(name + ',PDB_code,' + newline + '\n')
    f.close()
#PopHotMusic_experimental OK
for fil in glob.glob('PopHotMusic_experimental/*.hot'):
    name = str.split(fil,'.')[0].split('/')[1].replace('pdb','')
    f=open(fil, 'r')
    lines = f.readlines()
    for i in range(17,len(lines)):
        l = lines[i].split()
        l[2]= dico_aa[l[2]]
        l[3]= dico_aa[l[3]]
        newline = ",".join(l)
        hotmusic.write(name + ',Experimental,' + newline + '\n')
    f.close()
#PopHotMusic_modbase OK
for fil in glob.glob('PopHotMusic_modbase/*.hot'):
    name = str.split(fil,'.')[0].split('/')[1]
    f=open(fil, 'r')
    lines = f.readlines()
    for i in range(17,len(lines)):
        l = lines[i].split()
        l[1]= dico_aa[l[1]]
        l[2]= dico_aa[l[2]]
        newline = ",".join(l)
        hotmusic.write(name + ',Modbase_model,NA,' + newline + '\n')
    f.close()
hotmusic.close()

# Create a new table with only the variants
hot=open('hotmusic.csv','r')
struct_variants=open('structural_DIDAvariants_hotpop.csv','r')
R=open ('hotmusic_results.csv','w')
c1 = csv.reader(hot)
c2 = csv.reader(struct_variants,delimiter='\t')
c3 = csv.writer(R)
H = list(c1)
results_row = ['VariantID','UniprotID','PDBID','Type','Chain','Secondary_structure','Solvent_accessibility','DTm']
c3.writerow(results_row)
iterc2 = iter(c2)
next(iterc2)
for variants in iterc2:
    #PDB_code
    if len(variants[3])==4:
        PDB = variants[3].lower()
    #Experimental,Theoritical_model_1FV4 and Modbase
    else:
        PDB = variants[3].split('_')[-1]
    template = int(variants[5].split('-')[0])
    target = int(variants[6].split('-')[0])
    shift = target - template
    aa_wt = (variants[2].split('(')[1].split(')')[0])[0]
    aa_mut = (variants[2].split('(')[1].split(')')[0])[-1:]
    objet = re.search(r'[0-9]+',variants[2])
    positionPROT = int(objet.group())
    positionPDB = int(objet.group()) - shift   #position - shift
    #FOR PDB_code and Experimental types
    if variants[3].split('_')[0]!='Modbase':
        for row in H:
            if row[0] == PDB and row[3]==str(positionPDB) and row[4]==aa_wt and row[5]==aa_mut:
                results_row= [variants[0],variants[1],row[0],row[1],row[2],row[6],row[7],row[8]]
                c3.writerow(results_row)
    else:
        if variants[1] == 'P14222':
            key = 'P14222_6variants'
        elif variants[1] == 'Q07837':
            key = 'Q07837_variant56_59'
        else:
            key = variants[1] + 'variant' + variants[0]
        for row in H:
            if row[0]==key and row[3]==str(positionPROT) and row[4]==aa_wt and row[5]==aa_mut:
                results_row= [variants[0],variants[1],PDB,row[1],row[2],row[6],row[7],row[8]]
                c3.writerow(results_row)
hot.close()
struct_variants.close()
R.close()

