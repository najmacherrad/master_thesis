# ELM
#Compare results between wild type and mutant

# coding=utf-8
import numpy as np
import pandas as pd
import csv
from scipy import stats
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

#Importer les fichiers

file_wt = 'elmresults_wt.csv' #531 lines
file_mut = 'elmresults_mut.csv' #453 lines

#WILD TYPE
from collections import defaultdict
results = open('elmresults_wt.csv','r')
r1 = csv.reader(results,delimiter='\t')
per_id = defaultdict(list)
next(r1, None)  # skip the header row
for id,gene,effect,motif,startstop,comm in r1:
    per_id[id.strip()].append(motif.strip()+','+startstop.strip()+','+comm.strip())

csv_file =  open('_elm_grouped_results_wt.csv', 'w') #181 variants
writer = csv.writer(csv_file)
head_row = ['ID','infos']
writer.writerow(head_row)
for key, value in per_id.items():
    join_value=','.join(value)
    writer.writerow([key,join_value])

csv_file.close()

#MUTANT
from collections import defaultdict
results = open('elmresults_mut.csv','r')
r1 = csv.reader(results,delimiter='\t')
per_id = defaultdict(list)
next(r1, None)  # skip the header row
for id,gene,effect,motif,startstop,comm in r1:
    per_id[id.strip()].append(motif.strip()+','+startstop.strip()+','+comm.strip())

csv_file = open('_elm_grouped_results_mut.csv', 'wb') #174 variants
writer = csv.writer(csv_file)
writer.writerow(['ID','infos'])
for key, value in per_id.items():
    join_value=','.join(value)
    writer.writerow([key,join_value])

csv_file.close()

#--------------------------------------------------------------------------------------------
# DICTIONARIES
# Wild type 
results = open(file_wt,'r') 
r1 = csv.reader(results,delimiter='\t')
per_id_wt = defaultdict(list)
next(r1, None)  # skip the header row
for id,gene,effect,motif,startstop,comm in r1:
    per_id_wt[id.strip()].append([gene,effect,motif.strip(),startstop.strip(),comm.strip()])

# Mutant
results = open(file_mut,'r') 
r1 = csv.reader(results,delimiter='\t')
per_id_mut = defaultdict(list)
next(r1, None)  # skip the header row
for id,gene,effect,motif,startstop,comm in r1:
    per_id_mut[id.strip()].append([gene,effect,motif.strip(),startstop.strip(),comm.strip()])

# Statistics
keys_wt = set(per_id_wt.keys()) #181
keys_mut = set(per_id_mut.keys()) #174
intersection = keys_wt & keys_mut #157 qui ont des motifs avant et après mutation


#---------------------------------------------------------------------------------------------
#IDENTICAL
liste_identical=[]
for idwt in per_id_wt:
    for idmut in per_id_mut:
        if idwt==idmut and per_id_wt[idwt]==per_id_mut[idmut]:
            liste_identical.append(idwt)  #40 qui sont exactement identiques


fileI = open('elm_identical.csv','w')
c3 = csv.writer(fileI)
results_row = ['ID','Gene_name','Variant_effect','motif_wt_and_mut','startstop','comm']
c3.writerow(results_row)
for id in liste_identical:
    for value in per_id_mut.get(id):
        results_row = [id,value[0],value[1],value[2],value[3],value[4]]
        c3.writerow(results_row)
fileI.close()

a = pd.read_csv('elm_classes.csv','\t')
b = pd.read_csv('elm_identical.csv',',')
b = b.rename(columns = {'motif_wt_and_mut':'ELMIdentifier'})
merged = b.merge(a, on='ELMIdentifier')
merged.to_csv('_ELM_identical.csv', index=False)

#---------------------------------------------------------------------------------------------
#motif modified (number or type of motif)
modif = intersection - set(liste_identical) #117 -> list of IDs who have modified motifs (number or type of motif)

# Wild type
fileM =open ('elm_modif_wt.csv','w')
c3 = csv.writer(fileM)
results_row = ['ID','Gene_name','Variant_effect','motif_wt','startstop','comm']
c3.writerow(results_row)
for m in modif:
    for value in per_id_wt.get(m):
        results_row = [m,value[0],value[1],value[2],value[3],value[4]]
        c3.writerow(results_row)
fileM.close()

# Mutant
fileM =open ('elm_modif_mut.csv','w')
c3 = csv.writer(fileM)
results_row = ['ID','Gene_name','Variant_effect','motif_mut','startstop','comm']
c3.writerow(results_row)
for m in modif:
    for value in per_id_mut.get(m):
        results_row = [m,value[0],value[1],value[2],value[3],value[4]]
        c3.writerow(results_row)
fileM.close()

# join with ELM CLASSES (257)
a = pd.read_csv('elm_classes.csv','\t')
b = pd.read_csv('elm_modif_wt.csv',',')
b = b.rename(columns = {'motif_wt':'ELMIdentifier'})
bclass = pd.merge(b,a,how='left', on='ELMIdentifier') #391 lines
bclass = bclass.rename(columns = {'ELMIdentifier':'ELMIdentifier_wt'})

c = pd.read_csv('elm_modif_mut.csv',',')
c = c.rename(columns = {'motif_mut':'ELMIdentifier'})
cclass = pd.merge(c,a,how='left', on='ELMIdentifier') #343
cclass = cclass.rename(columns = {'ELMIdentifier':'ELMIdentifier_mut'})

# merge or just concatenate
result = pd.concat([bclass, cclass]) #734 lines
result = result[['ID','Gene_name','Variant_effect','ELMIdentifier_wt','ELMIdentifier_mut','startstop','comm','Accession','FunctionalSiteName','Description','Regex','Probability','#Instances','#Instances_in_PDB']]
result.to_csv('_ELM_modif.csv', index=False)

#---------------------------------------------------------------------------------------------
#LOSS OF FUNCTION
diff1 = keys_wt - intersection #24 variants that are in WT and not in MUTANT

fileM =open ('elm_lossfunction.csv','w')
c3 = csv.writer(fileM)
results_row = ['ID','Gene_name','Variant_effect','motif_wt','startstop','comm']
c3.writerow(results_row)
for d1 in diff1:
    for value in per_id_wt.get(d1):
        results_row = [d1,value[0],value[1],value[2],value[3],value[4]]
        c3.writerow(results_row)
fileM.close()

#GAIN OF FUNCTION
diff2 = keys_mut - intersection #17 variants that are in MUTANT and not in WT

fileM =open ('elm_gainfunction.csv','w')
c3 = csv.writer(fileM)
results_row = ['ID','Gene_name','Variant_effect','motif_mut','startstop','comm']
c3.writerow(results_row)
for d2 in diff2:
    for value in per_id_mut.get(d2):
        results_row = [d2,value[0],value[1],value[2],value[3],value[4]]
        c3.writerow(results_row)
fileM.close()

# join with ELM CLASSES (257)
a = pd.read_csv('elm_classes.csv','\t')
b = pd.read_csv('elm_lossfunction.csv',',')
b = b.rename(columns = {'motif_wt':'ELMIdentifier'})
merged = b.merge(a, on='ELMIdentifier')
merged.to_csv('_ELM_lossfunction.csv', index=False)

a = pd.read_csv('elm_classes.csv','\t')
b = pd.read_csv('elm_gainfunction.csv',',')
b = b.rename(columns = {'motif_mut':'ELMIdentifier'})
merged = b.merge(a, on='ELMIdentifier')
merged.to_csv('_ELM_gainfunction.csv', index=False)


#---------------------------------------------------------------------------------------------
# BARCHART statistics

#----------------------------
# MOTIF PRESENT OR NOT
#----------------------------
N = 2
ind = np.arange(N)
width = 0.30
fig, ax = plt.subplots()
wt = (182,59)
rects1 = ax.bar(ind, wt, width, color='grey')
mut = (175,62)
rects2 = ax.bar(ind + width, mut, width, color='r')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('yes', 'no'))
ax.set_xlabel('ELM motifs')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Deleterious DIDA mutants'))
fig.savefig('barplot_ELM.png')

stats.chi2_contingency(np.column_stack((wt,mut))) #(0.10043241698613907, 0.75131133259493388, 1, array([[ 179.99372385,  177.00627615],[  61.00627615,   59.99372385]]))

#----------------------------
#Type of motifs
#(CLV, DEG, DOC, LIG, MOD, TRG)
#----------------------------
type_WT,type_MUT=[],[]
for kw in keys_wt:
    listWT = per_id_wt[kw]
    for lw in listWT:
        type_WT.append(lw[2].split('_')[0])

for km in keys_mut:
    listMUT = per_id_mut[km]
    for lm in listMUT:
        type_MUT.append(lm[2].split('_')[0])

N = 6
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
wt = (type_WT.count('CLV'),type_WT.count('DEG'),type_WT.count('DOC'),type_WT.count('LIG'),type_WT.count('MOD'),type_WT.count('TRG'))
rects1 = ax.bar(ind, wt, width, color='grey')
mut = (type_MUT.count('CLV'),type_MUT.count('DEG'),type_MUT.count('DOC'),type_MUT.count('LIG'),type_MUT.count('MOD'),type_MUT.count('TRG'))
rects2 = ax.bar(ind + width, mut, width, color='r')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('CLV', 'DEG', 'DOC', 'LIG', 'MOD', 'TRG'))
ax.set_xlabel('Types of ELM motifs')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Deleterious DIDA \n mutants'),loc = 'upper left')
fig.savefig('barplot_ELM_motifs.png')

stats.chi2_contingency(np.column_stack((wt,mut))) #(11.813018301693811, 0.037441657429395876, 5, array([[  31.29878049,   26.70121951],[   9.17378049,    7.82621951],[  61.51829268,   52.48170732],[ 223.94817073,  191.05182927],[ 173.22256098,  147.77743902],[  31.83841463,   27.16158537]]))

#MODIF
dico_all={} #198
for kw in keys_wt:
    type_WT = []
    listWT = per_id_wt[kw]
    for lw in listWT:
        type_WT.append(lw[2].split('_')[0])
    wt = (type_WT.count('CLV'),type_WT.count('DEG'),type_WT.count('DOC'),type_WT.count('LIG'),type_WT.count('MOD'),type_WT.count('TRG'))
    dico_all[kw]=[wt,(0,0,0,0,0,0)]

for km in keys_mut:
    type_MUT = []
    listMUT = per_id_mut[km]
    for lm in listMUT:
        type_MUT.append(lm[2].split('_')[0])
    mut = (type_MUT.count('CLV'),type_MUT.count('DEG'),type_MUT.count('DOC'),type_MUT.count('LIG'),type_MUT.count('MOD'),type_MUT.count('TRG'))
    if km in dico_all.keys():
        dico_all[km][1] = mut
    else:
        dico_all[km] = [(0,0,0,0,0,0),mut]

f = open('_ELM_type_all.csv','w')
f.write('ID\twildtype\tmutant\n')
for key,value in dico_all.iteritems():
    f.write(str(key) + '\t' + str(value[0]) + '\t' +  str(value[1]) + '\n')

f.close()

#----------------------------
#LOSS/GAIN OF FUNCTION
type_loss=[]
for kloss in diff1:
    list_loss = per_id_wt[kloss]
    for lloss in list_loss:
        type_loss.append(lloss[2].split('_')[0])

type_gain=[]
for kgain in diff2:
    list_gain = per_id_mut[kgain]
    for lgain in list_gain:
        type_gain.append(lgain[2].split('_')[0])


N = 6
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
wt = (type_loss.count('CLV'),type_loss.count('DEG'),type_loss.count('DOC'),type_loss.count('LIG'),type_loss.count('MOD'),type_loss.count('TRG'))
rects1 = ax.bar(ind, wt, width, color='grey')
mut = (type_gain.count('CLV'),type_gain.count('DEG'),type_gain.count('DOC'),type_gain.count('LIG'),type_gain.count('MOD'),type_gain.count('TRG'))
rects2 = ax.bar(ind + width, mut, width, color='r')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('CLV', 'DEG', 'DOC', 'LIG', 'MOD', 'TRG'))
ax.set_xlabel('Types of ELM motifs')
plt.ylim(0,18)
ax.legend((rects1[0], rects2[0]), ('Loss of Function', 'Gain of Function'))
fig.savefig('barplot_ELM_loss_gain.png')

stats.chi2_contingency(np.column_stack((wt,mut))) #(15.911395113600994, 0.0071016509069735199, 5, array([[ 11.76923077,   5.23076923],[  1.38461538,   0.61538462],[  4.84615385,   2.15384615],[ 18.        ,   8.        ],[  9.69230769,   4.30769231],[  8.30769231,   3.69230769]]))

#----------------------------
#IDENTICAL
type_identical=[]
for kid in liste_identical:
    list_id = per_id_wt[kid]
    for lid in list_id:
        type_identical.append(lid[2].split('_')[0])

id = (type_identical.count('CLV'),type_identical.count('DEG'),type_identical.count('DOC'),type_identical.count('LIG'),type_identical.count('MOD'),type_identical.count('TRG'))
rects3 = ax.bar(ind + width, id, width, color='grey')

#----------------------------
#MODIF
dico_modif={}
for kmod in modif:
    type_mod_wt = []
    type_mod_mut = []
    list_mod_wt = per_id_wt[kmod]
    list_mod_mut = per_id_mut[kmod]
    for lmwt in list_mod_wt:
        type_mod_wt.append(lmwt[2].split('_')[0])
    for lmmut in list_mod_mut:
        type_mod_mut.append(lmmut[2].split('_')[0])
    nb_type_wt = (type_mod_wt.count('CLV'),type_mod_wt.count('DEG'),type_mod_wt.count('DOC'),type_mod_wt.count('LIG'),type_mod_wt.count('MOD'),type_mod_wt.count('TRG'))
    nb_type_mut = (type_mod_mut.count('CLV'),type_mod_mut.count('DEG'),type_mod_mut.count('DOC'),type_mod_mut.count('LIG'),type_mod_mut.count('MOD'),type_mod_mut.count('TRG'))
    dico_modif[kmod]=[nb_type_wt,nb_type_mut]

#CLV DEG DOC LIG MOD TRG
f = open('_ELM_type_modif.csv','w')
f.write('ID\twildtype\tmutant\n')
for key,value in dico_modif.iteritems():
    f.write(str(key) + '\t' + str(value[0]) + '\t' +  str(value[1]) + '\n')
f.close()

#----------------------------#----------------------------
#----------------------------#----------------------------
#DISEASE TABLE
#----------------------------#----------------------------
#----------------------------#----------------------------
a = pd.read_csv('_ELM_type_all.csv','\t')
b = pd.read_csv('didavariantskey.csv','\t')
c = a.merge(b,on='ID')
c.to_csv('_ELM_table.csv', index=False)

