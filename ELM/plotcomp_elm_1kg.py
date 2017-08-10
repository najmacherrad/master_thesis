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

#Import files
file_wt = 'elmresultsNEW_wt.csv' #3757 variants (11'039 lines)
file_mut = 'elmresultsNEW_1kg.csv' #3744 (10'818 lines)

#Wild type
from collections import defaultdict
results = open('elmresultsNEW_wt.csv','r')
r1 = csv.reader(results,delimiter='\t')
per_id = defaultdict(list)
next(r1, None)  # skip the header row
for id,gene,effect,motif,startstop,comm in r1:
    per_id[id.strip()].append(motif.strip()+','+startstop.strip()+','+comm.strip())

csv_file =  open('_elm_grouped_results_wt.csv', 'w')
writer = csv.writer(csv_file)
head_row = ['ID','infos']
writer.writerow(head_row)
for key, value in per_id.items():
    join_value=','.join(value)
    writer.writerow([key,join_value])

csv_file.close()


#Mutant
from collections import defaultdict
results = open('elmresultsNEW_1kg.csv','r')
r1 = csv.reader(results,delimiter='\t')
per_id = defaultdict(list)
next(r1, None)  # skip the header row
for id,gene,effect,motif,startstop,comm in r1:
    per_id[id.strip()].append(motif.strip()+','+startstop.strip()+','+comm.strip())

csv_file = open('_elm_grouped_results_1kg.csv', 'wb')
writer = csv.writer(csv_file)
writer.writerow(['ID','infos'])
for key, value in per_id.items():
    join_value=','.join(value)
    writer.writerow([key,join_value])

csv_file.close()

#--------------------------------------------------------------------------------------------
# DICTIONARY
# Wild types
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
keys_wt = set(per_id_wt.keys()) #3'757
keys_mut = set(per_id_mut.keys()) #3'744
intersection = keys_wt & keys_mut #3'312 with motifs before and after mutation


#------------------------------------------------------------------------------------------------
#IDENTICAL
liste_identical=[]
for idwt in per_id_wt:
    for idmut in per_id_mut:
        if idwt==idmut and per_id_wt[idwt]==per_id_mut[idmut]:
            liste_identical.append(idwt)  #837 exactely identical


fileI = open('elm_identical_1kg.csv','w')
c3 = csv.writer(fileI)
results_row = ['ID','Gene_name','Variant_effect','motif_wt_and_mut','startstop','comm']
c3.writerow(results_row)
for id in liste_identical:
    for value in per_id_mut.get(id):
        results_row = [id,value[0],value[1],value[2],value[3],value[4]]
        c3.writerow(results_row)

fileI.close()

a = pd.read_csv('elm_classes.csv','\t')
b = pd.read_csv('elm_identical_1kg.csv',',')
b = b.rename(columns = {'motif_wt_and_mut':'ELMIdentifier'})
merged = b.merge(a, on='ELMIdentifier')
merged.to_csv('_ELM_identical_1kg.csv', index=False)

#------------------------------------------------------------------------------------------------
#motif modified (number or type)
modif = intersection - set(liste_identical) #2475 list of IDs which have motif modified (number or type)

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
fileM =open ('elm_modif_1kg.csv','w')
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
bclass = pd.merge(b,a,how='left', on='ELMIdentifier')
bclass = bclass.rename(columns = {'ELMIdentifier':'ELMIdentifier_wt'})

c = pd.read_csv('elm_modif_1kg.csv',',')
c = c.rename(columns = {'motif_mut':'ELMIdentifier'})
cclass = pd.merge(c,a,how='left', on='ELMIdentifier')
cclass = cclass.rename(columns = {'ELMIdentifier':'ELMIdentifier_mut'})

# merge or just concatenate
result = pd.concat([bclass, cclass])
result = result[['ID','Gene_name','Variant_effect','ELMIdentifier_wt','ELMIdentifier_mut','startstop','comm','Accession','FunctionalSiteName','Description','Regex','Probability','#Instances','#Instances_in_PDB']]
result.to_csv('_ELM_modif_1kg.csv', index=False)

#------------------------------------------------------------------------------------------------
#LOSS OF FUNCTION
diff1 = keys_wt - intersection #445 variants that are in WT and not in MUTANT

fileM =open ('elm_lossfunction_1kg.csv','w')
c3 = csv.writer(fileM)
results_row = ['ID','Gene_name','Variant_effect','motif_wt','startstop','comm']
c3.writerow(results_row)
for d1 in diff1:
    for value in per_id_wt.get(d1):
        results_row = [d1,value[0],value[1],value[2],value[3],value[4]]
        c3.writerow(results_row)

fileM.close()

#GAIN OF FUNCTION
diff2 = keys_mut - intersection #432 variants that are in MUTANT and not in WT

fileM =open ('elm_gainfunction_1kg.csv','w')
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
b = pd.read_csv('elm_lossfunction_1kg.csv',',')
b = b.rename(columns = {'motif_wt':'ELMIdentifier'})
merged = b.merge(a, on='ELMIdentifier')
merged.to_csv('_ELM_lossfunction_1kg.csv', index=False)

a = pd.read_csv('elm_classes.csv','\t')
b = pd.read_csv('elm_gainfunction_1kg.csv',',')
b = b.rename(columns = {'motif_mut':'ELMIdentifier'})
merged = b.merge(a, on='ELMIdentifier')
merged.to_csv('_ELM_gainfunction_1kg.csv', index=False)


#------------------------------------------------------------------------------------------------
# BARCHART statistics

#----------------------------
# MOTIF PRESENT OR NOT
#----------------------------
N = 2
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
wt = (3757,2094)
rects1 = ax.bar(ind, wt, width, color='grey')
mut = (3744,2107)
rects2 = ax.bar(ind + width, mut, width, color='blue')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('yes', 'no'))
ax.set_xlabel('ELM motifs')
ax.legend((rects1[0], rects2[0]), ('Wild type residues', 'Neutral 1KGP mutants'))
fig.savefig('barplot_ELM_1kg.png')

stats.chi2_contingency(np.column_stack((wt,mut))) #(0.053474993304867921, 0.81712322382190838, 1, array([[ 3750.5,  3750.5],[ 2100.5,  2100.5]]))

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
rects2 = ax.bar(ind + width, mut, width, color='blue')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('CLV', 'DEG', 'DOC', 'LIG', 'MOD', 'TRG'))
ax.set_xlabel('Types of ELM motifs')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Neutral 1KGP mutants'),loc = 'upper left')
fig.savefig('barplot_ELM_motifs_1kg.png')

stats.chi2_contingency(np.column_stack((wt,mut))) #(88.653532685212184, 1.2885598135068658e-17, 5, array([[  730.31038111,   715.68961889],[  138.89028686,   136.10971314],[ 1295.97264034,  1270.02735966],[ 4432.36784554,  4343.63215446],[ 3802.56352656,  3726.43647344],[  638.89531958,   626.10468042]]))

#MODIF
dico_all={} #4189
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

f = open('_ELM_type_all_1kg.csv','w')
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
rects2 = ax.bar(ind + width, mut, width, color='blue')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('CLV', 'DEG', 'DOC', 'LIG', 'MOD', 'TRG'))
ax.set_xlabel('Types of ELM motifs')
ax.legend((rects1[0], rects2[0]), ('Loss of Function', 'Gain of Function'))
fig.savefig('barplot_ELM_loss_gain_1kg.png')

stats.chi2_contingency(np.column_stack((wt,mut))) #(70.075128548637792, 9.8855399805002182e-14, 5, array([[ 109.30238241,  101.69761759],[   9.84239462,    9.15760538],[ 105.15821625,   97.84178375],[ 340.33964569,  316.66035431],[ 212.90653635,  198.09346365],[  70.45082468,   65.54917532]]))

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
f = open('_ELM_type_modif_1kg.csv','w')
f.write('ID\twildtype\tmutant\n')
for key,value in dico_modif.iteritems():
    f.write(str(key) + '\t' + str(value[0]) + '\t' +  str(value[1]) + '\n')
f.close()

#--------------------------------------------------
# COMPARISON DIDA MUT VS 1KGP WT
#--------------------------------------------------
file_DIDAwt = 'elmresults_DIDAwt.csv'
file_DIDAmut = 'elmresults_DIDAmut.csv'

#DIDA
#wildtype
results = open(file_DIDAwt,'r')
r1 = csv.reader(results,delimiter='\t')
per_id_wt1 = defaultdict(list)
next(r1, None)  # skip the header row
for id,gene,effect,motif,startstop,comm in r1:
    per_id_wt1[id.strip()].append([gene,effect,motif.strip(),startstop.strip(),comm.strip()])

# Mutant
results = open(file_DIDAmut,'r')
r1 = csv.reader(results,delimiter='\t')
per_id_mut1 = defaultdict(list)
next(r1, None)  # skip the header row
for id,gene,effect,motif,startstop,comm in r1:
    per_id_mut1[id.strip()].append([gene,effect,motif.strip(),startstop.strip(),comm.strip()])

keys_wt1 = set(per_id_wt1.keys()) #181
keys_mut1 = set(per_id_mut1.keys())

type_WT1,type_MUT1=[],[]
for kw in keys_wt1:
    listWT1 = per_id_wt1[kw]
    for lw in listWT1:
        type_WT1.append(lw[2].split('_')[0])

for km in keys_mut1:
    listMUT1 = per_id_mut1[km]
    for lm in listMUT1:
        type_MUT1.append(lm[2].split('_')[0])


DIDAwt = (type_WT1.count('CLV'),type_WT1.count('DEG'),type_WT1.count('DOC'),type_WT1.count('LIG'),type_WT1.count('MOD'),type_WT1.count('TRG'))
DIDAwt_freq = [float(x)/237.0 for x in DIDAwt]

DIDAmut = (type_MUT1.count('CLV'),type_MUT1.count('DEG'),type_MUT1.count('DOC'),type_MUT1.count('LIG'),type_MUT1.count('MOD'),type_MUT1.count('TRG'))
DIDAmut_freq = [float(x)/237.0 for x in DIDAmut]

#1KGP
type_WT,type_MUT=[],[]
for kw in keys_wt:
    listWT = per_id_wt[kw]
    for lw in listWT:
        type_WT.append(lw[2].split('_')[0])

for km in keys_mut:
    listMUT = per_id_mut[km]
    for lm in listMUT:
        type_MUT.append(lm[2].split('_')[0])

wt = (type_WT.count('CLV'),type_WT.count('DEG'),type_WT.count('DOC'),type_WT.count('LIG'),type_WT.count('MOD'),type_WT.count('TRG'))
wt_freq = [float(x)/5851.0 for x in wt]

mut = (type_MUT.count('CLV'),type_MUT.count('DEG'),type_MUT.count('DOC'),type_MUT.count('LIG'),type_MUT.count('MOD'),type_MUT.count('TRG'))
mut_freq = [float(x)/5851.0 for x in mut]

#Prediction WILD TYPE
N = 6
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind, DIDAwt_freq, width, color='red')
rects2 = ax.bar(ind + width, wt_freq, width, color='blue')
ax.set_ylabel('Frequency')
ax.set_xticks(ind + width)
ax.set_xticklabels(('CLV', 'DEG', 'DOC', 'LIG', 'MOD', 'TRG'))
ax.set_xlabel('Types of ELM motifs')
ax.legend((rects1[0], rects2[0]), ('Wild types (DIDA)', 'Wild types (1KGP)'),loc = 'upper left')
fig.savefig('barplot_ELM_motifs_DIDA1kg_wt.png')

stats.chi2_contingency(np.column_stack((DIDAwt_freq,wt_freq)))
    # (0.004166091513110399, 0.99999994049894358, 5, array([[ 0.17407648,  0.14658636],[ 0.03608065,  0.0303828 ],[ 0.24484019,  0.20617508],[ 0.89962621,  0.75755742],[ 0.72984077,  0.61458447],[ 0.15604203,  0.13139991]]))

#Prediction MUTANT
N = 6
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind, DIDAmut_freq, width, color='red')
rects2 = ax.bar(ind + width, mut_freq, width, color='blue')
ax.set_ylabel('Frequency')
ax.set_xticks(ind + width)
ax.set_xticklabels(('CLV', 'DEG', 'DOC', 'LIG', 'MOD', 'TRG'))
ax.set_xlabel('Types of ELM motifs')
ax.legend((rects1[0], rects2[0]), ('Deleterious DIDA mutants', 'Neutral 1KGP mutants'),loc = 'upper left')
fig.savefig('barplot_ELM_motifs_DIDA1kg_mut.png')

stats.chi2_contingency(np.column_stack((DIDAmut_freq,mut_freq)))
    # (0.011563499081696788, 0.99999923830913207, 5, array([[ 0.08702232,  0.08417782],[ 0.02656772,  0.0256993 ],[ 0.23816999,  0.23038492],[ 0.8101333 ,  0.78365247],[ 0.65917001,  0.63762372],[ 0.09032908,  0.08737649]]))
