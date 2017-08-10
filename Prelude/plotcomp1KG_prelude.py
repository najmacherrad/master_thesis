# PRELUDE
# Compare results between wild type and mutant

# coding=utf-8
import numpy as np
import pandas as pd
import csv
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

def getColumn(filename, column,deli):
    results = csv.reader(open(filename), delimiter=deli)
    return [result[column] for result in results]

#Importer les fichiers
file_wt = 'preluderesultsNEW_wt.csv'
file_mut = 'preluderesultsNEW_1kg.csv'

#--------------------------------------------------
# Compare variant by variant if there is a change after mutation
#--------------------------------------------------
df1 = pd.read_csv(file_wt,'\t')
df2 = pd.read_csv(file_mut,'\t')
df3 = pd.concat([df1['ID'],df1['Gene_name'],df1['Variant_effect'],df1['prediction'],
                df2['prediction']], axis=1, 
                keys=['ID','Gene_name','Variant_effect','prediction_wt','prediction_mut'])
df3 = df3.where((pd.notnull(df3)), None)

df1List = df3['prediction_wt'].tolist()
df2List = df3['prediction_mut'].tolist()

newlist=[]
for i in range(0,len(df1List)):
    for j in range(0,len(df2List)):
        if i==j:
            if df1List[i]!=df2List[i]:
                newlist.append(str(df1List[i])+' to '+ str(df2List[i]))
            else:
                newlist.append('0')

# For environment prediction
df1List_envt = df1['prediction_envt'].tolist()
df2List_envt = df2['prediction_envt'].tolist()
newlist_envt=[]
for i in range(0,len(df1List_envt)):
    for j in range(0,len(df2List_envt)):
        if i==j:
            if df1List_envt[i]!=df2List_envt[i]:
                newlist_envt.append(str(df1List_envt[i])+' to '+ str(df2List_envt[i]))
            else:
                newlist_envt.append('0')

df3['difference']=''  
df3['difference'] = newlist
df3['difference_envt']=''
df3['difference_envt'] = newlist_envt
df3.to_csv('prelude_diff_pred1kg.csv', index=False)

# BARCHART statistique
df = pd.read_csv('prelude_diff_pred1kg.csv',',')
df = df[df.difference != '0']
dfdiff = df.groupby(['difference']).count()
dfdiff = dfdiff.drop(['Gene_name','Variant_effect','prediction_wt','prediction_mut','difference_envt'], 1)
dfdiff = dfdiff.drop(dfdiff.index[[28,29,30,31,32,33,34]])
dfdiff.plot(kind='bar',legend=False,color='k')
plt.ylabel('Number of variants')
plt.savefig('barplotdiffprelude1kg.png')


#--------------------------------------------------
#Prediction
struct_wt = getColumn(file_wt,3,'\t')
struct_mut = getColumn(file_mut,3,'\t')
struct_wt.pop(0)
struct_mut.pop(0)

missense_wt,missense_mut=[],[]
for i in range(0,len(struct_wt)):
    if struct_wt[i]=='NA':
        missense_wt.append(np.nan)
    else:
        missense_wt.append(struct_wt[i])
for i in range(0,len(struct_mut)):
    if struct_mut[i]=='NA':
        missense_mut.append(np.nan)
    else:
        missense_mut.append(struct_mut[i])


#Prediction
N = 7
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
wt = (missense_wt.count('A'),missense_wt.count('B'),missense_wt.count('C'),missense_wt.count('G'),missense_wt.count('E'),missense_wt.count('O'),missense_wt.count('P'))
rects1 = ax.bar(ind, wt, width, color='grey')
mut = (missense_mut.count('A'),missense_mut.count('B'),missense_mut.count('C'),missense_mut.count('G'),missense_mut.count('E'),missense_mut.count('O'),missense_mut.count('P'))
rects2 = ax.bar(ind + width, mut, width, color='blue')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C','G','E','O','P'))
ax.set_xlabel('Secondary structure prediction')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Neutral 1KGP mutants'))
fig.savefig('barplot_prelude_missense1kg.png')

#CHI SQUARE TEST
stats.chi2_contingency(np.column_stack((wt,mut))) #(209.94466575946325, 1.4470673433772526e-42, 6, array([[ 1789.10974873,  1826.89025127],[ 1456.6203264 ,  1487.3796736 ],[ 1173.11372075,  1197.88627925],[  289.44391676,   295.55608324],[  239.96632415,   245.03367585],[  118.25144633,   120.74855367],[  663.49451688,   677.50548312]]))

#--------------------------------------------------------------
#COMPARISON DELETERIOUS DIDA MUTANTS AND NEUTRAL 1KGP MUTANTS
#--------------------------------------------------------------
file_DIDAmut = 'preluderesults_DIDAmut.csv'

missense_DIDAmut = []
with open(file_DIDAmut, 'rb') as csvfile_mut:
    csvreader_mut = csv.reader(csvfile_mut, delimiter='\t')
    for row in csvreader_mut:
        if row[3]!='prediction':
            missense_DIDAmut.append(row[3])

DIDAmut = (missense_DIDAmut.count('A'),missense_DIDAmut.count('B'),missense_DIDAmut.count('C'),missense_DIDAmut.count('G'),missense_DIDAmut.count('E'),missense_DIDAmut.count('O'),missense_DIDAmut.count('P'))
DIDAmut_freq = [float(x)/241 for x in DIDAmut]

mut = (missense_mut.count('A'),missense_mut.count('B'),missense_mut.count('C'),missense_mut.count('G'),missense_mut.count('E'),missense_mut.count('O'),missense_mut.count('P'))
mut_freq = [float(x)/5851 for x in mut]

#Prediction MUTANT
N = 7
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind, DIDAmut_freq, width, color='red')
rects2 = ax.bar(ind + width, mut_freq, width, color='blue')
ax.set_ylabel('Frequency')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C','G','E','O','P'))
ax.set_xlabel('Secondary structure prediction')
ax.legend((rects1[0], rects2[0]), ('deleterious DIDA mutants', 'neutral 1KGP mutants'))
fig.savefig('barplot_prelude_DIDA1kg_mut.png')

stats.chi2_contingency(np.column_stack((DIDAmut_freq,mut_freq)))
# 0.022472769003365575, 0.99999976553953085, 6 , array([[ 0.33768539,  0.3447742 ],[ 0.25626832,  0.261648  ],[ 0.2213582 ,  0.22600503],[ 0.04139867,  0.04226773],[ 0.03528948,  0.03603029],[ 0.01259061,  0.01285492],[ 0.09540933,  0.09741219]])


