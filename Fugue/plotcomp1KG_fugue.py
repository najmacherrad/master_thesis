# FUGUE
#Compare results between wild type and mutant

# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

def getColumn(filename, column,deli):
    results = csv.reader(open(filename), delimiter=deli)
    return [result[column] for result in results]

#Importer les fichiers
file_wt = 'fugueresultsNEW_wt.csv'
file_mut = 'fugueresultsNEW_1kg.csv'

#--------------------------------------------------
# Compare variant by variant if there is a change after mutation
#--------------------------------------------------
df1 = pd.read_csv(file_wt,'\t')
df2 = pd.read_csv(file_mut,'\t')
df3 = pd.concat([df1['ID'],df1['Gene_name'],df1['Variant_effect'],df1['prediction'],df1['weight'],
                df2['prediction'],df2['weight']], axis=1, 
                keys=['ID','Gene_name','Variant_effect','prediction_wt','weight_wt','prediction_mut','weight_mut'])
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


#for environement
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
df3['weight']=''
df1['weight'] = df1['weight'].fillna(0)
df2['weight'] = df2['weight'].fillna(0)
df3['weight'] = df1['weight'].astype(int).astype(str) + ' to ' + df2['weight'].astype(int).astype(str)
df3['difference_envt']=''
df3['difference_envt'] = newlist_envt
df3['weight_envt']=''
df3['weight_envt']= df1['weight_envt'] + ' to ' + df2['weight_envt']
df3.to_csv('fugue_diff_pred1kg.csv', index=False)

# BARCHART statistique
df = pd.read_csv('fugue_diff_pred1kg.csv',',')
df = df[df.difference != '0']
dfdiff = df.groupby(['difference']).count()
dfdiff.plot(kind='bar',legend=False,color='k')
plt.ylabel('Number of variants')
plt.savefig('barplotdifffugue1kg.png')

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
N = 8
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
wt = (missense_wt.count('A'),missense_wt.count('B'),missense_wt.count('C'),missense_wt.count('G'),missense_wt.count('E'),missense_wt.count('O'),missense_wt.count('P'),missense_wt.count('x'))
rects1 = ax.bar(ind, wt, width, color='grey')
mut = (missense_mut.count('A'),missense_mut.count('B'),missense_mut.count('C'),missense_mut.count('G'),missense_mut.count('E'),missense_mut.count('O'),missense_mut.count('P'),missense_mut.count('x'))
rects2 = ax.bar(ind + width, mut, width, color='blue')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C','G','E','O','P','x'))
ax.set_xlabel('Conformation prediction')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Neutral 1KGP mutants'),loc='upper left')
fig.savefig('barplot_fugue_missense1kg.png')

#CHI SQUARE TEST
stats.chi2_contingency(np.column_stack((wt,mut))) #(68.292174366302532, 3.2683147438853704e-12, 7, array([[ 1013.5,  1013.5],[  418. ,   418. ],[  116.5,   116.5],[   48.5,    48.5],[   61. ,    61. ],[   25. ,    25. ],[   84. ,    84. ],[ 3970.5,  3970.5]]))

#--------------------------------------------------
#COMPARISON
#--------------------------------------------------
file_DIDAmut = 'fugueresults_DIDAmut.csv'

missense_DIDAmut = []
with open(file_DIDAmut, 'rb') as csvfile_mut:
    csvreader_mut = csv.reader(csvfile_mut, delimiter='\t')
    for row in csvreader_mut:
        if row[3]!='prediction':
            missense_DIDAmut.append(row[3])

DIDAmut = (missense_DIDAmut.count('A'),missense_DIDAmut.count('B'),missense_DIDAmut.count('C'),missense_DIDAmut.count('G'),missense_DIDAmut.count('E'),missense_DIDAmut.count('O'),missense_DIDAmut.count('P'),missense_DIDAmut.count('x'))
DIDAmut_freq = [float(x)/241 for x in DIDAmut]


mut = (missense_mut.count('A'),missense_mut.count('B'),missense_mut.count('C'),missense_mut.count('G'),missense_mut.count('E'),missense_mut.count('O'),missense_mut.count('P'),missense_mut.count('x'))
mut_freq = [float(x)/5737 for x in mut]

#Prediction MUTANT
N = 8
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind, DIDAmut_freq, width, color='red')
rects2 = ax.bar(ind + width, mut_freq, width, color='blue')
ax.set_ylabel('Frequency')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C','G','E','O','P','x'))
ax.set_xlabel('Conformation prediction')
ax.legend((rects1[0], rects2[0]), ('deleterious DIDA mutants', 'neutral 1KGP mutants'),loc='upper left')
fig.savefig('barplot_fugue_DIDA1kg_mut.png')

stats.chi2_contingency(np.column_stack((DIDAmut_freq,mut_freq))) #(0.02385403503723927, 0.99999998421700365, 7, array([[  1.81979174e-01,   1.81979174e-01],[  7.78538091e-02,   7.78538091e-02],[  2.31680212e-02,   2.31680212e-02],[  9.97166967e-03,   9.97166967e-03],[  2.44029981e-03,   2.44029981e-03],[  6.10074952e-04,   6.10074952e-04],[  7.66951368e-03,   7.66951368e-03],[  6.96307437e-01,   6.96307437e-01]]))

