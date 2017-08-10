# PRELUDE
# Compare results between wild type and mutant

# coding=utf-8
import numpy as np
import pandas as pd
import csv
from scipy import stats
import matplotlib.pyplot as plt

def getColumn(filename, column,deli):
    results = csv.reader(open(filename), delimiter=deli)
    return [result[column] for result in results]

#Importer les fichiers
file_wt = 'preluderesults_wt.csv'
file_mut = 'preluderesults_mut.csv'

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
df3.to_csv('prelude_diff_pred.csv', index=False)

# BARCHART statistique
df = pd.read_csv('prelude_diff_pred.csv',',')
df = df[df.difference != '0']
dfdiff = df.groupby(['difference']).count()
dfdiff = dfdiff.drop(['Gene_name','Variant_effect','prediction_wt','prediction_mut','difference_envt'], 1)
dfdiff.plot(kind='bar',legend=False,color='k')
plt.ylabel('Number of variants')
plt.savefig('barplotdiffprelude.png')


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
rects2 = ax.bar(ind + width, mut, width, color='r')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C','G','E','O','P'))
ax.set_xlabel('Secondary structure prediction')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Deleterious DIDA mutants'))
fig.savefig('barplot_prelude_missense.png')

#CHI SQUARE TEST
stats.chi2_contingency(np.column_stack((wt,mut))) #(2.4873908105954228, 0.86987615706529087, 6, array([[ 88.5,  88.5],[ 55.5,  55.5],[ 48.5,  48.5],[ 10.5,  10.5],[ 13.5,  13.5],[  4.5,   4.5],[ 20. ,  20. ]]))

