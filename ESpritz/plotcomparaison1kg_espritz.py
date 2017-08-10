#Espritz
#Compare results between wild type and mutant

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
                
#Import files
file_wt = 'espritzresultsNEW_wt.csv'
file_mut = 'espritzresultsNEW_1kg.csv'

#--------------------------------------------------
# Compare variant by variant if there is a change after mutation
#--------------------------------------------------
df1 = pd.read_csv(file_wt,'\t')
df2 = pd.read_csv(file_mut,'\t')
df3 = pd.concat([df1['ID'],df1['Gene_name'],df1['Variant_effect'],df1['prediction'],df1['confidence_label'],
                df2['prediction'],df2['confidence_label']], axis=1,
                keys=['ID','Gene_name','Variant_effect','prediction_wt','score_wt','prediction_mut','score_mut'])
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

#environnement
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
df3['confidence_label']=''
df1['confidence_label'] = df1['confidence_label'].fillna(0)
df2['confidence_label'] = df2['confidence_label'].fillna(0)
df3['confidence_label'] = df1['confidence_label'].astype(int).astype(str) + ' to ' + df2['confidence_label'].astype(int).astype(str)
df3['difference_envt']=''
df3['difference_envt'] = newlist_envt
df3['envt_score']=''
df3['envt_score']= df1['confidence_label_envt'] + ' to ' + df2['confidence_label_envt']
df3.to_csv('espritz_diff_disorder1kg.csv', index=False)

# BARCHART statistique
df = pd.read_csv('espritz_diff_disorder1kg.csv',',')
df = df[df.difference != '0']
dfdiff = df.groupby(['difference']).count()
dfdiff = dfdiff.drop(['Gene_name','Variant_effect','prediction_wt','prediction_mut','score_wt','score_mut','confidence_label','difference_envt','envt_score'], 1)
dfdiff.plot(kind='bar',legend=False,color='k',width = 0.2)
plt.ylabel('Number of variants')
plt.savefig('barplotdiffespritz1kg.png')

#--------------------------------------------------
#Prediction
pred_wt = getColumn(file_wt,3,'\t')
pred_mut = getColumn(file_mut,3,'\t')
pred_wt.pop(0)
pred_mut.pop(0)

missense_wt,missense_mut=[],[]
for i in range(0,len(pred_wt)):
    if pred_wt[i]=='NA':
        missense_wt.append(np.nan)
    else:
        missense_wt.append(pred_wt[i])
for i in range(0,len(pred_mut)):
    if pred_mut[i]=='NA':
        missense_mut.append(np.nan)
    else:
        missense_mut.append(pred_mut[i])

N = 2
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
wt = (missense_wt.count('O'),missense_wt.count('D'))
rects1 = ax.bar(ind, wt, width, color='grey')
mut = (missense_mut.count('O'),missense_mut.count('D'))
rects2 = ax.bar(ind + width, mut, width, color='b')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('O', 'D'))
ax.set_xlabel('Disorder prediction')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Neutral 1KGP mutants'))
fig.savefig('barplot_espritz1kg_missense.png')

#--------------------------------------------------
# COMPARISON DIDA MUT VS 1KGP WT
#--------------------------------------------------
file_DIDAmut = 'espritzresults_mut.csv'

missense_DIDAmut = []
with open(file_DIDAmut, 'rb') as csvfile_mut:
    csvreader_mut = csv.reader(csvfile_mut, delimiter='\t')
    for row in csvreader_mut:
        if row[3]!='disorder_prediction':
            missense_DIDAmut.append(row[3])

DIDAmut = (missense_DIDAmut.count('O'),missense_DIDAmut.count('D'))
DIDAmut_freq = [float(x)/241 for x in DIDAmut]

mut = (missense_mut.count('O'),missense_mut.count('D'))
mut_freq = [float(x)/5512 for x in mut]

#Prediction MUTANT
N = 2
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind, DIDAmut_freq, width, color='red')
rects2 = ax.bar(ind + width, mut_freq, width, color='blue')
ax.set_ylabel('Frequency')
ax.set_xticks(ind + width)
ax.set_xticklabels(('O', 'D'))
ax.set_xlabel('Disorder prediction')
ax.legend((rects1[0], rects2[0]), ('Deleterious DIDA mutants', 'Neutral 1KGP mutants'))
fig.savefig('barplot_espritz1kg_DIDA1kg_mut.png')

stats.chi2_contingency(np.column_stack((DIDAmut_freq,mut_freq)))# ((2.8991612092878767, 0.088625658868333151, 1, array([[ 0.80166284,  0.80166284],[ 0.19833716,  0.19833716]]))

