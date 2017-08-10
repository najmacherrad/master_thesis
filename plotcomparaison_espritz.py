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
file_wt = 'espritzresults_wt.csv'
file_mut = 'espritzresults_mut.csv'

#--------------------------------------------------
# Compare variant by variant if there is a change after mutation
#--------------------------------------------------
df1 = pd.read_csv(file_wt,'\t')
df2 = pd.read_csv(file_mut,'\t')
df3 = pd.concat([df1['ID'],df1['Gene_name'],df1['Variant_effect'],df1['disorder_prediction'],df1['confidence_label'],
                df2['disorder_prediction'],df2['confidence_label']], axis=1, 
                keys=['ID','Gene_name','Variant_effect','prediction_wt','confidence_label_wt','prediction_mut','confidence_label_mut'])
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
df1List_envt = df1['disorder_envt'].tolist()
df2List_envt = df2['disorder_envt'].tolist()
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
df3.to_csv('espritz_diff_disorder.csv', index=False)

# BARCHART statistique
df = pd.read_csv('espritz_diff_disorder.csv',',')
df = df[df.difference != '0']
dfdiff = df.groupby(['difference']).count()
dfdiff = dfdiff.drop(['Gene_name','Variant_effect','prediction_wt','prediction_mut','confidence_label_wt','confidence_label_mut','confidence_label','difference_envt','envt_score'], 1)
dfdiff.plot(kind='bar',width = 0.2,legend=False,color='k')
plt.ylabel('Number of variants')
plt.savefig('barplotdiffespritz.png')



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
rects2 = ax.bar(ind + width, mut, width, color='r')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('O', 'D'))
ax.set_xlabel('Disorder prediction')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Deleterious DIDA mutants'))
fig.savefig('barplot_espritz_missense.png')

stats.chi2_contingency(np.column_stack((wt,mut))) #(0.055121937273065157, 0.8143790140802849, 1, array([[ 196.5,  196.5],[  44.5,   44.5]]))
