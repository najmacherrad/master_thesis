# NetSurfP
#Compare results between wild type and mutant

# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from scipy import stats
from pylab import plot, show, savefig, xlim, figure, \
                hold, ylim, legend, boxplot, setp, axes

def getColumn(filename, column,deli):
    results = csv.reader(open(filename), delimiter=deli)
    return [result[column] for result in results]

#Import files
file_wt = 'netsurfpresults_wt.csv'
file_mut = 'netsurfpresults_mut.csv'

#-----------------------------------------------------------------------------
# RSA
#-----------------------------------------------------------------------------
#----------------
# SCATTER PLOT
RSA_wt = getColumn(file_wt,4,'\t')
RSA_mut = getColumn(file_mut,4,'\t')
RSA_wt.pop(0)
RSA_mut.pop(0)

x,y=[],[]
for i in range(0,len(RSA_wt)):
    if RSA_wt[i]=='NA':
        x.append(np.nan)
    else:
        x.append(float(RSA_wt[i]))
for i in range(0,len(RSA_mut)):
    if RSA_mut[i]=='NA':
        y.append(np.nan)
    else:
        y.append(float(RSA_mut[i]))

             
fig = plt.figure()
a=b=[0,0.2,0.3,0.4,0.5,0.6,0.9]
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,0.9)
plt.ylim(0,0.9)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('RSA_wtVSmut.jpg')

#----------------
# PROBABILITY DENSITY CURVE
fig = figure()
mu1, std1 = stats.norm.fit(x)
mu2, std2 = stats.norm.fit(y)
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k',label='Wild types (fit results: mu=%.2f,std=%.2f)'%(mu1, std1))
plt.plot(x2, p2, 'r',label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Solvent accessibility predicted values')
plt.ylabel('Frequency')
plt.xlim(0,0.9)
plt.ylim(0,4)
plt.legend(loc='upper right')
fig.savefig('histo_netsurfp_missense_wtVSmut.png')

# STATS
miss=[]
[miss.append(x - y) for x, y in zip(x, y)]
#KOLMOGOROV-SMINORV TEST:
stats.kstest(miss,'norm') # (D,pvalue) = (0.44913569824019062, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) # (T, pvalue) = (10720.0, 0.01473848472842257)
#So we reject H0 -> There is a significant difference between wt and mut

#-----------------------------------------------------------------------------
# RSA ENVIRONNEMENT
#-----------------------------------------------------------------------------
#-----------------
# SCATTER PLOT
RSA_wt = getColumn(file_wt,5,'\t')
RSA_mut = getColumn(file_mut,5,'\t')
RSA_wt.pop(0)
RSA_mut.pop(0)

x,y=[],[]
for i in range(0,len(RSA_wt)):
    if RSA_wt[i]=='NA':
        x.append(np.nan)
    else:
        x.append(float(RSA_wt[i]))
for i in range(0,len(RSA_mut)):
    if RSA_mut[i]=='NA':
        y.append(np.nan)
    else:
        y.append(float(RSA_mut[i]))


fig = plt.figure()
a=b=[0,0.2,0.3,0.4,0.5,0.6,0.9]
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,0.9)
plt.ylim(0,0.9)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('RSA_envt_wtVSmut.jpg')

#----------------
# PROBABILITY DENSITY CURVE
fig = figure()
mu1, std1 = stats.norm.fit(x)
mu2, std2 = stats.norm.fit(y)
bins = np.linspace(0, 99, 30)
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1) #Probability density function
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k',label='Wild types (fit results: mu=%.2f,std=%.2f)'%(mu1, std1))
plt.plot(x2, p2, 'r',label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Solvent accessibility predicted values')
plt.ylabel('Frequency')
plt.xlim(0,0.9)
plt.ylim(0,5)
plt.legend(loc='upper right')
fig.savefig('histo_netsurfp_missense_envt_wtVSmut.png')

# STATS
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') #(D,pvalue) = (0.47876635892857411, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) #-> (T, pvalue) = (13107.5, 0.17394709845400314)
#So we do not reject H0 -> There is no significant difference between wt and mut

#-----------------------------------------------------------------------------
# OUTLIERS FOR RSA (12)
#-----------------------------------------------------------------------------
RSA_wt = getColumn(file_wt,4,'\t')
RSA_mut = getColumn(file_mut,4,'\t')
RSA_wt.pop(0)
RSA_mut.pop(0)
RSAe_wt = getColumn(file_wt,5,'\t')
RSAe_mut = getColumn(file_mut,5,'\t')
RSAe_wt.pop(0)
RSAe_mut.pop(0)
variant_liste = getColumn(file_wt,0,'\t')
output = open('netsurfp_outliers.csv','w')
output.write('ID,RSA_wt,RSA_mut,difference,RSA_envt_wt,RSA_envt_mut,difference_envt\n')
for i in range(0,len(RSA_wt)):
    for j in range(0,len(RSA_mut)):
        if i==j:
            if RSA_wt[i]!='NA'and RSA_mut[j]!='NA':
                if (abs(float(RSA_wt[i])-float(RSA_mut[j]))) > 0.1:
                    output.write(variant_liste[i+1] + ',' + RSA_wt[i] + ',' + RSA_mut[j] + ',' + str(abs(float(RSA_wt[i])-float(RSA_mut[j]))) + ',' + RSAe_wt[i] + ',' + RSAe_mut[i] + ',' + str(abs(float(RSAe_wt[i])-float(RSAe_mut[j]))) + '\n')

output.close()

#-----------------------------------------------------------------------------
# RSA depending on Z-score
#-----------------------------------------------------------------------------
#-----------------
# SCATTER PLOT
Zscore_wt = getColumn(file_wt,6,'\t')
Zscore_mut = getColumn(file_mut,6,'\t')
Zscore_wt.pop(0)
Zscore_mut.pop(0)
RSA_wt = getColumn(file_wt,4,'\t')
RSA_mut = getColumn(file_mut,4,'\t')
RSA_wt.pop(0)
RSA_mut.pop(0)

ID = getColumn(file_wt,0,'\t')
ID.pop(0)

x_pos,x_neg,y_pos,y_neg=[],[],[],[]
IDwt_pos,IDwt_neg = [],[]
for i in range(0,len(RSA_wt)):
    if float(Zscore_wt[i])>=0:
        x_pos.append(float(RSA_wt[i]))
        IDwt_pos.append(ID[i])
    else:
        x_neg.append(float(RSA_wt[i]))
        IDwt_neg.append(ID[i])

IDmut_pos,IDmut_neg = [],[]
for i in range(0,len(RSA_mut)):
    if ID[i] in IDwt_pos:
        y_pos.append(float(RSA_mut[i]))
        IDmut_pos.append(ID[i])
    else:
        y_neg.append(float(RSA_mut[i]))
        IDmut_neg.append(ID[i])


# Z-score > 0 for wild types
fig = plt.figure()
a=b=[0,0,0.9]
plt.scatter(x_pos, y_pos,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,0.9)
plt.ylim(0,0.9)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('RSA_wtVSmut_zscore_pos.jpg')

#outliers (4)
output = open('netsurfp_outliers_zscore_pos.csv','w')
output.write('ID,RSA_wt,RSA_mut,difference\n')
for i in range(0,len(x_pos)):
    for j in range(0,len(y_pos)):
        if i==j:
            if (abs(float(x_pos[i])-float(y_pos[j]))) > 0.1:
                output.write(IDwt_pos[i] + ',' + str(x_pos[i]) + ',' + str(y_pos[j]) + ',' + str(abs(float(x_pos[i])-float(y_pos[j]))) + '\n')

output.close()

# Z-score < 0 fot wild types
fig = plt.figure()
a=b=[0,0,0.9]
plt.scatter(x_neg, y_neg,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,0.9)
plt.ylim(0,0.9)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('RSA_wtVSmut_zscore_neg.jpg')


#-----------------------------------------------------------------------------
# SECONDARY STRUCTURE
#-----------------------------------------------------------------------------
#Proba_a_helix
pa_wt = getColumn(file_wt,9,'\t')
pa_mut = getColumn(file_mut,9,'\t')
pa_wt.pop(0)
pa_mut.pop(0)

xa,ya=[],[]
for i in range(0,len(pa_wt)):
    if pa_wt[i]=='NA':
        xa.append(np.nan)
    else:
        xa.append(pa_wt[i])
for i in range(0,len(pa_mut)):
    if pa_mut[i]=='NA':
        ya.append(np.nan)
    else:
        ya.append(pa_mut[i])

             
fig = plt.figure()
a=b=[0,1]
plt.scatter(xa, ya,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('Proba_a_helix_wtVSmut.jpg')

#---------------------------
#Proba_b_strand
pb_wt = getColumn(file_wt,10,'\t')
pb_mut = getColumn(file_mut,10,'\t')
pb_wt.pop(0)
pb_mut.pop(0)

xb,yb=[],[]
for i in range(0,len(pb_wt)):
    if pb_wt[i]=='NA':
        xb.append(np.nan)
    else:
        xb.append(pb_wt[i])
for i in range(0,len(pb_mut)):
    if pb_mut[i]=='NA':
        yb.append(np.nan)
    else:
        yb.append(pb_mut[i])
             
fig = plt.figure()
a=b=[0,1]
plt.scatter(xb, yb,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('Proba_b_strand_wtVSmut.jpg')

#-----------------------------
#Proba_coil
pc_wt = getColumn(file_wt,11,'\t')
pc_mut = getColumn(file_mut,11,'\t')
pc_wt.pop(0)
pc_mut.pop(0)

xc,yc=[],[]
for i in range(0,len(pc_wt)):
    if pc_wt[i]=='NA':
        xc.append(np.nan)
    else:
        xc.append(pc_wt[i])
for i in range(0,len(pc_mut)):
    if pc_mut[i]=='NA':
        yc.append(np.nan)
    else:
        yc.append(pc_mut[i])
             
fig = plt.figure()
a=b=[0,1]
plt.scatter(xc, yc,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('Proba_coil_wtVSmut.jpg')

#-----------------------------
# BAR PLOTS
# Change probability by letter A (alpha-helix), B (beta-strand) or C (coil)
struct_wt = []
for i in range(0,len(xa)):
    proba = max(xa[i],xb[i],xc[i])
    if proba == xa[i]:
        struct_wt.append('A')
    elif proba == xb[i]:
        struct_wt.append('B')
    elif proba == xc[i]:
        struct_wt.append('C')

struct_mut = []
for i in range(0,len(ya)):
    proba = max(ya[i],yb[i],yc[i])
    if proba == ya[i]:
        struct_mut.append('A')
    elif proba == yb[i]:
        struct_mut.append('B')
    elif proba == yc[i]:
        struct_mut.append('C')

# difference
df1 = pd.read_csv(file_wt,'\t')
df1['struct_wt']=struct_wt
df2 = pd.read_csv(file_mut,'\t')
df2['struct_mut']=struct_mut
df3 = pd.concat([df1['ID'],df1['Gene_name'],df1['Variant_effect'],df1['struct_wt'],
                 df2['struct_mut']], axis=1)
df1List = df3['struct_wt'].tolist()
df2List = df3['struct_mut'].tolist()

newlist=[]
for i in range(0,len(df1List)):
    for j in range(0,len(df2List)):
        if i==j:
            if df1List[i]!=df2List[i]:
                newlist.append(str(df1List[i])+' to '+ str(df2List[i]))
            else:
                newlist.append('0')

df3['difference']=''
df3['difference'] = newlist
df3 = df3[df3.difference != '0']
df3.to_csv('netsurfp_diff_struct.csv', index=False)

dfdiff = df3.groupby(['difference']).count()
dfdiff = dfdiff.drop(['Gene_name','Variant_effect','struct_wt','struct_mut'], 1)
dfdiff.plot(kind='bar',legend=False,color='k')
plt.ylabel('Number of variants')
plt.savefig('barplotdiff_struct_netsurfp.png')


# Comparison wild types VS Deleterious DIDA mutants
N = 3
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
wt = (struct_wt.count('A'),struct_wt.count('B'),struct_wt.count('C'))
rects1 = ax.bar(ind, wt, width, color='grey')
mut = (struct_mut.count('A'),struct_mut.count('B'),struct_mut.count('C'))
rects2 = ax.bar(ind + width, mut, width, color='r')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C'))
ax.set_xlabel('Secondary structure prediction')
ax.legend((rects1[0], rects2[0]), ('Wild types', 'Deleterious DIDA mutants'),loc='upper center')
fig.savefig('barplot_netsurfp_struct_missense.png')

stats.chi2_contingency(np.column_stack((wt,mut)))
    #(0.70700945913170021, 0.70222267137184713, 2, array([[ 108.5,  108.5],[  38. ,   38. ],[  94.5,   94.5]]))
