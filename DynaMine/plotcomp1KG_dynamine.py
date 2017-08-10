#DynaMine
#Compare results between wild type and mutant

# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from scipy import stats
import pylab
from pylab import plot, show, savefig, xlim, figure, \
                hold, ylim, legend, boxplot, setp, axes

def getColumn(filename, column,deli):
    results = csv.reader(open(filename), delimiter=deli)
    return [result[column] for result in results]

#Importer les fichiers
file_wt = 'dynamineresultsNEW_wt.csv'
file_mut = 'dynamineresultsNEW_1kg.csv'

#-----------------------------------------------------------------------------
# FELIBILITY S2
#-----------------------------------------------------------------------------
#--------------
# SCATTER PLOT
pred_wt = getColumn(file_wt,3,'\t') 
pred_mut = getColumn(file_mut,3,'\t') 
pred_wt.pop(0)
pred_mut.pop(0)

x,y=[],[]
for i in range(0,len(pred_wt)):
    if pred_wt[i]=='NA':
        x.append(np.nan)
    else:
        x.append(float(pred_wt[i]))
for i in range(0,len(pred_mut)):
    if pred_mut[i]=='NA':
        y.append(np.nan)
    else:
        y.append(float(pred_mut[i]))
             
fig = plt.figure()
a=b=[0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2]
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,1.2)
plt.ylim(0,1.2)
plt.xlabel('Wild types')
plt.ylabel('Neutral 1KGP mutants')
fig.savefig('pred_wtVS1kg.jpg')

#----------------------------
# PROBABILITY DENSITY CURVE
fig = figure()
mu1, std1 = stats.norm.fit(x)
mu2, std2 = stats.norm.fit(y)
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, 1.2, 100)
x2 = np.linspace(xmin2, 1.2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k',label='Wild types (fit results: mu=%.2f,std=%.2f)'%(mu1, std1))
plt.plot(x2, p2, 'b',label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Flexibility S2 predicted values')
plt.ylabel('Frequency')
plt.xlim(0,1.2)
plt.ylim(0,5)
plt.legend(loc='upper right')
fig.savefig('histo_missense_wtVS1kg.png')

# STATS
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') # (D,pvalue) = (0.45873047720322435, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) #-> (T, pvalue) = (5982726.0, 6.1610257571015192e-82)
#So we reject H0 -> There is a significant difference between wt and mut

#-----------------------------------------------------------------------------
# FELIBILITY S2 ENVIRONMENT
#-----------------------------------------------------------------------------
#--------------
# SCATTER PLOT
predenvt_wt = getColumn(file_wt,4,'\t')
predenvt_mut = getColumn(file_mut,4,'\t')
predenvt_wt.pop(0)
predenvt_mut.pop(0)

x,y=[],[]
for i in range(0,len(predenvt_wt)):
    if predenvt_wt[i]=='NA':
        x.append(np.nan)
    else:
        x.append(float(predenvt_wt[i]))
for i in range(0,len(predenvt_mut)):
    if predenvt_mut[i]=='NA':
        y.append(np.nan)
    else:
        y.append(float(predenvt_mut[i]))

fig = plt.figure()
a=b=[0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2]
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,1.2)
plt.ylim(0,1.2)
plt.xlabel('Wild types')
plt.ylabel('Neutral 1KGP mutants')
fig.savefig('predenvt_wtVS1kg.jpg')

#----------------------------
# PROBABILITY DENSITY CURVE
fig = figure()
mu1, std1 = stats.norm.fit(x)
mu2, std2 = stats.norm.fit(y)
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, 1.2, 100)
x2 = np.linspace(xmin2, 1.2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k',label='Wild types (fit results: mu=%.2f,std=%.2f)'%(mu1, std1))
plt.plot(x2, p2, 'b',label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Flexibility S2 predicted values')
plt.ylabel('Frequency')
plt.xlim(0,1.2)
plt.ylim(0,5)
plt.legend(loc='upper right')
fig.savefig('histo_missense_wtVS1kg_envt.png')


# STATS
miss2=[]
[miss2.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss2,'norm') # (D,pvalue) = (0.46756747979693603, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss2) # (T, pvalue) = (6232277.0, 2.6686982056384321e-70)
#So we reject H0 -> There is a significant difference between wt and mut

#-----------------------------------------------------------------------------
# OUTLIERS FOR FLEXIBILITY (24)
#-----------------------------------------------------------------------------
dyn_wt = getColumn(file_wt,3,'\t')
dyn_mut = getColumn(file_mut,3,'\t')
dyn_wt.pop(0)
dyn_mut.pop(0)
dyne_wt = getColumn(file_wt,4,'\t')
dyne_mut = getColumn(file_mut,4,'\t')
dyne_wt.pop(0)
dyne_mut.pop(0)
variant_liste = getColumn(file_wt,0,'\t')
output = open('dynamine_outliers_1kg.csv','w')
output.write('ID,pred_wt,pred_mut,difference,pred_envt_wt,pred_envt_mut,difference_envt\n')
for i in range(0,len(dyn_wt)):
    for j in range(0,len(dyn_mut)):
        if i==j:
            if dyn_wt[i]!='NA'and dyn_mut[j]!='NA':
                if (abs(float(dyn_wt[i])-float(dyn_mut[j]))) > 0.25:
                    output.write(variant_liste[i+1] + ',' + dyn_wt[i] + ',' + dyn_mut[j] + ',' + str(abs(float(dyn_wt[i])-float(dyn_mut[j]))) + ',' + dyne_wt[i] + ',' + dyne_mut[i] + ',' + str(abs(float(dyne_wt[i])-float(dyne_mut[j]))) + '\n')

output.close()

#-----------------------------------------------------------------------------
# FLEXIBILITY : COMPARISON deleterious DIDA mutants VS neutral 1KGP mutants
#-----------------------------------------------------------------------------
file_DIDAmut = 'dynamineresults_mut.csv'

pred_DIDA = getColumn(file_DIDAmut,3,'\t')
pred_1kg = getColumn(file_mut,3,'\t')
pred_DIDA.pop(0)
pred_1kg.pop(0)

xpred,ypred=[],[]
for i in range(0,len(pred_DIDA): #241
    if pred_DIDA[i]=='NA':
        xpred.append(np.nan)
    else:
        xpred.append(float(pred_DIDA[i]))
for i in range(0,len(pred_1kg)): #5846
    if pred_1kg[i]=='NA':
        ypred.append(np.nan)
    else:
        ypred.append(float(pred_1kg[i]))

fig = figure()
mu1, std1 = stats.norm.fit(xpred)
mu2, std2 = stats.norm.fit(ypred)
bins = np.linspace(0, 1.2, 35)
plt.hist(xpred,bins,normed=True,alpha=0.3, color='red',label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1))
plt.hist(ypred,bins,normed=True,alpha=0.3, label='neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),color='blue')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'b', linewidth=2)
plt.xlabel('Flexibility S2 predicted values')
plt.ylabel('Frequency')
plt.xlim(0,1.2)
plt.ylim(0,5.5)
plt.legend(loc='upper left')
fig.savefig('histo_dynamine_DIDAVS1kg.png')

#MANN-WHITNEY:
stats.ranksums(xpred,ypred) # (U,p-value) = (2.9860551466029612, 0.0028260167537576602)
# Reject H0
# The distributions of two sets of variables have a difference


#-----------------------------------------------------------------------------
# FLEXIBILITY ENVT: COMPARISON deleterious DIDA mutants VS neutral 1KGP mutants
#-----------------------------------------------------------------------------
pred_DIDA = getColumn(file_DIDAmut,4,'\t')
pred_1kg = getColumn(file_mut,4,'\t')
pred_DIDA.pop(0)
pred_1kg.pop(0)

xpred,ypred=[],[]
for i in range(0,len(pred_DIDA)): #241
    if pred_DIDA[i]=='NA':
        xpred.append(np.nan)
    else:
        xpred.append(float(pred_DIDA[i]))
for i in range(0,len(pred_1kg)): #5846
    if pred_1kg[i]=='NA':
        ypred.append(np.nan)
    else:
        ypred.append(float(pred_1kg[i]))

fig = figure()
mu1, std1 = stats.norm.fit(xpred)
mu2, std2 = stats.norm.fit(ypred)
bins = np.linspace(0, 1.2, 35)
plt.hist(xpred, bins, alpha=0.3, label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,color='red')
plt.hist(ypred, bins, alpha=0.3, label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,color='blue')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'b', linewidth=2)
plt.xlabel('Flexibility S2 predicted values')
plt.ylabel('Frequency')
plt.legend(loc='upper left')
plt.xlim(0,1.2)
plt.ylim(0,5.5)
fig.savefig('histo_dynamineenvt_DIDA1kg.png')


#MANN-WHITNEY:
stats.ranksums(xpred,ypred) # (U,p-value) = (2.6912760375831581, 0.0071179272598568873)
# Reject H0
# The distributions of two sets of variables have a difference
