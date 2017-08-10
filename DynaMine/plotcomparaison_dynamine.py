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
file_wt = 'dynamineresults_wt.csv'
file_mut = 'dynamineresults_mut.csv'

#-----------------------------------------------------------------------------
# FELIBILITY S2
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
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
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('pred_wtVSmut.jpg')

#-----------------------------------------------------------------------------
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
plt.plot(x2, p2, 'r',label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Flexibility S2 predicted values')
plt.ylabel('Frequency')
plt.xlim(0,1.2)
plt.ylim(0,5)
plt.legend(loc='upper right')
fig.savefig('histo_missense_wtVSmut.png')

# STATS
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') #(D,pvalue) = (0.45919964220308951, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) # (T, pvalue) = (10126.5, 3.9248500732405356e-05)
#So we reject H0 -> There is a significant difference between wt and mut


#-----------------------------------------------------------------------------
# FLEXIBILITY ENVIRONMENT
#-----------------------------------------------------------------------------
#-----------------
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
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('predenvt_wtVSmut.jpg')

#-----------------------------------------------------------------------------
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
plt.plot(x2, p2, 'r',label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Flexibility S2 predicted values')
plt.ylabel('Frequency')
plt.xlim(0,1.2)
plt.ylim(0,5)
plt.legend(loc='upper right')
fig.savefig('histo_missense_wtVSmut_envt.png')


# STATS
miss2=[]
[miss2.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss2,'norm') # (D,pvalue) = (0.47448511340469512, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss2) # (T, pvalue) = (10105.0, 7.5535095218014586e-05)
#So we reject H0 -> There is a significant difference between wt and mut

