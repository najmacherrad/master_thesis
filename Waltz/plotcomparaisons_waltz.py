# Waltz
# Compare results between wild type and mutant

# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from scipy import stats
from pylab import plot, show, savefig, xlim, figure, \
                hold, ylim, legend, boxplot, setp, axes
import pylab
from numpy import *

def getColumn(filename, column,deli):
    results = csv.reader(open(filename), delimiter=deli)
    return [result[column] for result in results]

#import files
file_wt = 'waltzresults_wt.csv'
file_mut = 'waltzresults_mut.csv'
 
#------------------------------------
# AGGREGATION
#------------------------------------
#--------------------------------------
# SCATTER PLOT
pred_wt = getColumn(file_wt,3,'\t') 
pred_mut = getColumn(file_mut,3,'\t') 
pred_wt.pop(0)
pred_mut.pop(0)

x,y=[],[]
for i in range(0,len(pred_wt)): #max=98.662207
    if pred_wt[i]=='NA':
        x.append(np.nan)
    else:
        x.append(float(pred_wt[i]))
for i in range(0,len(pred_mut)): #max=99.665552
    if pred_mut[i]=='NA':
        y.append(np.nan)
    else:
        y.append(float(pred_mut[i]))
             
fig = plt.figure()
a=b=[0,100]
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(-1,101)
plt.ylim(-1,101)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('waltz_wtVSmut.jpg')
                    
#----------------
# PROBABILITY DENSITY CURVE
fig = figure()
mu1, std1 = stats.norm.fit(x)
mu2, std2 = stats.norm.fit(y)
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, 100, 100)
x2 = np.linspace(xmin2, 100, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k',label='Wild types (fit results: mu=%.2f,std=%.2f)'%(mu1, std1))
plt.plot(x2, p2, 'r',label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Aggregation conformation predicted values (amylogenic regions)')
plt.ylabel('Frequency')
plt.xlim(0,100)
#plt.ylim(0,0.0)
plt.legend(loc='upper right')
fig.savefig('histwaltz_missense.png')

#missense_wt - missense_mut
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') # (D,pvalue) = (0.3552063996073398, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) # (T, pvalue) = (4898.0, 0.29548245005836105)
#So we do not reject H0 -> There is no significant difference between wt and mut

#--------------------------------------
# AGGREGATION ENVIRONMENT
#--------------------------------------
#--------------------------------------
# SCATTER PLOT
pred_wt = getColumn(file_wt,4,'\t') 
pred_mut = getColumn(file_mut,4,'\t') 
pred_wt.pop(0)
pred_mut.pop(0)

x,y=[],[]
for i in range(0,len(pred_wt)): #max=98.662207
    if pred_wt[i]=='NA':
        x.append(np.nan)
    else:
        x.append(float(pred_wt[i]))
for i in range(0,len(pred_mut)): #max=98.996656
    if pred_mut[i]=='NA':
        y.append(np.nan)
    else:
        y.append(float(pred_mut[i]))
             
fig = plt.figure()
a=b=[0,100]
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(-1,101)
plt.ylim(-1,101)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('waltz_envt_wtVSmut.jpg')

#--------------------------------------
# HISTOGRAM
fig = figure()
mu1, std1 = stats.norm.fit(x)
mu2, std2 = stats.norm.fit(y)
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, 100, 100)
x2 = np.linspace(xmin2, 100, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k',label='Wild types (fit results: mu=%.2f,std=%.2f)'%(mu1, std1))
plt.plot(x2, p2, 'r',label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Aggregation conformation predicted values (amylogenic regions)')
plt.ylabel('Frequency')
plt.xlim(0,100)
plt.ylim(0,0.06)
plt.legend(loc='upper right')
fig.savefig('histwaltzenvt_missense.png')

#missense_wt - missense_mut
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') # (D,pvalue) = (0.34964202670995748, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) #-> (T, pvalue) = (8711.0, 0.55024961096028457)
#So we do not reject H0 -> There is no significant difference between wt and mut

#-----------------------------------------------------------------------------
# OUTLIERS FOR AGGREGATION ()
#-----------------------------------------------------------------------------
pred_wt = getColumn(file_wt,3,'\t')
pred_mut = getColumn(file_mut,3,'\t')
pred_wt.pop(0)
pred_mut.pop(0)
pred_envt_wt = getColumn(file_wt,4,'\t')
pred_envt_mut = getColumn(file_mut,4,'\t')
pred_envt_wt.pop(0)
pred_envt_mut.pop(0)
variant_liste = getColumn(file_wt,0,'\t')
output = open('waltz_outliers.csv','w')
output.write('ID,agg_wt,agg_mut,difference,agg_envt_wt,agg_envt_mut,difference_envt\n')
for i in range(0,len(pred_wt)):
    for j in range(0,len(pred_mut)):
        if i==j:
            if pred_wt[i]!='NA'and pred_mut[j]!='NA':
                if (abs(float(pred_wt[i])-float(pred_mut[j]))) > 20:
                    output.write(variant_liste[i+1] + ',' + pred_wt[i] + ',' + pred_mut[j] + ',' + str(abs(float(pred_wt[i])-float(pred_mut[j]))) +  ',' + pred_envt_wt[i] + ',' + pred_envt_mut[i] + ',' + str(abs(float(pred_envt_wt[i])-float(pred_envt_mut[j]))) + '\n')

output.close()


#-------------------------------------------------------------------------------
#COMPARISON WITH NETSURFP RSA
#-------------------------------------------------------------------------------
W_wt = pd.read_csv(file_wt,'\t')
W_mut = pd.read_csv(file_mut,'\t')
W_wt['DWaltz'] = ''
W_wt['DWaltz'] = W_wt.aggregation - W_mut.aggregation
W_wt['DWaltz_envt'] = ''
W_wt['DWaltz_envt'] = W_wt.aggregation_envt - W_mut.aggregation_envt
W_wt = W_wt.drop(['aggregation','aggregation_envt'], 1)
W_wt.to_csv('waltzresults_compare.csv', index=False)

#RESIDUE
waltz = getColumn('waltzresults_compare.csv',3,',')
waltz.pop(0)
netsurfp = getColumn('netsurfpresults_compare.csv',3,',')
netsurfp.pop(0)

x,y=[],[]
for i in range(0,len(netsurfp)): #min=-0.183 and max=0.302
    if netsurfp[i]=='':
        x.append(np.nan)
    else:
        x.append(float(netsurfp[i]))
for i in range(0,len(waltz)): #min=-98.862207 and max=98.327759
    if waltz[i]=='':
        y.append(np.nan)
    else:
        y.append(float(waltz[i]))

fig = plt.figure()
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.grid('on')
plt.xlim(-0.4,0.4)
plt.ylim(-100,100)
plt.xlabel('delta(Solvent accessibility prediction) by NetSurfP')
plt.ylabel('delta(Aggregation conformation prediction) by Waltz')
fig.savefig('WaltzVSnetsurfp.jpg')

#ENVIRONMENT
waltz_envt = getColumn('waltzresults_compare.csv',4,',')
waltz_envt.pop(0)
netsurfp_envt = getColumn('netsurfpresults_compare.csv',4,',')
netsurfp_envt.pop(0)

x,y=[],[]
for i in range(0,len(netsurfp_envt)): #min=-0.183 and max=0.302
    if netsurfp_envt[i]=='':
        x.append(np.nan)
    else:
        x.append(float(netsurfp_envt[i]))
for i in range(0,len(waltz_envt)): #min=-98.862207 and max=98.327759
    if waltz_envt[i]=='':
        y.append(np.nan)
    else:
        y.append(float(waltz_envt[i]))

fig = plt.figure()
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.grid('on')
plt.xlim(-0.4,0.4)
plt.ylim(-100,100)
plt.xlabel('delta(Solvent accessibility prediction) by NetSurfP')
plt.ylabel('delta(Aggregation conformation prediction) by Waltz')
fig.savefig('WaltzVSnetsurfp_envt.jpg')

