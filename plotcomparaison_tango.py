# TANGO
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
import operator
from numpy import *
import petl as etl
import operator
import pandas as pd
import re
import csv

def getColumn(filename, column,deli):
    results = csv.reader(open(filename), delimiter=deli)
    return [result[column] for result in results]

#Import files
file_wt = 'tangoresults_wt.csv'
file_mut = 'tangoresults_mut.csv'

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  AGGREGATION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#--------------------------------------
# SCATTER PLOT
pred_wt = getColumn(file_wt,3,'\t')
pred_mut = getColumn(file_mut,3,'\t')
pred_wt.pop(0)
pred_mut.pop(0)

x,y=[],[]
for i in range(0,len(pred_wt)): #max=99.986
    if pred_wt[i]=='NA':
        x.append(np.nan)
    else:
        x.append(float(pred_wt[i]))
for i in range(0,len(pred_mut)): #max=99.986
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
fig.savefig('agg_wtVSmut.jpg')

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
plt.xlabel('Aggregation conformation predicted values')
plt.ylabel('Frequency')
plt.xlim(0,100)
plt.ylim(0,0.025)
plt.legend(loc='upper right')
fig.savefig('histagg_missense.png')


#missense_wt - missense_mut
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') # (D,pvalue) = (0.32541604726493034, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) # (T, pvalue) = (1965.0, 0.075068140452097282)
#So we do not reject H0 -> There is no significant difference between wt and mut


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  AGGREGATION ENVIRONMENT
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#--------------------------------------
# SCATTER PLOT
pred_wt = getColumn(file_wt,4,'\t')
pred_mut = getColumn(file_mut,4,'\t')
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
a=b=[0,100]
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(-1,101)
plt.ylim(-1,101)
plt.xlabel('Wild types')
plt.ylabel('Deleterious DIDA mutants')
fig.savefig('aggenvt_wtVSmut.jpg')

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
plt.xlabel('Aggregation conformation predicted values')
plt.ylabel('Frequency')
plt.xlim(0,100)
plt.ylim(0,0.030)
plt.legend(loc='upper right')
fig.savefig('histaggcenvt_missense.png')

#missense_wt - missense_mut
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') # (D,pvalue) = (0.27621231080632236, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) # (T, pvalue) = (3945.0, 0.46774505385307708)
#So we do not reject H0 -> There is no significant difference between wt and mut


#-----------------------------------------------------------------------------
# OUTLIERS FOR AGGREGATION (25)
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
output = open('aggregation_outliers.csv','w')
output.write('ID,agg_wt,agg_mut,difference,agg_envt_wt,agg_envt_mut,difference_envt\n')
for i in range(0,len(pred_wt)):
    for j in range(0,len(pred_mut)):
        if i==j:
            if pred_wt[i]!='NA'and pred_mut[j]!='NA':
                if (abs(float(pred_wt[i])-float(pred_mut[j]))) > 20:
                    output.write(variant_liste[i+1] + ',' + pred_wt[i] + ',' + pred_mut[j] + ',' + str(abs(float(pred_wt[i])-float(pred_mut[j]))) + ',' + pred_envt_wt[i] + ',' + pred_envt_mut[i] + ',' + str(abs(float(pred_envt_wt[i])-float(pred_envt_mut[j]))) + '\n')

output.close()


#-------------------------------------------------------------------------------
# COMPARISON WITH NETSURFP RSA
#-------------------------------------------------------------------------------
tango_wt = pd.read_csv(file_wt,'\t')
tango_mut = pd.read_csv(file_mut,'\t')
tango_wt['DTango'] = ''
tango_wt['DTango'] = tango_wt.Aggregation - tango_mut.Aggregation
tango_wt['DTango_envt'] = tango_wt.Aggregation_envt - tango_mut.Aggregation_envt
tango_wt = tango_wt.drop(['Aggregation','Aggregation_envt'], 1)
tango_wt.to_csv('tangoresults_compare.csv', index=False)

N_wt = pd.read_csv('netsurfpresults_wt.csv','\t')
N_mut = pd.read_csv('netsurfpresults_mut.csv','\t')
N_wt['DNetsurfp'] = N_wt.RSA - N_mut.RSA
N_wt['DNetsurfp_envt'] = N_wt.RSA_envt - N_mut.RSA_envt
N_wt = N_wt.drop(['Class_assignment','RSA','RSA_envt','Z_score','Z_score_envt','ASA','Proba_a_helix','Proba_B_strand','Proba_Coil'], 1)
N_wt.to_csv('netsurfpresults_compare.csv', index=False)

tango = getColumn('tangoresults_compare.csv',3,',')
tango_envt =  getColumn('tangoresults_compare.csv',4,',')
netsurfp = getColumn('netsurfpresults_compare.csv',3,',')
netsurfp_envt = getColumn('netsurfpresults_compare.csv',4,',')
tango.pop(0)
netsurfp.pop(0)
tango_envt.pop(0)
netsurfp_envt.pop(0)

#RESIDUE
x,y=[],[]
for i in range(0,len(netsurfp)): #min=-0.183 and max=0.302
    if netsurfp[i]=='':
        x.append(np.nan)
    else:
        x.append(float(netsurfp[i]))
for i in range(0,len(tango)): #min=-98.23 and max=98.842
    if tango[i]=='':
        y.append(np.nan)
    else:
        y.append(float(tango[i]))
             
fig = plt.figure()
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.grid('on')
plt.xlim(-0.4,0.4)
plt.ylim(-100,100)
plt.xlabel('delta(Solvent accessibility prediction) by NetSurfP')
plt.ylabel('delta(Aggregation conformation prediction) by Tango')
fig.savefig('TangoVSnetsurfp.jpg')

#ENVIRONMENT
x,y=[],[]
for i in range(0,len(netsurfp_envt)): #min=-0.183 and max=0.302
    if netsurfp_envt[i]=='':
        x.append(np.nan)
    else:
        x.append(float(netsurfp_envt[i]))
for i in range(0,len(tango_envt)): #min=-98.23 and max=98.842
    if tango_envt[i]=='':
        y.append(np.nan)
    else:
        y.append(float(tango_envt[i]))

fig = plt.figure()
plt.scatter(x, y,edgecolor = 'none', c= 'k')
plt.grid('on')
plt.xlim(-0.4,0.4)
plt.ylim(-100,100)
plt.xlabel('delta(Solvent accessibility prediction) by NetSurfP')
plt.ylabel('delta(Aggregation conformation prediction) by Tango')
fig.savefig('TangoVSnetsurfp_envt.jpg')

