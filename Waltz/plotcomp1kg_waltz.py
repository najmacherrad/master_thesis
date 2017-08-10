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

def getColumn(filename, column,deli):
    results = csv.reader(open(filename), delimiter=deli)
    return [result[column] for result in results]

#import files
file_wt = 'waltzresultsNEW_wt.csv'
file_mut = 'waltzresultsNEW_1kg.csv'

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
plt.ylabel('Neutral 1KGP mutants')
fig.savefig('waltz_wtVS1kg.jpg')

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
plt.plot(x2, p2, 'b',label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Aggregation conformation predicted values (amylogenic regions)')
plt.ylabel('Frequency')
plt.xlim(0,100)
plt.ylim(0,0.045)
plt.legend(loc='upper right')
fig.savefig('histwaltz_missense1kg.png')

#missense_wt - missense_mut
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') # (D,pvalue) =  (0.37733219518495709, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) # (T, pvalue) =  (1948626.0, 0.00024177854220634842)
#So we reject H0 -> There is a significant difference between wt and mut

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
plt.ylabel('Neutral 1KGP mutants')
fig.savefig('waltz_envt_wtVS1kg.jpg')

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
plt.plot(x2, p2, 'b',label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Aggregation conformation predicted values (amylogenic regions)')
plt.ylabel('Frequency')
plt.xlim(0,100)
plt.ylim(0,0.06)
plt.legend(loc='upper right')
fig.savefig('histwaltzenvt_missense_1kg.png')

# STATS
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') # (D,pvalue) = (0.36028431579335313, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) # (T, pvalue) = (3834070.0, 2.138887727291201e-17)
#Donc il y a une différence significative entre wild type et mutant

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
output = open('waltz_outliers_1kg.csv','w')
output.write('ID,agg_wt,agg_mut,difference,agg_envt_wt,agg_envt_mut,difference_envt\n')
for i in range(0,len(pred_wt)):
    for j in range(0,len(pred_mut)):
        if i==j:
            if pred_wt[i]!='NA'and pred_mut[j]!='NA':
                if (abs(float(pred_wt[i])-float(pred_mut[j]))) > 20:
                    output.write(variant_liste[i+1] + ',' + pred_wt[i] + ',' + pred_mut[j] + ',' + str(abs(float(pred_wt[i])-float(pred_mut[j]))) + ',' + pred_envt_wt[i] + ',' + pred_envt_mut[i] + ',' + str(abs(float(pred_envt_wt[i])-float(pred_envt_mut[j]))) + '\n')

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
W_wt.to_csv('waltzresults_compare1kg.csv', index=False)

T = open('waltzresults_compare1kg.csv','r')
N = open('netsurfpresults_compare1kg.csv','r')
T2 = open('_waltzresults_compare1kg.csv','w')
N2 = open('_netsurfpresults_compare1kg.csv','w')
c1 = csv.reader(T, delimiter=',')
c2 = csv.reader(N, delimiter=',')
c3 = csv.writer(T2)
c4 = csv.writer(N2)
waltz = list(c1) #5'847
netsurfp = list(c2) #2'517
for i in waltz:
    for j in netsurfp:
        if i[0] == j[0]:
            c3.writerow(i)
            c4.writerow(j)

N.close()
T.close()
N2.close()
T2.close()

#RESIDUE
waltz = getColumn('_waltzresults_compare1kg.csv',3,',')
netsurfp = getColumn('_netsurfpresults_compare1kg.csv',3,',')
waltz.pop(0)
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
fig.savefig('WaltzVSnetsurfp_1kg.jpg')

#ENVIRONMENT
waltz_envt = getColumn('_waltzresults_compare1kg.csv',4,',')
netsurfp_envt = getColumn('_netsurfpresults_compare1kg.csv',4,',')
waltz_envt.pop(0)
netsurfp_envt.pop(0)

x,y=[],[]
for i in range(0,len(netsurfp_envt)): #min=-0.127 and max=0.08275
    if netsurfp_envt[i]=='':
        x.append(np.nan)
    else:
        x.append(float(netsurfp_envt[i]))
for i in range(0,len(waltz_envt)): #min=-93.645485 and max=97.993311
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
fig.savefig('WaltzVSnetsurfp_1kg_envt.jpg')

#-----------------------------------------------------------------------------
# AGGREGATION : COMPARISON deleterious DIDA mutants VS neutral 1KGP mutants
#-----------------------------------------------------------------------------
file_DIDAmut = 'waltzresults_mut.csv'

pred_DIDA = getColumn(file_DIDAmut,3,'\t')
pred_1kg = getColumn(file_mut,3,'\t')
pred_DIDA.pop(0)
pred_1kg.pop(0)

xpred,ypred=[],[]
for i in range(0,len(pred_DIDA)):
    if pred_DIDA[i]=='NA':
        xpred.append(np.nan)
    else:
        xpred.append(float(pred_DIDA[i]))
for i in range(0,len(pred_1kg)):
    if pred_1kg[i]=='NA':
        ypred.append(np.nan)
    else:
        ypred.append(float(pred_1kg[i]))


fig = figure()
mu1, std1 = stats.norm.fit(xpred)
mu2, std2 = stats.norm.fit(ypred)
bins = np.linspace(0, 100, 35)
plt.hist(xpred,bins,normed=True,alpha=0.3, color='r',label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1))
plt.hist(ypred,bins,normed=True,alpha=0.3, label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),color='blue')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'b', linewidth=2)
plt.xlabel('Aggregation conformation predicted values (amylogenic regions)')
plt.ylabel('log(Frequency)')
plt.xlim(0,100)
plt.ylim(1e-24,1e4)
plt.legend(loc='upper right')
plt.yscale('log')
fig.savefig('histo_waltz_DIDAVS1kg.png')


#MANN-WHITNEY:
stats.ranksums(xpred,ypred) # (U,p-value) = (2.8180561149873737, 0.0048315365001516033)
# Reject H0
# The distributions of two sets of variables have a difference

#-----------------------------------------------------------------------------
# AGGREGATION ENVT : COMPARISON deleterious DIDA mutants VS neutral 1KGP mutants
#-----------------------------------------------------------------------------
pred_DIDA = getColumn(file_DIDAmut,4,'\t')
pred_1kg = getColumn(file_mut,4,'\t')
pred_DIDA.pop(0)
pred_1kg.pop(0)

xpred,ypred=[],[]
for i in range(0,len(pred_DIDA)):
    if pred_DIDA[i]=='NA':
        xpred.append(np.nan)
    else:
        xpred.append(float(pred_DIDA[i]))
for i in range(0,len(pred_1kg)):
    if pred_1kg[i]=='NA':
        ypred.append(np.nan)
    else:
        ypred.append(float(pred_1kg[i]))


fig = figure()
mu1, std1 = stats.norm.fit(xpred)
mu2, std2 = stats.norm.fit(ypred)
bins = np.linspace(0, 100, 35)
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
plt.xlabel('Aggregation conformation predicted values (amylogenic regions)')
plt.ylabel('log(Frequency)')
plt.xlim(0,100)
plt.ylim(1e-24,1e4)
plt.legend(loc='upper right')
plt.yscale('log')
fig.savefig('histo_waltzenvt_DIDA1kg.png')

#MANN-WHITNEY:
stats.ranksums(xpred,ypred) # (U,p-value) = (3.0503709286051883, 0.0022855889607669021)
# Reject H0
# The distributions of two sets of variables have a difference
