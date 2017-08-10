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

#Importfiles
file_wt = 'netsurfpresultsNEW2_wt.csv'
file_mut = 'netsurfpresultsNEW2_1kg.csv'

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
plt.ylabel('Neutral 1KGP mutants')
fig.savefig('RSA_wtVS1kg.jpg')

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
plt.plot(x2, p2, 'b',label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Solvent accessibility predicted values')
plt.ylabel('Frequency')
plt.xlim(0,0.9)
plt.ylim(0,4)
plt.legend(loc='upper right')
fig.savefig('histo_netsurfp_missense_wtVS1kg.png')


#missense_wt - missense_mut
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') #(D,pvalue) = (0.42761364158461712, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) #-> (T, pvalue) = (1403683.5, 0.020490035411006691)
#So we reject H0 -> There is a significant difference between wt and mut

#-----------------------------------------------------------------------------
# RSA ENVIRONNEMENT
#-----------------------------------------------------------------------------
#-----------------
# SCATTER PLOT
#-----------------------------------------------------------------------------
#RSA_envt
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
plt.ylabel('Neutral 1KGP mutants')
fig.savefig('RSA_envt_wtVS1kg.jpg')

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
plt.plot(x2, p2, 'b',label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Solvent accessibility predicted values')
plt.ylabel('Frequency')
plt.xlim(0,0.9)
plt.ylim(0,5)
plt.legend(loc='upper right')
fig.savefig('histo_netsurfp_missense_envt_wtVS1kg.png')

# STATS
miss=[]
[miss.append(a_i - b_i) for a_i, b_i in zip(x, y)]
#KOLMOGOROV-SMINORV:
stats.kstest(miss,'norm') # (D,pvalue) = (0.45460452749063657, 0.0)
#So we reject H0 -> not normal distribution
#WILCOXON TEST:
stats.wilcoxon(miss) #-> (T, pvalue) = (1548668.0, 0.43701657073338696)
#So we do not reject H0 -> There is no significant difference between wt and mut

#-----------------------------------------------------------------------------
# OUTLIERS FOR RSA (270)
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
variant_liste.pop(0)
output = open('netsurfp_outliers_1kg.csv','w')
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
a=b=[0,0,0.8]
plt.scatter(x_pos, y_pos,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,0.9)
plt.ylim(0,0.9)
plt.xlabel('Wild types')
plt.ylabel('Neutral 1KGP mutants')
fig.savefig('RSA_wtVS1kg_zscore_pos.jpg')

#outliers (41)
output = open('netsurfp1kg_outliers_zscore_pos.csv','w')
output.write('ID,RSA_wt,RSA_mut,difference\n')
for i in range(0,len(x_pos)):
    for j in range(0,len(y_pos)):
        if i==j:
            if (abs(float(x_pos[i])-float(y_pos[j]))) > 0.1:
                output.write(IDwt_pos[i] + ',' + str(x_pos[i]) + ',' + str(y_pos[j]) + ',' + str(abs(float(x_pos[i])-float(y_pos[j]))) + '\n')

output.close()

# Z-score < 0 fot wild types
fig = plt.figure()
a=b=[0,0,0.8]
plt.scatter(x_neg, y_neg,edgecolor = 'none', c= 'k')
plt.plot(a,b,'r-')
plt.grid('on')
plt.xlim(0,0.8)
plt.ylim(0,0.8)
plt.xlabel('Wild type residues')
plt.ylabel('Mutant residues')
fig.savefig('RSA_wtVS1kg_zscore_neg.jpg')


#-----------------------------------------------------------------------------
# RSA : COMPARISON deleterious DIDA mutants VS neutral 1KGP mutants
#-----------------------------------------------------------------------------
file_DIDA = 'netsurfpresults_mut_DIDA.csv'

RSA_DIDA = getColumn(file_DIDA,4,'\t')
RSA_1kg = getColumn(file_mut,4,'\t')
RSA_DIDA.pop(0)
RSA_1kg.pop(0)

xRSA,yRSA=[],[]
for i in range(0,len(RSA_DIDA)): #241
    if RSA_DIDA[i]=='NA':
        xRSA.append(np.nan)
    else:
        xRSA.append(float(RSA_DIDA[i]))
for i in range(0,len(RSA_1kg)): #2516
    if RSA_1kg[i]=='NA':
        yRSA.append(np.nan)
    else:
        yRSA.append(float(RSA_1kg[i]))


fig = figure()
mu1, std1 = stats.norm.fit(xRSA)
mu2, std2 = stats.norm.fit(yRSA)
bins = np.linspace(-0.3, 1, 35)
plt.hist(xRSA, bins, alpha=0.3, label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,color='red')
plt.hist(yRSA, bins, alpha=0.3, label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,color='blue')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'b', linewidth=2)
plt.xlabel('Solvent accessibility predicted values')
plt.ylabel('Frequency')
plt.ylim(0,5)
plt.xlim(-0.3,1)
plt.legend(loc='upper right')
fig.savefig('histoRSA_DIDA1kg.png')


#MANN-WHITNEY:
stats.ranksums(xRSA,yRSA) # (U,p-value) = (-5.995280821744239, 2.0313410214210638e-09)
# Reject H0
# The distributions of two sets of variables have a difference


#-----------------------------------------------------------------------------
# RSA ENVIRONMENT: COMPARISON deleterious DIDA mutants VS neutral 1KGP mutants
#-----------------------------------------------------------------------------
RSA_DIDA = getColumn(file_DIDA,5,'\t')
RSA_1kg = getColumn(file_mut,5,'\t')
RSA_DIDA.pop(0)
RSA_1kg.pop(0)

xRSA,yRSA=[],[]
for i in range(0,len(RSA_DIDA)): #241
    if RSA_DIDA[i]=='NA':
        xRSA.append(np.nan)
    else:
        xRSA.append(float(RSA_DIDA[i]))
for i in range(0,len(RSA_1kg)): #2516
    if RSA_1kg[i]=='NA':
        yRSA.append(np.nan)
    else:
        yRSA.append(float(RSA_1kg[i]))

fig = figure()
mu1, std1 = stats.norm.fit(xRSA)
mu2, std2 = stats.norm.fit(yRSA)
bins = np.linspace(-0.3, 1, 35)
plt.hist(xRSA, bins, alpha=0.3, label='Deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,color='red')
plt.hist(yRSA, bins, alpha=0.3, label='Neutral 1KGP mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,color='blue')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'b', linewidth=2)
plt.xlabel('Solvent accessibility predicted values')
plt.ylabel('Frequency')
plt.ylim(0,5)
plt.xlim(-0.3,1)
plt.legend(loc='upper right')
fig.savefig('histoRSAenvt_DIDA1kg.png')

#MANN-WHITNEY:
stats.ranksums(xRSA,yRSA) # (U,p-value) = (-7.4005610929180445, 1.356102615569394e-13)
# Reject H0
# The distributions of two sets of variables have a difference



