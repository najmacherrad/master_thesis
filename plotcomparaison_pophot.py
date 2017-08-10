# Pop/Hot

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

#Import files
file_pop = 'popmusic_results.csv'
file_hot = 'hotmusic_results.csv'
file_netsurfp = 'Netsurfpresults.csv' #241 variants

# Faire new column pour secondary structure
# T/S/C = coil , E/B = beta , H/G = alpha
pop = pd.read_csv(file_pop,'\t') #72 lines
hot = pd.read_csv(file_hot,'\t') #72 lines
pophot = pop.merge(hot)

def f(row):
    if row['Secondary_structure'] == 'H' or row['Secondary_structure'] == 'G':
        val = 'A'
    elif row['Secondary_structure'] == 'E' or row['Secondary_structure'] == 'B':
        val = 'B'
    else:
        val = 'C'
    return val

pophot['sec_struct_type'] = pophot.apply(f, axis=1)
netsurfp = pd.read_csv(file_netsurfp,'\t') #255 lines
netsurfp = netsurfp.drop(['Gene_name','Variant_effect','RSA_envt_wt','Class_assignment_mut','RSA_mut','RSA_envt'],1)
netsurfp = netsurfp.rename(columns={'ID':'VariantID'})
netsurfp = netsurfp.rename(columns={'RSA_wt':'RSA_netsurfp'})
netsurfp = netsurfp.rename(columns={'Class_assignment_wt':'Class_netsurfp'})
merged_pophot = pophot.merge(netsurfp, on='VariantID')
merged_pophot.to_csv('pophotmusic_compare.csv', index=False)


# -------------------------------------------------------
# ANALYSIS OF CHANGES IN THERMODYNAMIC AND THERMOSTABILITY
# -------------------------------------------------------
#DIGENIC DELETERIOUS DIDA
DDG1 = getColumn('pophotmusic_compare.csv',7,',')
DTm1 = getColumn('pophotmusic_compare.csv',8,',')
DDG1.pop(0)
DTm1.pop(0)

DDG1_,DTm1_=[],[]
for i in range(0,len(DDG1)):
    if DDG1[i]=='':
        DDG1_.append(np.nan)
    else:
        DDG1_.append(float(DDG1[i]))
for i in range(0,len(DTm1)):
    if DTm1[i]=='':
        DTm1_.append(np.nan)
    else:
        DTm1_.append(float(DTm1[i]))


# -------------------------------------------
# COMPARISON WITH DELETERIOUS HUMSAVAR/SWISSVAR
# -------------------------------------------
DDG2d = getColumn('3BIO_delet.csv',6,',')
DTm2d = getColumn('3BIO_delet.csv',5,',')
DDG2d.pop(0)
DTm2d.pop(0)

DDG2d_,DTm2d_=[],[]
for i in range(0,len(DDG2d)):
    if DDG2d[i]=='':
        DDG2d_.append(np.nan)
    else:
        DDG2d_.append(float(DDG2d[i]))
for i in range(0,len(DTm2d)):
    if DTm2d[i]=='':
        DTm2d_.append(np.nan)
    else:
        DTm2d_.append(float(DTm2d[i]))

#DDG
fig = figure()
mu1, std1 = stats.norm.fit(DDG1_)
mu2, std2 = stats.norm.fit(DDG2d_)
bins = np.linspace(-6, 6, 50)
plt.hist(DDG1_, bins, alpha=0.3, label='Digenic deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,histtype='stepfilled',color='r')
plt.hist(DDG2d_, bins, alpha=0.3, label='Monogenic deleterious HumsaVar/SwissVar mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,histtype='stepfilled',color='r',hatch='/')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'r--', linewidth=2)
plt.xlabel('Thermodynamic stability DDG')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
plt.ylim(0,1)
fig.savefig('histoDDG_DIDAdeleter.png')

#MANN-WHITNEY:
stats.ranksums(DDG1_,DDG2d_) # (U,p-value) = (-5.3324556260236875, 9.6893451446496889e-08)

#DTm
fig = figure()
mu1, std1 = stats.norm.fit(DTm1_)
mu2, std2 = stats.norm.fit(DTm2d_)
bins = np.linspace(-15, 15, 50)
plt.hist(DTm1_, bins, alpha=0.3, label='Digenic deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,histtype='stepfilled',color='r')
plt.hist(DTm2d_, bins, alpha=0.3, label='Monogenic deleterious HumsaVar/SwissVar mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,histtype='stepfilled',color='r',hatch='/')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'r--', linewidth=2)
plt.xlabel('Thermostability DTm')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
plt.ylim(0,0.6)
fig.savefig('histoDTm_DIDAdeleter.png')

#MANN-WHITNEY:
stats.ranksums(DTm1_,DTm2d) # (U,p-value) = (4.3200335636790577, 1.5600547121398696e-05)


# -------------------------------------------
# COMPARISON WITH NEUTRAL HUMSAVAR/SWISSVAR
# -------------------------------------------
DDG2n = getColumn('3BIO_neut.csv',6,',')
DTm2n = getColumn('3BIO_neut.csv',5,',')
DDG2n.pop(0)
DTm2n.pop(0)

DDG2n_,DTm2n_=[],[]
for i in range(0,len(DDG2n)):
    if DDG2n[i]=='':
        DDG2n_.append(np.nan)
    else:
        DDG2n_.append(float(DDG2n[i]))
for i in range(0,len(DTm2n)):
    if DTm2n[i]=='':
        DTm2n_.append(np.nan)
    else:
        DTm2n_.append(float(DTm2n[i]))

#DDG
fig = figure()
mu1, std1 = stats.norm.fit(DDG1_)
mu2, std2 = stats.norm.fit(DDG2n_)
bins = np.linspace(-6, 6, 50)
plt.hist(DDG1_, bins, alpha=0.3, label='Digenic deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,histtype='stepfilled',color='red')
plt.hist(DDG2n_, bins, alpha=0.3, label='Monogenic neutral HumsaVar/SwissVar mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,histtype='stepfilled',color='blue',hatch='/')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'b--', linewidth=2)
plt.xlabel('Thermodynamic stability DDG')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
plt.ylim(0,1)
fig.savefig('histoDDG_DIDAneutral.png')

#MANN-WHITNEY:
stats.ranksums(DDG1_,DDG2n_) # (U,p-value) = (-0.73797835568169945, 0.46052760196902542)

#DTm
fig = figure()
mu1, std1 = stats.norm.fit(DTm1_)
mu2, std2 = stats.norm.fit(DTm2n_)
bins = np.linspace(-15, 15, 50)
plt.hist(DTm1_, bins, alpha=0.3, label='Digenic deleterious DIDA mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,histtype='stepfilled',color='red')
plt.hist(DTm2n_, bins, alpha=0.3, label='Monogenic neutral HumsaVar/SwissVar mutants \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,histtype='stepfilled',color='blue',hatch='/')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'r', linewidth=2)
plt.plot(x2, p2, 'b--', linewidth=2)
plt.xlabel('Thermostability DTm')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
plt.ylim(0,0.6)
fig.savefig('histoDTm_DIDAneutral.png')

#MANN-WHITNEY:
stats.ranksums(DTm1_,DTm2n_) # (U,p-value) = (0.22989610759092682, 0.81817250102183126)


# -------------------------------------------
# SOLVENT ACCESSIBILITY REAL VALUES
# -------------------------------------------

#-------------------------------------------
# deleterious DIDA mutants VS deleterious HumsaVar/SwissVar mutants
accDIDA = getColumn('pophotmusic_compare.csv',6,',')
accDel = getColumn('3BIO_delet.csv',4,',')
accDIDA.pop(0)
accDel.pop(0)

accD,accL=[],[]
for i in range(0,len(accDIDA)): #72
    if accDIDA[i]=='':
        accD.append(np.nan)
    else:
        accD.append(float(accDIDA[i])/100.0)
for i in range(0,len(accDel)): #4164
    if accDel[i]=='':
        accL.append(np.nan)
    else:
        accL.append(float(accDel[i])/100.0)

fig = figure()
mu1, std1 = stats.norm.fit(accD)
mu2, std2 = stats.norm.fit(accL)
bins = np.linspace(0, 6, 200)
plt.hist(accD, bins, alpha=0.3, label='Wild types (DIDA) \n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,color='grey')
plt.hist(accL, bins, alpha=0.3, label='Wild types (deleterious HumsaVar/SwissVar) \n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,color='grey',hatch='/')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k', linewidth=2)
plt.plot(x2, p2, 'k--', linewidth=2)
plt.xlabel('Solvent Accessibility real values')
plt.ylabel('Frequency')
#plt.ylim(0,10)
plt.xlim(0,1)
plt.legend(loc='upper right')
fig.savefig('histo_accVSdel.png')

#MANN-WHITNEY:
stats.ranksums(accD,accL) # (U,p-value) = (4.3427769767821776, 1.4069293980496851e-05)
# Reject H0
# The distributions of two sets of variables have a difference

#-------------------------------------------
# deleterious DIDA mutants VS neutral HumsaVar/SwissVar mutants
accDIDA = getColumn('pophotmusic_compare.csv',6,',')
accNeut = getColumn('3BIO_neut.csv',4,',')
accDIDA.pop(0)
accNeut.pop(0)

accD,accN=[],[]
for i in range(0,len(accDIDA)): #72
    if accDIDA[i]=='':
        accD.append(np.nan)
    else:
        accD.append(float(accDIDA[i])/100.0)
for i in range(0,len(accNeut)): #1342
    if accNeut[i]=='':
        accN.append(np.nan)
    else:
        accN.append(float(accNeut[i])/100.0)

fig = figure()
mu1, std1 = stats.norm.fit(accD)
mu2, std2 = stats.norm.fit(accN)
bins = np.linspace(0, 6, 200)
plt.hist(accD, bins, alpha=0.3, label='Wild types (DIDA)\n(fit results: mu=%.2f,std=%.2f)'%(mu1, std1),normed=True,color='grey')
plt.hist(accN, bins, alpha=0.3, label='Wild types (neutral HumsaVar/SwissVar)\n(fit results: mu=%.2f,std=%.2f)'%(mu2, std2),normed=True,color='grey',hatch='-')
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k', linewidth=2)
plt.plot(x2, p2, 'k--', linewidth=2)
plt.xlabel('Solvent Accessibility real values')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
plt.xlim(0,1)
fig.savefig('histo_accVSneut.png')

#MANN-WHITNEY:
stats.ranksums(accD,accN) # (U,p-value) = (-0.86922059236054039, 0.38472648620148686)
# Do not reject H0
# The distributions of two sets of variables have no difference

# -------------------------------------------
# wild types DIDA
# REAL VS PREDICTION
# solvent accessibility real values VS RSA predicted values by netsurfp
pop = getColumn('pophotmusic_compare.csv',6,',')
netsurfp = getColumn('pophotmusic_compare.csv',11,',')
pop.pop(0)
netsurfp.pop(0)

pop_,net_=[],[]
for i in range(0,len(netsurfp)):
    if netsurfp[i]=='':
        net_.append(np.nan)
    else:
        net_.append(float(netsurfp[i]))
for i in range(0,len(pop)):
    if pop[i]=='':
        pop_.append(np.nan)
    else:
        pop_.append(float(pop[i])/100.0)

fig = plt.figure()
a=b=[0,1]
plt.scatter(net_, pop_,edgecolor = 'none',c='k')
plt.plot(a,b,'r-')
plt.axis('on')
plt.grid('on')
plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('Solvent accessibility predicted values by NetSurfP')
plt.ylabel('Solvent accessibility real values')
fig.savefig('PopHotVSnetsurfp.jpg')

#----------------
# PROBABILITY DENSITY CURVE
fig = figure()
mu1, std1 = stats.norm.fit(pop_)
mu2, std2 = stats.norm.fit(net_)
xmin1, xmax1 = plt.xlim()
xmin2, xmax2 = plt.xlim()
x1 = np.linspace(xmin1, xmax1, 100)
x2 = np.linspace(xmin2, xmax2, 100)
p1 = stats.norm.pdf(x1, mu1, std1)
p2 = stats.norm.pdf(x2, mu2, std2)
plt.plot(x1, p1, 'k--', label='Real values\n(fit results: mu=%.2f,std=%.2f)'% (mu1, std1))
plt.plot(x2, p2, 'k', label='Predicted values by NetSurfP\n(fit results:mu=%.2f,std=%.2f)'%(mu2, std2))
plt.xlabel('Solvent accessibility')
plt.ylabel('Frequency')
plt.xlim(0,1)
#plt.ylim(0,1)
plt.legend(loc='upper right')
fig.savefig('histo_accVSRSA.png')

diff=[]
[diff.append(a_i - b_i) for a_i, b_i in zip(pop_, net_)]
#KOLMOGOROV-SMINORV:
stats.kstest(diff,'norm') #(D,pvalue) = (0.33500812974131394, 1.0441387732207374e-07)
#So we reject H0 -> not normal distribution
#WILCOXON:
stats.wilcoxon(diff) #(T, pvalue) = (1251.0, 0.723686128044988)
#So we reject H0 -> There is a significant difference between real and predicted


#------------------------------------------------------------------------------------
# SECONDARY STRUCTURE
#------------------------------------------------------------------------------------
# Faire new column pour secondary structure
# T/S/C = coil , E/B = beta , H/G = alpha
deletHum = pd.read_csv('3BIO_delet_struct.csv','\t') #4164
neutHum = pd.read_csv('3BIO_neut_struct.csv','\t') #1342

def f(row):
    if row['secStruc'] == 'H' or row['secStruc'] == 'G':
        val = 'A'
    elif row['secStruc'] == 'E' or row['secStruc'] == 'B':
        val = 'B'
    else:
        val = 'C'
    return val

deletHum['secStruc_2'] = deletHum.apply(f, axis=1)
deletHum.to_csv('3BIO_delet_struct_2.csv', index=False)

neutHum['secStruc_2'] = neutHum.apply(f, axis=1)
neutHum.to_csv('3BIO_neut_struct_2.csv', index=False)

#-------------------------------------------
# deleterious DIDA mutants VS deleterious HumsaVar/SwissVar mutants
structDIDA = getColumn('pophotmusic_compare.csv',9,',')
structDel = getColumn('3BIO_delet_struct_2.csv',9,',')
structDIDA.pop(0)
structDel.pop(0)

N = 3
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
DIDA_struct = (structDIDA.count('A')/72.0,structDIDA.count('B')/72.0,structDIDA.count('C')/72.0)
rects1 = ax.bar(ind, DIDA_struct, width, color='grey')
delet_struct = (structDel.count('A')/4164.0,structDel.count('B')/4164.0,structDel.count('C')/4164.0)
rects2 = ax.bar(ind + width, delet_struct, width, color='grey',hatch='/')
ax.set_ylabel('Frequency')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C'))
ax.set_xlabel('Secondary structure real values')
plt.ylim(0,0.6)
ax.legend((rects1[0], rects2[0]), ('Wild types (DIDA)', 'Wild types (deleterious HumsaVar/SwissVar)'))
fig.savefig('barplot_structVSdelet.png')

stats.chi2_contingency(np.column_stack((DIDA_struct,delet_struct))) #(0.017735112259088252, 0.9911716446889175, 2, array([[ 0.33345341,  0.33345341],[ 0.23098783,  0.23098783],[ 0.43555876,  0.43555876]]))

#-------------------------------------------
# deleterious DIDA mutants VS neutral HumsaVar/SwissVar mutants
structDIDA = getColumn('pophotmusic_compare.csv',9,',')
structNeut = getColumn('3BIO_neut_struct_2.csv',9,',')
structDIDA.pop(0)
structNeut.pop(0)

N = 3
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
DIDA_struct = (structDIDA.count('A')/72.0,structDIDA.count('B')/72.0,structDIDA.count('C')/72.0)
rects1 = ax.bar(ind, DIDA_struct, width, color='grey')
neut_struct = (structNeut.count('A')/1342.0,structNeut.count('B')/1342.0,structNeut.count('C')/1342.0)
rects2 = ax.bar(ind + width, neut_struct, width, color='grey',hatch='-')
ax.set_ylabel('Frequency')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C'))
ax.set_xlabel('Secondary structure real values')
plt.ylim(0,0.6)
ax.legend((rects1[0], rects2[0]), ('Wild types (DIDA)', 'Wild types (neutral HumsaVar/SwissVar)'))
fig.savefig('barplot_structVSneut.png')

stats.chi2_contingency(np.column_stack((DIDA_struct,neut_struct))) #(0.0093352043114804193, 0.99534327417050061, 2, array([[ 0.35668157,  0.35668157],[ 0.20489733,  0.20489733],[ 0.4384211 ,  0.4384211 ]]))

# -------------------------------------------
# wild types DIDA
# REAL VS PREDICTION
secstruct = pd.read_csv('pophotmusic_compare.csv',',')
netsurfpSTRUCT = pd.read_csv('netsurfp_diff_struct.csv',',')
netsurfpSTRUCT = netsurfpSTRUCT.rename(columns = {'ID':'VariantID'})
merged = netsurfpSTRUCT.merge(secstruct,on='VariantID')
merged = merged.drop({'struct_mut','difference','UniprotID','PDBID','Type','Chain','Secondary_structure','Solvent_accessibility','DDG','DTm','Class_netsurfp','RSA_netsurfp','dsysmap'},1)

df1List = merged['sec_struct_type'].tolist()
df2List = merged['struct_wt'].tolist()

newlist=[]
for i in range(0,len(df1List)):
    for j in range(0,len(df2List)):
        if i==j:
            if df1List[i]!=df2List[i]:
                newlist.append(str(df1List[i])+' to '+ str(df2List[i]))
            else:
                newlist.append('0')

merged['difference']=''
merged['difference'] = newlist
merged.to_csv('secstruct_count.csv', index=False)
merged = merged[merged.difference != '0']
merged.to_csv('secstruct_diff_realVSnetsurfp.csv', index=False)

dfdiff = merged.groupby(['difference']).count()
dfdiff = dfdiff.drop(['Gene_name','Variant_effect','struct_wt','sec_struct_type'], 1)
dfdiff.plot(kind='bar',legend=False,color='k')
plt.ylabel('Number of variants')
plt.savefig('barplotdiff_secstruct.png')

# Comparison real VS predicted
N = 3
ind = np.arange(N)  # the x locations for the groups
width = 0.30       # the width of the bars
fig, ax = plt.subplots()
real = (df1List.count('A'),df1List.count('B'),df1List.count('C'))
rects1 = ax.bar(ind, real, width, color='grey',hatch='.')
pred = (df2List.count('A'),df2List.count('B'),df2List.count('C'))
rects2 = ax.bar(ind + width, pred, width, color='grey')
ax.set_ylabel('Number of variants')
ax.set_xticks(ind + width)
ax.set_xticklabels(('A', 'B','C'))
ax.set_xlabel('Secondary structure prediction')
plt.ylim(0,45)
ax.legend((rects1[0], rects2[0]), ('Real for wild types (DIDA)', 'Prediction by NetSurfP for wild types (DIDA)'))
fig.savefig('barplot_realVSnetsurfp_struct.png')

stats.chi2_contingency(np.column_stack((real,pred)))
    #(1.1544011544011543, 0.56146795476357725, 2, array([[ 24. ,  24. ],[ 16.5,  16.5],[ 31.5,  31.5]]))

