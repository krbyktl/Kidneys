# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 13:14:19 2020

@author: leokt
"""
#import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


#import data
fileloc = "J:\Depot - dDAVP-time course - Kirby\Analysis\\201112_Time_Phospho.xlsx"
df_phos = pd.read_excel(fileloc, sheet_name="Norm using Total Median Norm Fc")


#isolate peptide probabilities
prob_col = df_phos['Phospho (STY) Peptide Probabilities']
df_phos['STY Probabilities'] = np.arange(len(prob_col))

import re

for i in range(len(prob_col)):
    temp = re.findall(r'\d+', prob_col[i])
    probstr = ""
    for t in range(len(temp)):
        if temp[t] == '0':
            probstr = probstr + temp[t] + '.'
        elif temp[t] == '1':
            if t != len(temp)-1:
                probstr = probstr + temp[t] + ';'
            else:
                probstr = probstr + temp[t]
        else:
            if t != len(temp)-1:
                probstr += temp[t] + ';'
            else:
                probstr += temp[t]
    df_phos['STY Probabilities'][i] = probstr
    
#Concatenate gene symbol and site


#plot monophosphorylated peptides
#15 min 
#fileloc1 = "J:\Depot - dDAVP-time course - Kirby\Analysis\\201114_Phospho_15min.xlsx"
df_phos15 = pd.read_excel(fileloc, sheet_name="15 min monophos")

fig, ax = plt.subplots(figsize=(12,10))
outr = df_phos15.query('`minus log p` > 1.303 & `OVERALL RATIO_1` > 0.15135757*2')
outl = df_phos15.query('`minus log p` > 1.303 & `OVERALL RATIO_1` < -0.15135757*2')
outr = outr.append(outl)
outr.index = range(len(outr.index))
sns.scatterplot(data = df_phos15, x="OVERALL RATIO_1",  y="minus log p",
                         color = 'lightgrey')
outr['Phosphosites'] = outr['Gene Symbol'] + '_' + outr['Position in Protein (A)']
sns.scatterplot(data = outr,  x="OVERALL RATIO_1",  y="minus log p",
                hue = 'Phosphosites')
plt.legend(bbox_to_anchor=(1.4, 1), loc='upper right', ncol = 2, fontsize=9)
sns.despine()
plt.ylabel('-$log_{10}$(p-value)', fontsize = 18)
plt.xlabel('$log_2$(dDAVP/control)', fontsize = 18)
num_pep = len(df_phos15['Gene Symbol'])
sig_pep = len(outr['Gene Symbol'])
plt.title('Volcano Plot (15 min) n = ' + str(num_pep) + ', n* = ' + str(sig_pep), fontsize = 24)
ax.axhline(1.303, color = 'blue', ls='--')
plt.ylim(0,4.5)
plt.xlim(-2,2)
plt.text(-1.8,1.35,'1.303')
ax.axvline(0.15135757*2, color = 'blue', ls='--')
plt.text(0.313,4.1,'0.303', rotation=90)
ax.axvline(-0.15135757*2, color = 'blue', ls='--')
plt.text(-0.29,4.1,'-0.303', rotation=90)
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Scripts\Figures\\Volcano phos 15",
            bbox_inches='tight', dpi = 1200)
        

#5 min
#fileloc2 = "J:\Depot - dDAVP-time course - Kirby\Analysis\\201116_Phospho_5min.xlsx"
df_phos5 = pd.read_excel(fileloc, sheet_name="5 min monophos")       

fig, ax = plt.subplots(figsize=(12,10))
outr = df_phos5.query('`minus log p` > 1.303 & `OVERALL RATIO_1` > 0.130652376*2')
outl = df_phos5.query('`minus log p` > 1.303 & `OVERALL RATIO_1` < -0.130652376*2')
outr = outr.append(outl)
outr.index = range(len(outr.index))
sns.scatterplot(data = df_phos5, x="OVERALL RATIO_1",  y="minus log p",
                         color = 'lightgrey')
outr['Phosphosites'] = outr['Gene Symbol'] + '_' + outr['Position in Protein (A)']
sns.scatterplot(data = outr,  x="OVERALL RATIO_1",  y="minus log p",
                hue = 'Phosphosites')
plt.legend(bbox_to_anchor=(1.4, 1), loc='upper right', ncol = 2, fontsize=9)
sns.despine()
plt.ylabel('-$log_{10}$(p-value)', fontsize = 18)
plt.xlabel('$log_2$(dDAVP/control)', fontsize = 18)
num_pep = len(df_phos5['Gene Symbol'])
sig_pep = len(outr['Gene Symbol'])
plt.title('Volcano Plot (5 min) n = ' + str(num_pep) + ', n* = ' + str(sig_pep), fontsize = 24)
ax.axhline(1.303, color = 'blue', ls='--')
plt.ylim(0,4.5)
plt.xlim(-2,2)
plt.text(-1.8,1.35,'1.303')
ax.axvline(0.130652376*2, color = 'blue', ls='--')
plt.text(0.272,4.1,'0.261', rotation=90)
ax.axvline(-0.130652376*2, color = 'blue', ls='--')
plt.text(-0.252,4.1,'-0.261', rotation=90)
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Scripts\Figures\\Volcano phos 5",
            bbox_inches='tight', dpi = 1200)
        
#2 min
#fileloc3 = "J:\Depot - dDAVP-time course - Kirby\Analysis\\201116_Phospho_2min.xlsx"
df_phos2 = pd.read_excel(fileloc, sheet_name="2 min monophos")       

fig, ax = plt.subplots(figsize=(12,10))
outr = df_phos2.query('`minus log p` > 1.303 & `OVERALL RATIO_1` > 0.106992368*2')
outl = df_phos2.query('`minus log p` > 1.303 & `OVERALL RATIO_1` < -0.106992368*2')
outr = outr.append(outl)
outr.index = range(len(outr.index))
sns.scatterplot(data = df_phos2, x="OVERALL RATIO_1",  y="minus log p",
                         color = 'lightgrey')
outr['Phosphosites'] = outr['Gene Symbol'] + '_' + outr['Position in Protein (A)']
sns.scatterplot(data = outr,  x="OVERALL RATIO_1",  y="minus log p",
                hue = 'Phosphosites')
plt.legend(bbox_to_anchor=(1.4, 1), loc='upper right', ncol = 2, fontsize=9)
sns.despine()
plt.ylabel('-$log_{10}$(p-value)', fontsize = 18)
plt.xlabel('$log_2$(dDAVP/control)', fontsize = 18)
num_pep = len(df_phos2['Gene Symbol'])
sig_pep = len(outr['Gene Symbol'])
plt.title('Volcano Plot (2 min) n = ' + str(num_pep) + ', n* = ' + str(sig_pep), fontsize = 24)
ax.axhline(1.303, color = 'blue', ls='--')
plt.ylim(0,4.5)
plt.xlim(-2,2)
plt.text(-1.8,1.35,'1.303')
ax.axvline(0.106992368*2, color = 'blue', ls='--')
plt.text(0.224,4.1,'0.214', rotation=90)
ax.axvline(-0.106992368*2, color = 'blue', ls='--')
plt.text(-0.204,4.1,'-0.214', rotation=90)
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Scripts\Figures\\Volcano phos 2",
            bbox_inches='tight', dpi = 1200)

#1 min
#fileloc4 = "J:\Depot - dDAVP-time course - Kirby\Analysis\\201116_Phospho_1min.xlsx"
df_phos1 = pd.read_excel(fileloc, sheet_name="1 min monophos")       

fig, ax = plt.subplots(figsize=(12,10))
plt_T1 = sns.scatterplot(data = df_phos1, x="OVERALL RATIO_1",  y="minus log p",
                         color = 'lightgrey')
outr = df_phos1.query('`minus log p` > 1.303 & `OVERALL RATIO_1` > 0.10278922*2')
outl = df_phos1.query('`minus log p` > 1.303 & `OVERALL RATIO_1` < -0.10278922*2')
outr = outr.append(outl)
outr.index = range(len(outr.index))
sns.scatterplot(data = df_phos1, x="OVERALL RATIO_1",  y="minus log p",
                         color = 'lightgrey')
outr['Phosphosites'] = outr['Gene Symbol'] + '_' + outr['Position in Protein (A)']
sns.scatterplot(data = outr,  x="OVERALL RATIO_1",  y="minus log p",
                hue = 'Phosphosites')
plt.legend(bbox_to_anchor=(1.4, 1), loc='upper right', ncol = 2, fontsize=9)
sns.despine()
plt.ylabel('-$log_{10}$(p-value)', fontsize = 18)
plt.xlabel('$log_2$(dDAVP/control)', fontsize = 18)
num_pep = len(df_phos2['Gene Symbol'])
sig_pep = len(outr['Gene Symbol'])
plt.title('Volcano Plot (1 min) n = ' + str(num_pep) + ', n* = ' + str(sig_pep), fontsize = 24)
ax.axhline(1.303, color = 'blue', ls='--')
plt.ylim(0,4.5)
plt.xlim(-2,2)
plt.text(-1.8,1.35,'1.303')
ax.axvline(0.10278922*2, color = 'blue', ls='--')
plt.text(0.212,4.1,'0.202', rotation=90)
ax.axvline(-0.10278922*2, color = 'blue', ls='--')
plt.text(-0.192,4.1,'-0.202', rotation=90)
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Scripts\Figures\\Volcano phos 1",
            bbox_inches='tight', dpi = 1200)

        








