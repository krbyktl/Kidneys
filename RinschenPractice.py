# -*- coding: utf-8 -*-

"""
Spyder Editor

This is a temporary script file.

@author: leokt
Practice data pushing with mpkCCD phosphoproteome database from Rinschen et al 2010

"""

#import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#import data
path = "/Users/kirbyleo/Documents/NIH/Datasets/"
file = "mpkCCD_phosphoproteomic_db2.csv"
fileloc = path + file
df = pd.read_csv(fileloc)
df_filt = df.iloc[15:]


#generate dataframe
seq = df_filt["mpkCCD Cell Phosphoproteomic Database (mCPD)"]
gene = df_filt["Unnamed: 2"]
peptide = df_filt["Unnamed: 3"]
protein = df_filt["Unnamed: 4"]
phos = df_filt["Unnamed: 5"]
exp1 = df_filt["Unnamed: 7"]
for a in range(0,len(exp1)):
    if exp1.iloc[a] == '>10':
        exp1.iloc[a] = '10'
    if exp1.iloc[a] == '<-10':
        exp1.iloc[a] = '-10'
    if float(exp1.iloc[a]) > 20:
        exp1.iloc[a] = '20'
    if float(exp1.iloc[a]) < -20:
        exp1.iloc[a] = '-20'
exp2 = df_filt["Unnamed: 8"]
for b in range(0,len(exp2)):
    if exp2.iloc[b] == '>10':
        exp2.iloc[b] = '10'
    if exp2.iloc[b] == '<-10':
        exp2.iloc[b] = '-10'
    if float(exp2.iloc[b]) > 20:
        exp2.iloc[b] = '20'
    if float(exp2.iloc[b]) < -20:
        exp2.iloc[b] = '-20'
exp3 = df_filt["Unnamed: 9"]
for c in range(0,len(exp3)):
    if exp3.iloc[c] == '>10':
        exp3.iloc[c] = '10'
    if exp3.iloc[c] == '<-10':
        exp3.iloc[c] = '-10'
    if float(exp3.iloc[c]) > 20:
        exp3.iloc[c] = '20'
    if float(exp3.iloc[c]) < -20:
        exp3.iloc[c] = '-20'

dic = {"RefSeq number": seq, "Gene Symbol": gene, "Peptide sequence": peptide,
       "Protein name": protein, "Phosphosite(s)": phos, "log2(dDAVP/control): Exp. 1": exp1, 
       "log2(dDAVP/control): Exp. 2": exp2, "log2(dDAVP/control): Exp. 3": exp3}
iden = pd.DataFrame(data = dic)

iden["Exp Avg"] = iden[["log2(dDAVP/control): Exp. 1","log2(dDAVP/control): Exp. 2",
    "log2(dDAVP/control): Exp. 3"]].astype(float).mean(axis=1, skipna=True)

#plot figure 1
dist_changes = sns.distplot(iden["Exp Avg"], bins = 120, kde = False, color = "black")
plt.xlabel("$log_2$(dDAVP/control)")
plt.ylabel("Number of Unique Peptides")
plt.xlim(-8,8)
plt.savefig("RinschenFig1", dpi = 300)

#amino acid sorting
amino_acids = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']



