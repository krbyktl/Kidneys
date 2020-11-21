# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 20:24:49 2020

@author: leokt
"""

#import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


#import data
fileloc = "J:\Depot - dDAVP-time course - Kirby\Analysis\\201110_Total_Protein.xlsx"
df_cntl = pd.read_excel(fileloc, sheet_name="median norm", header=2)

#variables

df_cntl['T1_log_1'] = np.log2(df_cntl['V1_Replicate1']/df_cntl['C1_Replicate1'])
df_cntl['T1_log_2'] = np.log2(df_cntl['V1_Replicate2']/df_cntl['C1_Replicate2'])
df_cntl['T1_log_3'] = np.log2(df_cntl['V1_Replicate3']/df_cntl['C1_Replicate3'])
df_cntl['T2_log_1'] = np.log2(df_cntl['V2_Replicate1']/df_cntl['C2_Replicate1'])
df_cntl['T2_log_2'] = np.log2(df_cntl['V2_Replicate2']/df_cntl['C2_Replicate2'])
df_cntl['T2_log_3'] = np.log2(df_cntl['V2_Replicate3']/df_cntl['C2_Replicate3'])
df_cntl['T3_log_1'] = np.log2(df_cntl['V3_Replicate1']/df_cntl['C3_Replicate1'])
df_cntl['T3_log_2'] = np.log2(df_cntl['V3_Replicate2']/df_cntl['C3_Replicate2'])
df_cntl['T3_log_3'] = np.log2(df_cntl['V3_Replicate3']/df_cntl['C3_Replicate3'])
df_cntl['T4_log_1'] = np.log2(df_cntl['V4_Replicate1']/df_cntl['C4_Replicate1'])
df_cntl['T4_log_2'] = np.log2(df_cntl['V4_Replicate2']/df_cntl['C4_Replicate2'])
df_cntl['T4_log_3'] = np.log2(df_cntl['V4_Replicate3']/df_cntl['C4_Replicate3'])


df_cntl['T1_log'] = df_cntl[['T1_log_1','T1_log_2','T1_log_3']].mean(axis=1)
df_cntl['T2_log'] = df_cntl[['T2_log_1','T2_log_2','T2_log_3']].mean(axis=1)
df_cntl['T3_log'] = df_cntl[['T3_log_1','T3_log_2','T3_log_3']].mean(axis=1)
df_cntl['T4_log'] = df_cntl[['T4_log_1','T4_log_2','T4_log_3']].mean(axis=1)


#abundance of unqiue peptides non-phos data
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1,4,figsize=(20,5))
plt.subplots_adjust(wspace=0, hspace=0)
dist_T1 = sns.distplot(df_cntl['T1_log'], axlabel = "$log_2$(dDAVP/control)", bins = 180, kde = False, color = "black", ax = ax1)
ax1.set_ylabel("Number of Unique Proteins")
ax1.set_title("dDAVP vs. Vehicle 1 min")
dist_T2 = sns.distplot(df_cntl['T2_log'], axlabel = "$log_2$(dDAVP/control)", bins = 180, kde = False, color = "blue", ax = ax2)
ax2.set_yticklabels([])
ax2.set_title("dDAVP vs. Vehicle 2 min")
dist_T3 = sns.distplot(df_cntl['T3_log'], axlabel = "$log_2$(dDAVP/control)", bins = 180, kde = False, color = "purple", ax = ax3)
ax3.set_yticklabels([])
ax3.set_title("dDAVP vs. Vehicle 5 min")
dist_T4 = sns.distplot(df_cntl['T4_log'], axlabel = "$log_2$(dDAVP/control)", bins = 180, kde = False, color = "red", ax = ax4)
ax4.set_yticklabels([])
ax4.set_title("dDAVP vs. Vehicle 15 min")
plt.setp((ax1, ax2, ax3, ax4), ylim = (0,2500))
plt.savefig("abundance_distribution", dpi = 300)
plt.close()

#total heatmap time course
df_hm = pd.DataFrame({'T1':df_cntl['T1_log'], 'T2':df_cntl['T2_log'],
                      'T3':df_cntl['T3_log'],'T4':df_cntl['T4_log']})
df_hm = pd.DataFrame(df_hm.values,
                     columns = ['1','2','5','15'],
                     index = df_cntl['Gene Symbol'])
fig, ax = plt.subplots(figsize=(11, 9))
hm = sns.heatmap(df_hm, cmap="Blues", 
                 cbar_kws={'label': '$log_2$(dDAVP/control)'})
ax.set_yticklabels([])
plt.xlabel('Time (min)')
plt.ylabel('')
plt.title('Time Course Heatmap Total Abundance', loc='left')
plt.savefig("abundance_heatmap", dpi = 1200)
plt.close()

#volcano plot of peptide abundance FDR vs. ratio(dDAVP/vehicle)
#15 min
df_T4 = pd.read_excel(fileloc, sheet_name="15 min")

fig, ax = plt.subplots(figsize=(12,10))
plt_T4 = sns.scatterplot(data = df_T4, x="OVERALL RATIO",  y="minus log p",
                         color = 'lightgrey')
outr = df_T4.query('`minus log p` > 1.303 & `OVERALL RATIO` > 0.054397916*2')
outl = df_T4.query('`minus log p` > 1.303 & `OVERALL RATIO` < -0.054397916*2')
outr.index = range(len(outr.index))
outl.index = range(len(outl.index))
plt.scatter(outr['OVERALL RATIO'], outr['minus log p'], color = "red")
plt.scatter(outl['OVERALL RATIO'], outl['minus log p'], color = "red")
for i in range(len(outr)):
    ax.annotate(outr['Gene Symbol'][i], 
                 xy = (outr['OVERALL RATIO'][i], outr['minus log p'][i]),
                 fontsize = 10,
                 textcoords="offset points")
for i in range(len(outl)):
    ax.annotate(outl['Gene Symbol'][i], 
                 xy = (outl['OVERALL RATIO'][i], outl['minus log p'][i]),
                 fontsize = 10,
                 textcoords="offset points")
sns.despine()
plt.ylabel('-$log_{10}$(p-value)', fontsize = 18)
plt.xlabel('$log_2$(dDAVP/control)', fontsize = 18)
plt.title("Volcano Plot (15 min) n = 5444", fontsize = 24)
plt.ylim(0,4.5)
plt.xlim(-1.3,1.3)
ax.axhline(1.303, color = 'blue', ls='--')
plt.text(-1.2,1.35,'1.303')
ax.axvline(0.054397916*2, color = 'blue', ls='--')
plt.text(0.12,4.1,'0.109', rotation=90)
ax.axvline(-0.054397916*2, color = 'blue', ls='--')
plt.text(-0.10,4.1,'-0.109', rotation=90)
plt.savefig("Volcano total 15", dpi = 1200)

#5 min
df_T3 = pd.read_excel(fileloc, sheet_name="5 min")

fig, ax = plt.subplots(figsize=(12,10))
plt_T3 = sns.scatterplot(data = df_T3, x="OVERALL RATIO",  y="minus log p",
                         color = 'lightgrey')
outr = df_T3.query('`minus log p` > 1.303 & `OVERALL RATIO` > 0.050789416*2')
outl = df_T3.query('`minus log p` > 1.303 & `OVERALL RATIO` < -0.050789416*2')
outr.index = range(len(outr.index))
outl.index = range(len(outl.index))
plt.scatter(outr['OVERALL RATIO'], outr['minus log p'], color = "red")
plt.scatter(outl['OVERALL RATIO'], outl['minus log p'], color = "red")
for i in range(len(outr)):
    ax.annotate(outr['Gene Symbol'][i], 
                 xy = (outr['OVERALL RATIO'][i], outr['minus log p'][i]),
                 fontsize = 10,
                 textcoords="offset points")
for i in range(len(outl)):
    ax.annotate(outl['Gene Symbol'][i], 
                 xy = (outl['OVERALL RATIO'][i], outl['minus log p'][i]),
                 fontsize = 10,
                 textcoords="offset points")
sns.despine()
plt.ylabel('-$log_{10}$(p-value)', fontsize = 18)
plt.xlabel('$log_2$(dDAVP/control)', fontsize = 18)
plt.title("Volcano Plot (5 min) n = 5444", fontsize = 24)
plt.ylim(0,4.5)
plt.xlim(-1.3,1.3)
ax.axhline(1.303, color = 'blue', ls='--')
plt.text(-1.2,1.35,'1.303')
ax.axvline(0.050789416*2, color = 'blue', ls='--')
plt.text(0.12,4.1,'0.102', rotation=90)
ax.axvline(-0.050789416*2, color = 'blue', ls='--')
plt.text(-0.10,4.1,'-0.102', rotation=90)
plt.savefig("Volcano total 5", dpi = 1200)

#2 min
df_T2 = pd.read_excel(fileloc, sheet_name="2 min")

fig, ax = plt.subplots(figsize=(12,10))
plt_T2 = sns.scatterplot(data = df_T2, x="OVERALL RATIO",  y="minus log p",
                         color = 'lightgrey')
outr = df_T2.query('`minus log p` > 1.303 & `OVERALL RATIO` > 0.044202462*2')
outl = df_T2.query('`minus log p` > 1.303 & `OVERALL RATIO` < -0.044202462*2')
outr.index = range(len(outr.index))
outl.index = range(len(outl.index))
plt.scatter(outr['OVERALL RATIO'], outr['minus log p'], color = "red")
plt.scatter(outl['OVERALL RATIO'], outl['minus log p'], color = "red")
for i in range(len(outr)):
    ax.annotate(outr['Gene Symbol'][i], 
                 xy = (outr['OVERALL RATIO'][i], outr['minus log p'][i]),
                 fontsize = 10,
                 textcoords="offset points")
for i in range(len(outl)):
    ax.annotate(outl['Gene Symbol'][i], 
                 xy = (outl['OVERALL RATIO'][i], outl['minus log p'][i]),
                 fontsize = 10,
                 textcoords="offset points")
sns.despine()
plt.ylabel('-$log_{10}$(p-value)', fontsize = 18)
plt.xlabel('$log_2$(dDAVP/control)', fontsize = 18)
plt.title("Volcano Plot (2 min) n = 5444", fontsize = 24)
plt.ylim(0,4.5)
plt.xlim(-1.3,1.3)
ax.axhline(1.303, color = 'blue', ls='--')
plt.text(-1.2,1.35,'1.303')
ax.axvline(0.044202462*2, color = 'blue', ls='--')
plt.text(0.11,4.1,'0.088', rotation=90)
ax.axvline(-0.044202462*2, color = 'blue', ls='--')
plt.text(-0.09,4.1,'-0.088', rotation=90)
plt.savefig("Volcano total 2", dpi = 1200)

#1 min
df_T1 = pd.read_excel(fileloc, sheet_name="1 min")

fig, ax = plt.subplots(figsize=(12,10))
plt_T1 = sns.scatterplot(data = df_T1, x="OVERALL RATIO",  y="minus log p",
                         color = 'lightgrey')
outr = df_T1.query('`minus log p` > 1.303 & `OVERALL RATIO` > 0.05433744*2')
outl = df_T1.query('`minus log p` > 1.303 & `OVERALL RATIO` < -0.05433744*2')
outr.index = range(len(outr.index))
outl.index = range(len(outl.index))
plt.scatter(outr['OVERALL RATIO'], outr['minus log p'], color = "red")
plt.scatter(outl['OVERALL RATIO'], outl['minus log p'], color = "red")
for i in range(len(outr)):
    ax.annotate(outr['Gene Symbol'][i], 
                 xy = (outr['OVERALL RATIO'][i], outr['minus log p'][i]),
                 fontsize = 10,
                 textcoords="offset points")
for i in range(len(outl)):
    ax.annotate(outl['Gene Symbol'][i], 
                 xy = (outl['OVERALL RATIO'][i], outl['minus log p'][i]),
                 fontsize = 10,
                 textcoords="offset points")
sns.despine()
plt.ylabel('-$log_{10}$(p-value)', fontsize = 18)
plt.xlabel('$log_2$(dDAVP/control)', fontsize = 18)
plt.title("Volcano Plot (1 min) n = 5444", fontsize = 24)
plt.ylim(0,4.5)
plt.xlim(-1.3,1.3)
ax.axhline(1.303, color = 'blue', ls='--')
plt.text(-1.2,1.35,'1.303')
ax.axvline(0.05433744*2, color = 'blue', ls='--')
plt.text(0.12,4.1,'0.109', rotation=90)
ax.axvline(-0.05433744*2, color = 'blue', ls='--')
plt.text(-0.10,4.1,'-0.109', rotation=90)
plt.savefig("Volcano total 1", dpi = 1200)





