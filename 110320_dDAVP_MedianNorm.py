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
import statsmodels as sm


#import data
path = "C:/Users/leokt/Documents/Time course/"
file = "dDAVP_Time_Proteins_MedianNorm.xlsx"
fileloc = path + file
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
ax1.set_ylabel("Number of Unique Peptides Total Abundance")
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

#volcano plot of peptide abundance FDR vs. ratio(dDAVP/vehicle)
ttest_T1 = []
ttest_T2 = []
ttest_T3 = []
ttest_T4 = []


for i in range(len(df_cntl['V1_Replicate1'])):
    c1
    c2
    c3
    v1
    v2
    v3
    sm.stats.weightstats.ttest_ind()

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

#heatmap transcription factors?





