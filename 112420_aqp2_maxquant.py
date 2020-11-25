# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 13:58:48 2020

@author: leokt
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

fileloc = "J:\Depot - dDAVP-time course - Kirby\Analysis\\Aqp2_Msms.xlsx"
df_aqp2 = pd.read_excel(fileloc, sheet_name="phospho_msms")

#normalization of intensities
df_TMT1 = df_aqp2[df_aqp2['Replicate'] == 'TMT1']
df_TMT1 = df_TMT1.rename(columns = {'Reporter intensity corrected 1':'C1',
                         'Reporter intensity corrected 2':'C2',
                         'Reporter intensity corrected 3':'C3',
                         'Reporter intensity corrected 4':'C4',
                         'Reporter intensity corrected 8':'V1',
                         'Reporter intensity corrected 9':'V2',
                         'Reporter intensity corrected 10':'V3',
                         'Reporter intensity corrected 11':'V4'})
TMT1 = ['C1','C2','C3','C4','V1','V2','V3','V4']
norm_factor1 = [1.01, 0.91, 0.96, 1.03, 1, 0.95, 1, 1]
for i in range(len(TMT1)):
    df_TMT1[TMT1[i]] = norm_factor1[i]*df_TMT1[TMT1[i]]
    
df_TMT2 = df_aqp2[df_aqp2['Replicate'] == 'TMT2']
df_TMT2 = df_TMT2.rename(columns = {'Reporter intensity corrected 1':'V1',
                         'Reporter intensity corrected 2':'V2',
                         'Reporter intensity corrected 3':'V3',
                         'Reporter intensity corrected 4':'V4',
                         'Reporter intensity corrected 8':'C1',
                         'Reporter intensity corrected 9':'C2',
                         'Reporter intensity corrected 10':'C3',
                         'Reporter intensity corrected 11':'C4'})
TMT2 = ['V1','V2','V3','V4','C1','C2','C3','C4']
norm_factor2 = [0.96, 1.05, 0.96, 1.13, 1, 1, 0.92, 1.07]
for i in range(len(TMT2)):
    df_TMT2[TMT2[i]] = norm_factor2[i]*df_TMT2[TMT2[i]]
    
df_TMT3 = df_aqp2[df_aqp2['Replicate'] == 'TMT3']
df_TMT3 = df_TMT3.rename(columns = {'Reporter intensity corrected 1':'C4',
                         'Reporter intensity corrected 2':'C3',
                         'Reporter intensity corrected 3':'C2',
                         'Reporter intensity corrected 4':'C1',
                         'Reporter intensity corrected 8':'V4',
                         'Reporter intensity corrected 9':'V3',
                         'Reporter intensity corrected 10':'V2',
                         'Reporter intensity corrected 11':'V1'})
TMT3 = ['C4','C3','C2','C1','V4','V3','V2','V1']
norm_factor3 = [0.99, 1.20, 0.96, 0.95, 1.21, 1.01, 1, 1]
for i in range(len(TMT3)):
    df_TMT3[TMT3[i]] = norm_factor3[i]*df_TMT3[TMT3[i]]

#single phosphosites


