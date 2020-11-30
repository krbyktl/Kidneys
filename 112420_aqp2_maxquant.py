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

#exclude total
df_aqp2 = df_aqp2[~df_aqp2['Raw file'].str.contains('Total')]

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
"""
df_TMT1['T1'] = np.log2(df_TMT1['V1']/df_TMT1['C1'])
df_TMT1['T2'] = np.log2(df_TMT1['V2']/df_TMT1['C2'])
df_TMT1['T3'] = np.log2(df_TMT1['V3']/df_TMT1['C3'])
df_TMT1['T4'] = np.log2(df_TMT1['V4']/df_TMT1['C4'])
"""

    
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
"""
df_TMT2['T1'] = np.log2(df_TMT2['V1']/df_TMT2['C1'])
df_TMT2['T2'] = np.log2(df_TMT2['V2']/df_TMT2['C2'])
df_TMT2['T3'] = np.log2(df_TMT2['V3']/df_TMT2['C3'])
df_TMT2['T4'] = np.log2(df_TMT2['V4']/df_TMT2['C4'])
"""
    
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
"""
df_TMT3['T1'] = np.log2(df_TMT3['V1']/df_TMT3['C1'])
df_TMT3['T2'] = np.log2(df_TMT3['V2']/df_TMT3['C2'])
df_TMT3['T3'] = np.log2(df_TMT3['V3']/df_TMT3['C3'])
df_TMT3['T4'] = np.log2(df_TMT3['V4']/df_TMT3['C4'])
"""    

#TMT1 analysis
df_TMT1_1 = df_TMT1[df_TMT1['Phospho (STY)'] == 1]
df_TMT1_S256 = pd.concat([df_TMT1_1[df_TMT1_1['Modified sequence'] == '_QS(Phospho (STY))VELHSPQSLPR_'],
                          df_TMT1_1[df_TMT1_1['Modified sequence'] == '_RQS(Phospho (STY))VELHSPQSLPR_'],
                          df_TMT1_1[df_TMT1_1['Modified sequence'] == '_RRQS(Phospho (STY))VELHSPQSLPR_']])
df_TMT1_S261 = pd.concat([df_TMT1_1[df_TMT1_1['Modified sequence'] == '_QSVELHS(Phospho (STY))PQSLPR_'],
                          df_TMT1_1[df_TMT1_1['Modified sequence'] == '_RQSVELHS(Phospho (STY))PQSLPR_']])
df_TMT1_S264 = df_TMT1_1[df_TMT1_1['Modified sequence'] == '_QSVELHSPQS(Phospho (STY))LPR_']

df_TMT1_2 = df_TMT1[df_TMT1['Phospho (STY)'] == 2]
df_TMT1_S256_S261 = pd.concat([df_TMT1_2[df_TMT1_2['Modified sequence'] == '_QS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_'],
                               df_TMT1_2[df_TMT1_2['Modified sequence'] == '_RQS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_'],
                               df_TMT1_2[df_TMT1_2['Modified sequence'] == '_RRQS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_']])
df_TMT1_S261_S264 = df_TMT1_2[df_TMT1_2['Modified sequence'] == '_QSVELHS(Phospho (STY))PQS(Phospho (STY))LPR_']
df_TMT1_S256_S264 = pd.concat([df_TMT1_2[df_TMT1_2['Modified sequence'] == '_RQS(Phospho (STY))VELHSPQS(Phospho (STY))LPR_'],
                               df_TMT1_2[df_TMT1_2['Modified sequence'] == '_QS(Phospho (STY))VELHSPQS(Phospho (STY))LPR_']])

df_TMT1_S256_S261_S264 = df_TMT1[df_TMT1['Phospho (STY)'] == 3]

#TMT2 analysis
df_TMT2_1 = df_TMT2[df_TMT2['Phospho (STY)'] == 1]
df_TMT2_S256 = pd.concat([df_TMT2_1[df_TMT2_1['Modified sequence'] == '_QS(Phospho (STY))VELHSPQSLPR_'],
                          df_TMT2_1[df_TMT2_1['Modified sequence'] == '_RQS(Phospho (STY))VELHSPQSLPR_'],
                          df_TMT2_1[df_TMT2_1['Modified sequence'] == '_RRQS(Phospho (STY))VELHSPQSLPR_']])
df_TMT2_S261 = pd.concat([df_TMT2_1[df_TMT2_1['Modified sequence'] == '_QSVELHS(Phospho (STY))PQSLPR_'],
                          df_TMT2_1[df_TMT2_1['Modified sequence'] == '_RQSVELHS(Phospho (STY))PQSLPR_']])
df_TMT2_S264 = df_TMT2_1[df_TMT2_1['Modified sequence'] == '_QSVELHSPQS(Phospho (STY))LPR_']

df_TMT2_2 = df_TMT2[df_TMT2['Phospho (STY)'] == 2]
df_TMT2_S256_S261 = pd.concat([df_TMT2_2[df_TMT2_2['Modified sequence'] == '_QS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_'],
                               df_TMT2_2[df_TMT2_2['Modified sequence'] == '_RQS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_'],
                               df_TMT2_2[df_TMT2_2['Modified sequence'] == '_RRQS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_']])
df_TMT2_S261_S264 = df_TMT2_2[df_TMT2_2['Modified sequence'] == '_QSVELHS(Phospho (STY))PQS(Phospho (STY))LPR_']
df_TMT2_S256_S264 = pd.concat([df_TMT2_2[df_TMT2_2['Modified sequence'] == '_RQS(Phospho (STY))VELHSPQS(Phospho (STY))LPR_'],
                               df_TMT2_2[df_TMT2_2['Modified sequence'] == '_QS(Phospho (STY))VELHSPQS(Phospho (STY))LPR_'],
                               df_TMT2_2[df_TMT2_2['Modified sequence'] == '_RRQS(Phospho (STY))VELHSPQS(Phospho (STY))LPR_']])

df_TMT2_S256_S261_S264 = df_TMT2[df_TMT2['Phospho (STY)'] == 3]

#TMT3 analysis
df_TMT3_1 = df_TMT3[df_TMT3['Phospho (STY)'] == 1]
df_TMT3_S256 = pd.concat([df_TMT3_1[df_TMT3_1['Modified sequence'] == '_QS(Phospho (STY))VELHSPQSLPR_'],
                          df_TMT3_1[df_TMT3_1['Modified sequence'] == '_RQS(Phospho (STY))VELHSPQSLPR_'],
                          df_TMT3_1[df_TMT3_1['Modified sequence'] == '_RRQS(Phospho (STY))VELHSPQSLPR_']])
df_TMT3_S261 = pd.concat([df_TMT3_1[df_TMT3_1['Modified sequence'] == '_QSVELHS(Phospho (STY))PQSLPR_'],
                          df_TMT3_1[df_TMT3_1['Modified sequence'] == '_RQSVELHS(Phospho (STY))PQSLPR_']])
df_TMT3_S264 = df_TMT3_1[df_TMT3_1['Modified sequence'] == '_QSVELHSPQS(Phospho (STY))LPR_']
#omit GLEPDTDWEER peptides

df_TMT3_2 = df_TMT3[df_TMT3['Phospho (STY)'] == 2]
df_TMT3_S256_S261 = pd.concat([df_TMT3_2[df_TMT3_2['Modified sequence'] == '_QS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_'],
                               df_TMT3_2[df_TMT3_2['Modified sequence'] == '_RQS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_'],
                               df_TMT3_2[df_TMT3_2['Modified sequence'] == '_RRQS(Phospho (STY))VELHS(Phospho (STY))PQSLPR_']])
df_TMT3_S261_S264 = df_TMT3_2[df_TMT3_2['Modified sequence'] == '_QSVELHS(Phospho (STY))PQS(Phospho (STY))LPR_']
df_TMT3_S256_S264 = pd.concat([df_TMT3_2[df_TMT3_2['Modified sequence'] == '_RQS(Phospho (STY))VELHSPQS(Phospho (STY))LPR_'],
                               df_TMT3_2[df_TMT3_2['Modified sequence'] == '_QS(Phospho (STY))VELHSPQS(Phospho (STY))LPR_']])

df_TMT3_S256_S261_S264 = df_TMT3[df_TMT3['Phospho (STY)'] == 3]

#value export comparison
S256 = pd.concat([df_TMT1_S256, df_TMT2_S256, df_TMT3_S256])
S261 = pd.concat([df_TMT1_S261, df_TMT2_S261, df_TMT3_S261])
S264 = pd.concat([df_TMT1_S264, df_TMT2_S264, df_TMT3_S264])
S256_S261 = pd.concat([df_TMT1_S256_S261, df_TMT2_S256_S261, df_TMT3_S256_S261])
S261_S264 = pd.concat([df_TMT1_S261_S264, df_TMT2_S261_S264, df_TMT3_S261_S264])
S256_S264 = pd.concat([df_TMT1_S256_S264, df_TMT2_S256_S264, df_TMT3_S256_S264])
S256_S261_S264 = pd.concat([df_TMT1_S256_S261_S264, df_TMT2_S256_S261_S264,
                           df_TMT3_S256_S261_S264])

#calculation between the replicates
final_calculations = pd.DataFrame(data = {'TMT1': [np.log2(df_TMT1_S256['V1'].sum()/df_TMT1_S256['C1'].sum()),
                                                   np.log2(df_TMT1_S256['V2'].sum()/df_TMT1_S256['C2'].sum()),
                                                   np.log2(df_TMT1_S256['V3'].sum()/df_TMT1_S256['C3'].sum()),
                                                   np.log2(df_TMT1_S256['V4'].sum()/df_TMT1_S256['C4'].sum()),
                                                   
                                                   np.log2(df_TMT1_S261['V1'].sum()/df_TMT1_S261['C1'].sum()),
                                                   np.log2(df_TMT1_S261['V2'].sum()/df_TMT1_S261['C2'].sum()),
                                                   np.log2(df_TMT1_S261['V3'].sum()/df_TMT1_S261['C3'].sum()),
                                                   np.log2(df_TMT1_S261['V4'].sum()/df_TMT1_S261['C4'].sum()),
                                                   
                                                   np.log2(df_TMT1_S264['V1'].sum()/df_TMT1_S264['C1'].sum()),
                                                   np.log2(df_TMT1_S264['V2'].sum()/df_TMT1_S264['C2'].sum()),
                                                   np.log2(df_TMT1_S264['V3'].sum()/df_TMT1_S264['C3'].sum()),
                                                   np.log2(df_TMT1_S264['V4'].sum()/df_TMT1_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT1_S256_S261['V1'].sum()/df_TMT1_S256_S261['C1'].sum()),
                                                   np.log2(df_TMT1_S256_S261['V2'].sum()/df_TMT1_S256_S261['C2'].sum()),
                                                   np.log2(df_TMT1_S256_S261['V3'].sum()/df_TMT1_S256_S261['C3'].sum()),
                                                   np.log2(df_TMT1_S256_S261['V4'].sum()/df_TMT1_S256_S261['C4'].sum()),
                                                   
                                                   np.log2(df_TMT1_S261_S264['V1'].sum()/df_TMT1_S261_S264['C1'].sum()),
                                                   np.log2(df_TMT1_S261_S264['V2'].sum()/df_TMT1_S261_S264['C2'].sum()),
                                                   np.log2(df_TMT1_S261_S264['V3'].sum()/df_TMT1_S261_S264['C3'].sum()),
                                                   np.log2(df_TMT1_S261_S264['V4'].sum()/df_TMT1_S261_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT1_S256_S264['V1'].sum()/df_TMT1_S256_S264['C1'].sum()),
                                                   np.log2(df_TMT1_S256_S264['V2'].sum()/df_TMT1_S256_S264['C2'].sum()),
                                                   np.log2(df_TMT1_S256_S264['V3'].sum()/df_TMT1_S256_S264['C3'].sum()),
                                                   np.log2(df_TMT1_S256_S264['V4'].sum()/df_TMT1_S256_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT1_S256_S261_S264['V1'].sum()/df_TMT1_S256_S261_S264['C1'].sum()),
                                                   np.log2(df_TMT1_S256_S261_S264['V2'].sum()/df_TMT1_S256_S261_S264['C2'].sum()),
                                                   np.log2(df_TMT1_S256_S261_S264['V3'].sum()/df_TMT1_S256_S261_S264['C3'].sum()),
                                                   np.log2(df_TMT1_S256_S261_S264['V4'].sum()/df_TMT1_S256_S261_S264['C4'].sum())
                                                   ],
                                          'TMT2': [np.log2(df_TMT2_S256['V1'].sum()/df_TMT2_S256['C1'].sum()),
                                                   np.log2(df_TMT2_S256['V2'].sum()/df_TMT2_S256['C2'].sum()),
                                                   np.log2(df_TMT2_S256['V3'].sum()/df_TMT2_S256['C3'].sum()),
                                                   np.log2(df_TMT2_S256['V4'].sum()/df_TMT2_S256['C4'].sum()),
                                                   
                                                   np.log2(df_TMT2_S261['V1'].sum()/df_TMT2_S261['C1'].sum()),
                                                   np.log2(df_TMT2_S261['V2'].sum()/df_TMT2_S261['C2'].sum()),
                                                   np.log2(df_TMT2_S261['V3'].sum()/df_TMT2_S261['C3'].sum()),
                                                   np.log2(df_TMT2_S261['V4'].sum()/df_TMT2_S261['C4'].sum()),
                                                   
                                                   np.log2(df_TMT2_S264['V1'].sum()/df_TMT2_S264['C1'].sum()),
                                                   np.log2(df_TMT2_S264['V2'].sum()/df_TMT2_S264['C2'].sum()),
                                                   np.log2(df_TMT2_S264['V3'].sum()/df_TMT2_S264['C3'].sum()),
                                                   np.log2(df_TMT2_S264['V4'].sum()/df_TMT2_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT2_S256_S261['V1'].sum()/df_TMT2_S256_S261['C1'].sum()),
                                                   np.log2(df_TMT2_S256_S261['V2'].sum()/df_TMT2_S256_S261['C2'].sum()),
                                                   np.log2(df_TMT2_S256_S261['V3'].sum()/df_TMT2_S256_S261['C3'].sum()),
                                                   np.log2(df_TMT2_S256_S261['V4'].sum()/df_TMT2_S256_S261['C4'].sum()),
                                                   
                                                   np.log2(df_TMT2_S261_S264['V1'].sum()/df_TMT2_S261_S264['C1'].sum()),
                                                   np.log2(df_TMT2_S261_S264['V2'].sum()/df_TMT2_S261_S264['C2'].sum()),
                                                   np.log2(df_TMT2_S261_S264['V3'].sum()/df_TMT2_S261_S264['C3'].sum()),
                                                   np.log2(df_TMT2_S261_S264['V4'].sum()/df_TMT2_S261_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT2_S256_S264['V1'].sum()/df_TMT2_S256_S264['C1'].sum()),
                                                   np.log2(df_TMT2_S256_S264['V2'].sum()/df_TMT2_S256_S264['C2'].sum()),
                                                   np.log2(df_TMT2_S256_S264['V3'].sum()/df_TMT2_S256_S264['C3'].sum()),
                                                   np.log2(df_TMT2_S256_S264['V4'].sum()/df_TMT2_S256_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT2_S256_S261_S264['V1'].sum()/df_TMT2_S256_S261_S264['C1'].sum()),
                                                   np.log2(df_TMT2_S256_S261_S264['V2'].sum()/df_TMT2_S256_S261_S264['C2'].sum()),
                                                   np.log2(df_TMT2_S256_S261_S264['V3'].sum()/df_TMT2_S256_S261_S264['C3'].sum()),
                                                   np.log2(df_TMT2_S256_S261_S264['V4'].sum()/df_TMT2_S256_S261_S264['C4'].sum())
                                                   ],
                                          'TMT3': [np.log2(df_TMT3_S256['V1'].sum()/df_TMT3_S256['C1'].sum()),
                                                   np.log2(df_TMT3_S256['V2'].sum()/df_TMT3_S256['C2'].sum()),
                                                   np.log2(df_TMT3_S256['V3'].sum()/df_TMT3_S256['C3'].sum()),
                                                   np.log2(df_TMT3_S256['V4'].sum()/df_TMT3_S256['C4'].sum()),
                                                   
                                                   np.log2(df_TMT2_S261['V1'].sum()/df_TMT2_S261['C1'].sum()),
                                                   np.log2(df_TMT2_S261['V2'].sum()/df_TMT2_S261['C2'].sum()),
                                                   np.log2(df_TMT2_S261['V3'].sum()/df_TMT2_S261['C3'].sum()),
                                                   np.log2(df_TMT2_S261['V4'].sum()/df_TMT2_S261['C4'].sum()),
                                                   
                                                   np.log2(df_TMT3_S264['V1'].sum()/df_TMT3_S264['C1'].sum()),
                                                   np.log2(df_TMT3_S264['V2'].sum()/df_TMT3_S264['C2'].sum()),
                                                   np.log2(df_TMT3_S264['V3'].sum()/df_TMT3_S264['C3'].sum()),
                                                   np.log2(df_TMT3_S264['V4'].sum()/df_TMT3_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT3_S256_S261['V1'].sum()/df_TMT3_S256_S261['C1'].sum()),
                                                   np.log2(df_TMT3_S256_S261['V2'].sum()/df_TMT3_S256_S261['C2'].sum()),
                                                   np.log2(df_TMT3_S256_S261['V3'].sum()/df_TMT3_S256_S261['C3'].sum()),
                                                   np.log2(df_TMT3_S256_S261['V4'].sum()/df_TMT3_S256_S261['C4'].sum()),
                                                   
                                                   np.log2(df_TMT3_S261_S264['V1'].sum()/df_TMT3_S261_S264['C1'].sum()),
                                                   np.log2(df_TMT3_S261_S264['V2'].sum()/df_TMT3_S261_S264['C2'].sum()),
                                                   np.log2(df_TMT3_S261_S264['V3'].sum()/df_TMT3_S261_S264['C3'].sum()),
                                                   np.log2(df_TMT3_S261_S264['V4'].sum()/df_TMT3_S261_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT3_S256_S264['V1'].sum()/df_TMT3_S256_S264['C1'].sum()),
                                                   np.log2(df_TMT3_S256_S264['V2'].sum()/df_TMT3_S256_S264['C2'].sum()),
                                                   np.log2(df_TMT3_S256_S264['V3'].sum()/df_TMT3_S256_S264['C3'].sum()),
                                                   np.log2(df_TMT3_S256_S264['V4'].sum()/df_TMT3_S256_S264['C4'].sum()),
                                                   
                                                   np.log2(df_TMT3_S256_S261_S264['V1'].sum()/df_TMT3_S256_S261_S264['C1'].sum()),
                                                   np.log2(df_TMT3_S256_S261_S264['V2'].sum()/df_TMT3_S256_S261_S264['C2'].sum()),
                                                   np.log2(df_TMT3_S256_S261_S264['V3'].sum()/df_TMT3_S256_S261_S264['C3'].sum()),
                                                   np.log2(df_TMT3_S256_S261_S264['V4'].sum()/df_TMT3_S256_S261_S264['C4'].sum())
                                                   ]})

from openpyxl import load_workbook
book = load_workbook(fileloc)
writer = pd.ExcelWriter(fileloc, engine = 'openpyxl')
writer.book = book
S256.to_excel(writer, sheet_name = 'S256')
S261.to_excel(writer, sheet_name = 'S261')
S264.to_excel(writer, sheet_name = 'S264')
S256_S261.to_excel(writer, sheet_name = 'S256_S261')
S261_S264.to_excel(writer, sheet_name = 'S261_S264')
S256_S264.to_excel(writer, sheet_name = 'S256_S264')
S256_S261_S264.to_excel(writer, sheet_name = 'S256_S261_S264')
final_calculations.to_excel(writer, sheet_name = 'calculations')

writer.save()
writer.close()
               
#plot the values after hand calculations
df_aqp2_new = pd.read_excel(fileloc, sheet_name="calculations")

fig, ax = plt.subplots(figsize=(12,10))
for key, group in df_aqp2_new.groupby('Phosphosites'):
    group.plot('Timepoint','OVERALL RATIO', yerr='SD',
                linewidth = 3, label = key, ax = ax)
sns.despine()
plt.title('Aqp2 Time course', fontsize = 24)
plt.ylabel('$log_2$(dDAVP/control)', fontsize = 18)
plt.xlabel('Time (min)', fontsize = 18)
plt.ylim(-2,1.5)
plt.legend(bbox_to_anchor=(1.35, 1), loc='upper right', fontsize=18)
ax.axhline(0, color = 'gray', ls='--')
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Scripts\Figures\\Aqp2_Phosphosites",
            bbox_inches='tight', dpi = 1200)

fig, ax = plt.subplots(figsize=(12,10))
sns.lineplot(data = df_aqp2_new, x = 'Timepoint',
             y = 'OVERALL RATIO', hue = 'Phosphosites', palette = 'mako',
             linewidth = 3)
sns.despine()
plt.title('Aqp2 Time course', fontsize = 24)
plt.ylabel('$log_2$(dDAVP/control)', fontsize = 18)
plt.xlabel('Time (min)', fontsize = 18)
plt.ylim(-2,1.5)
plt.legend(bbox_to_anchor=(1.35, 1), loc='upper right', fontsize=18)
ax.axhline(0, color = 'gray', ls='--')
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Scripts\Figures\\Aqp2_Phosphosites_nosd",
            bbox_inches='tight', dpi = 1200)











np.log2(df_TMT2_S256['V1'].sum())/(df_TMT2_S256['C1']),
log2(df_TMT3_S256['V1'].sum())/(df_TMT3_S256['C1']