# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:43:05 2021

@author: leokt
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

fileloc = "J:\Depot - dDAVP-time course - Kirby\BayesAnalysis\\colocalization mapping.xlsx"
cluster_coloc = pd.read_excel(fileloc, sheet_name = "norm cluster distr")

cluster_coloc = cluster_coloc.melt(id_vars=["Gene Symbol","Position","Centralized 13 aa","Cluster"],
                   var_name="Fractions",
                   value_name="Normalized abundance level")


g = sns.FacetGrid(cluster_coloc, col="Cluster", col_wrap=7,
                  margin_titles=True)
g.map(sns.barplot, "Fractions", "Normalized abundance level", color=".3",
      order=['1K','4K','17K','200Kp','200Ks'], ci=95)


plt.tight_layout()
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Figures\colocalization mapping clusters",
            dpi=600)

kinase_coloc = pd.read_excel(fileloc, sheet_name = "norm kinase distr")
cluster_coloc = pd.read_excel(fileloc, sheet_name = "norm cluster distr")

cluster_means = cluster_coloc.groupby(["Cluster"]).mean()
cluster_names = cluster_means.index
cluster_means = cluster_means.to_numpy()

kinase_names = kinase_coloc['Kinases'].to_numpy()
kinase_coloc = kinase_coloc.drop(['Kinases'],axis=1)
kinase_coloc = kinase_coloc.to_numpy()

coloc_dotscores = list(range(len(cluster_means)))
for u in range(len(cluster_means)):
    kinscores = list(range(len(kinase_coloc)))
    for p in range(len(kinase_coloc)):
        rank = np.dot(kinase_coloc[p],cluster_means[u])
        kinscores[p] = rank
    coloc_dotscores[u] = kinscores
colocdot_df = pd.DataFrame(np.transpose(coloc_dotscores), columns = cluster_names, index = kinase_names.tolist())       


fileloc4 = "J:\Depot - dDAVP-time course - Kirby\BayesAnalysis\\colocalization mapping.xlsx"
from openpyxl import load_workbook
book = load_workbook(fileloc4)
writer = pd.ExcelWriter(fileloc4, engine = 'openpyxl')
writer.book = book
colocdot_df.to_excel(writer, sheet_name = 'coloc dot ranking')
writer.save()
writer.close()






