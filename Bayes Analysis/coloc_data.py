# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:43:05 2021

@author: leokt
"""

# =============================================================================
# Data sourced from:
# Yang CR, Tongyoo P, Emamian M, Sandoval PC, Raghuram V, Knepper MA.
# Deep proteomic profiling of vasopressin-sensitive collecting duct 
# cells. I. Virtual Western blots and molecular weight distributions. 
# Am J Physiol Cell Physiol. 2015 Dec 15;309(12):C785-98. 
# doi: 10.1152/ajpcell.00213.2015
# =============================================================================


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

fileloc = "/Users/kirbyleo/Box Sync/Depot - dDAVP-time course - Kirby/BayesAnalysis/colocalization mapping V2.xlsx"
cluster_coloc = pd.read_excel(fileloc, sheet_name = "cluster distr")

cluster_coloc = cluster_coloc.melt(id_vars=["Gene Symbol","Cluster"],
                   var_name="Fractions",
                   value_name="Abundance level")


g = sns.FacetGrid(cluster_coloc, col="Cluster", col_wrap=7,
                  margin_titles=True)
g.map(sns.barplot, "Fractions", "Abundance level", color=".3",
      order=['1K','4K','17K','200Kp','200Ks'], ci=95)


plt.tight_layout()
plt.savefig("/Users/kirbyleo/Box Sync/Depot - dDAVP-time course - Kirby/Figures/colocalization mapping clusters v3",
            dpi=600)

#create dot-score

kinase_coloc = pd.read_excel(fileloc, sheet_name = "kinase distr")
cluster_coloc = pd.read_excel(fileloc, sheet_name = "cluster distr")

cluster_means = cluster_coloc.groupby(["Cluster"]).mean()
cluster_names = cluster_means.index
cluster_means = cluster_means.to_numpy()

kinase_names = kinase_coloc['Kinases'].to_numpy()
kinase_coloc = kinase_coloc.drop(['Kinases'],axis=1)
kinase_coloc = kinase_coloc.to_numpy()

from bayesmodules import dotscores
colocdot_df = dotscores(kinase_coloc, cluster_means, kinase_names, cluster_names)

#save file
fileloc4 = "J:\Depot - dDAVP-time course - Kirby\BayesAnalysis\\colocalization mapping V2.xlsx"
from openpyxl import load_workbook
book = load_workbook(fileloc4)
writer = pd.ExcelWriter(fileloc4, engine = 'openpyxl')
writer.book = book
colocdot_df.to_excel(writer, sheet_name = 'coloc dot ranking V2')
writer.save()
writer.close()






