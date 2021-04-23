# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:56:55 2021

@author: leokt
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#group i
fileloc = "J:\Depot - dDAVP-time course - Kirby\Analysis\\Data manipulation\Final Grouping\GroupI.xlsx"
#grp = pd.read_excel(fileloc, sheet_name="I")
grp_over = pd.read_excel(fileloc, sheet_name="I.A.2.c")

#grp_time = pd.DataFrame([grp[0],grp[1],grp[2],grp[5],grp[15]])
#grp_time = grp_time.rename(columns=grp['Unnamed: 1'])

grpo_time = pd.DataFrame([grp_over[0],grp_over[1],grp_over[2],grp_over[5],grp_over[15]])
grpo_time = grpo_time.rename(columns=grp_over['Unnamed: 1'])

#pantone blue: [.6,.7,.9]
#darker blue: [0,.6,1]
#light green: [.6,.9,.7]
#light red: [.9,.6,.7]
#magenta: [.8,0,.5]

#grp_plt = plt.plot(grp_time, color = [.5,.5,.5,.05])
grp_plt = plt.plot(grpo_time, color = [.8,0,.5])
grp_plt = sns.set_style('white')
grp_plt = sns.set_style('ticks')
grp_plt = sns.despine()
plt.ylabel('$log_2$(dDAVP/vehicle)', fontsize = 18)
plt.xlabel('Time (min)', fontsize = 18)
plt.axhline(0, color = 'gray', ls='--')
plt.axvline(0, color = 'gray', ls='--')
plt.ylim(-1.5,1)
plt.yticks(fontsize=12)
plt.xticks([0,5,10,15],fontsize=12)
plt.tight_layout()
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Figures\Figure X - Individual Groups\Python figures\Group IA2c_m",
            dpi=600)

#group ii
fileloc = "J:\Depot - dDAVP-time course - Kirby\Analysis\\Data manipulation\Final Grouping\GroupII.xlsx"
#grp = pd.read_excel(fileloc, sheet_name="II")
grp_over = pd.read_excel(fileloc, sheet_name="II.A.1.b")

#grp_time = pd.DataFrame([grp[0],grp[1],grp[2],grp[5],grp[15]])
#grp_time = grp_time.rename(columns=grp['Unnamed: 0'])

grpo_time = pd.DataFrame([grp_over[0],grp_over[1],grp_over[2],grp_over[5],grp_over[15]])
grpo_time = grpo_time.rename(columns=grp_over['Unnamed: 0'])

#[0,.5,.5]
#[0,.6,0]

#grp_plt = plt.plot(grp_time, color = [.5,.5,.5,.05])
grp_plt = plt.plot(grpo_time, color = [0,.5,.8])
grp_plt = sns.set_style('white')
grp_plt = sns.set_style('ticks')
grp_plt = sns.despine()
plt.ylabel('$log_2$(dDAVP/vehicle)', fontsize = 18)
plt.xlabel('Time (min)', fontsize = 18)
plt.axhline(0, color = 'gray', ls='--')
plt.axvline(0, color = 'gray', ls='--')
plt.ylim(-1,2)
plt.yticks(fontsize=12)
plt.xticks([0,5,10,15],fontsize=12)
plt.tight_layout()
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Figures\Figure X - Individual Groups\Python figures\Group IIA1b_b",
            dpi=600)


#group iii
fileloc = "J:\Depot - dDAVP-time course - Kirby\Analysis\\Data manipulation\Final Grouping\GroupIII.xlsx"
#grp = pd.read_excel(fileloc, sheet_name="III")
grp_over = pd.read_excel(fileloc, sheet_name="III.B.2")

#grp_time = pd.DataFrame([grp[0],grp[1],grp[2],grp[5],grp[15]])
#grp_time = grp_time.rename(columns=grp['Unnamed: 0'])

grpo_time = pd.DataFrame([grp_over[0],grp_over[1],grp_over[2],grp_over[5],grp_over[15]])
grpo_time = grpo_time.rename(columns=grp_over['Unnamed: 0'])

#

#grp_plt = plt.plot(grp_time, color = [.5,.5,.5,.05])
grp_plt = plt.plot(grpo_time, color = [.8,.5,0])
grp_plt = sns.set_style('white')
grp_plt = sns.set_style('ticks')
grp_plt = sns.despine()
plt.ylabel('$log_2$(dDAVP/vehicle)', fontsize = 18)
plt.xlabel('Time (min)', fontsize = 18)
plt.axhline(0, color = 'gray', ls='--')
plt.axvline(0, color = 'gray', ls='--')
plt.ylim(-1,1)
plt.yticks([-1.0,-0.5,0,0.5,1.0],fontsize=12)
plt.yticks(fontsize=12)
plt.xticks([0,5,10,15],fontsize=12)
plt.tight_layout()
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Figures\Figure X - Individual Groups\Python figures\Group IIIB2_o",
            dpi=600)



#group iv
fileloc = "J:\Depot - dDAVP-time course - Kirby\Analysis\\Data manipulation\Final Grouping\GroupIV.xlsx"
#grp = pd.read_excel(fileloc, sheet_name="IV")
grp_over = pd.read_excel(fileloc, sheet_name="IV.B")

#grp_time = pd.DataFrame([grp[0],grp[1],grp[2],grp[5],grp[15]])
#grp_time = grp_time.rename(columns=grp['Unnamed: 0'])

grpo_time = pd.DataFrame([grp_over[0],grp_over[1],grp_over[2],grp_over[5],grp_over[15]])
grpo_time = grpo_time.rename(columns=grp_over['Unnamed: 0'])

#blue: [0,.4,.6]
#grey: [.3,.3,.3]

#grp_plt = plt.plot(grp_time, color = [.5,.5,.5,.05])
grp_plt = plt.plot(grpo_time, color = [.3,.3,.3])
grp_plt = sns.set_style('white')
grp_plt = sns.set_style('ticks')
grp_plt = sns.despine()
plt.ylabel('$log_2$(dDAVP/vehicle)', fontsize = 18)
plt.xlabel('Time (min)', fontsize = 18)
plt.axhline(0, color = 'gray', ls='--')
plt.axvline(0, color = 'gray', ls='--')
plt.ylim(-0.6,1.4)
plt.yticks([-1.0,-0.5,0,0.5,1.0],fontsize=12)
plt.yticks(fontsize=12)
plt.xticks([0,5,10,15],fontsize=12)
plt.tight_layout()
plt.savefig("J:\Depot - dDAVP-time course - Kirby\Figures\Figure X - Individual Groups\Python figures\Group IVB",
            dpi=600)

