# -*- coding: utf-8 -*-
"""
Created on Tue May 11 23:02:42 2021

@author: leokt
"""

#create a list of random peptides
import pandas as pd
import numpy as np

AA_noJ = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
randpep_list = list(range(12))
for r in range(len(randpep_list)):
    rand_AA = list(range(1000))
    for w in range(len(rand_AA)):
        base = np.random.choice(AA_noJ).tolist()
        rand_AA[w] = base
    randpep_list[r] = rand_AA
randpep = np.transpose(np.asarray(randpep_list))
randpep = pd.DataFrame(randpep)
ns = (randpep.groupby([0]).size()/randpep.groupby([0]).size().sum()).to_dict()
nf = (randpep.groupby([1]).size()/randpep.groupby([1]).size().sum()).to_dict()
nfr = (randpep.groupby([2]).size()/randpep.groupby([2]).size().sum()).to_dict()
nt = (randpep.groupby([3]).size()/randpep.groupby([3]).size().sum()).to_dict()
ntw = (randpep.groupby([4]).size()/randpep.groupby([4]).size().sum()).to_dict()
no = (randpep.groupby([5]).size()/randpep.groupby([5]).size().sum()).to_dict()
o = (randpep.groupby([6]).size()/randpep.groupby([6]).size().sum()).to_dict()
tw = (randpep.groupby([7]).size()/randpep.groupby([7]).size().sum()).to_dict()
t = (randpep.groupby([8]).size()/randpep.groupby([8]).size().sum()).to_dict()
fr = (randpep.groupby([9]).size()/randpep.groupby([9]).size().sum()).to_dict()
f = (randpep.groupby([10]).size()/randpep.groupby([10]).size().sum()).to_dict()
s = (randpep.groupby([11]).size()/randpep.groupby([11]).size().sum()).to_dict()
freqrand = [ns,nf,nfr,nt,ntw,no,o,tw,t,fr,f,s]
for j in freqrand: 
    for k in range(len(AA_noJ)):
        if AA_noJ[k] not in j:
            j[AA_noJ[k]] = 0  
M = 20
#columns
N = 12
init = [ [ 0 for i in range(M) ] for j in range(N) ]
for l in range(M): 
    for m in range(N):
        init[m][l] = freqrand[m][AA_noJ[l]]      
randpep_array =  np.transpose(np.asarray(init))

import sugiyama_data

IC_randpep = randpep_array*np.log2(randpep_array/sugiyama_data.IMCDbkgd_array)
ravel_randpep = np.nan_to_num(np.ravel(IC_randpep))
corr_randpep = list(range(len(sugiyama_data.ravel_kinase)))
for y in range(len(sugiyama_data.ravel_kinase)):
    corr_randpep[y] = np.corrcoef(ravel_randpep,sugiyama_data.ravel_kinase[y])[0][1]