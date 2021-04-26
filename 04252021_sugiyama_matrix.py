# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 19:40:56 2021

@author: leokt
"""
import pandas as pd
import numpy as np

#import data and clean up
fileloc = "J:\Depot - dDAVP-time course - Kirby\BayesAnalysis\\clean sugiyama kinase substrates.xlsx"
ks_df = pd.read_excel(fileloc)

ks_df.rename(columns = {'(-6)':'negsix',
                        '(-5)':'negfive',
                        '(-4)':'negfour',
                        '(-3)':'negthree',
                        '(-2)':'negtwo',
                        '(-1)':'negone',
                        '(+1)':'posone',
                        '(+2)':'postwo',
                        '(+3)':'posthree',
                        '(+4)':'posfour',
                        '(+5)':'posfive',
                        '(+6)':'possix',}, inplace = True)

#divide dataframe per kinase
kinases = ks_df['Kinase'].unique()

AA = ['A','C','D','E','F','G','H','I','J','K','L','M','N','P','Q','R','S','T','V','W','Y']
#pos = {-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6}

kin_ind = list(range(len(kinases)))
freq_ind = list(range(len(kinases)))

for i in kin_ind:
    kin_ind[i] = (ks_df[ks_df['Kinase']==kinases[i]])

    #calculate frequencies in each position per kinase
    #initialize an empty matrix
    #rows
    M = 21
    #columns
    N= 13
    init = [ [ 0 for i in range(M) ] for j in range(N) ]
    
    ns = (kin_ind[i].groupby(['negsix']).size()/kin_ind[i].groupby(['negsix']).size().sum()).to_dict()
    nf = (kin_ind[i].groupby(['negfive']).size()/kin_ind[i].groupby(['negfive']).size().sum()).to_dict()
    nfr = (kin_ind[i].groupby(['negfour']).size()/kin_ind[i].groupby(['negfour']).size().sum()).to_dict()
    nt = (kin_ind[i].groupby(['negthree']).size()/kin_ind[i].groupby(['negthree']).size().sum()).to_dict()
    ntw = (kin_ind[i].groupby(['negtwo']).size()/kin_ind[i].groupby(['negtwo']).size().sum()).to_dict()
    no = (kin_ind[i].groupby(['negone']).size()/kin_ind[i].groupby(['negone']).size().sum()).to_dict()
    zero = (kin_ind[i].groupby([0]).size()/kin_ind[i].groupby([0]).size().sum()).to_dict()
    o = (kin_ind[i].groupby(['posone']).size()/kin_ind[i].groupby(['posone']).size().sum()).to_dict()
    tw = (kin_ind[i].groupby(['postwo']).size()/kin_ind[i].groupby(['postwo']).size().sum()).to_dict()
    t = (kin_ind[i].groupby(['posthree']).size()/kin_ind[i].groupby(['posthree']).size().sum()).to_dict()
    fr = (kin_ind[i].groupby(['posfour']).size()/kin_ind[i].groupby(['posfour']).size().sum()).to_dict()
    f = (kin_ind[i].groupby(['posfive']).size()/kin_ind[i].groupby(['posfive']).size().sum()).to_dict()
    s = (kin_ind[i].groupby(['possix']).size()/kin_ind[i].groupby(['possix']).size().sum()).to_dict()
    freq = [ns,nf,nfr,nt,ntw,no,zero,o,tw,t,fr,f,s]
    
    for j in freq: 
        for k in range(M):
            if AA[k] not in j:
                j[AA[k]] = 0

#map frequencies into a 21x13 matrix 
    for l in range(21): 
        for m in range(13):
            init[m][l] = freq[m][AA[l]]      
    print(init)
    freq_ind[i] =  np.transpose(np.asarray(init))
    
#exclude "J" and 0 position entry to convert to 20x12 matrix
freq_array_mod = list(range(len(kinases)))
for n in range(len(freq_ind)):
    mod_array = np.delete(freq_ind[n],8,0)
    mod_array = np.delete(mod_array,6,1)
    freq_array_mod[n] = mod_array

#placeholder dictionary
kinase_freq_dict = {k:v for k,v in zip(kinases,freq_array_mod)}

#introduce reference probabilities
fileloc1 = "J:\Depot - dDAVP-time course - Kirby\BayesAnalysis\\human PPSP background.xlsx"
Sbkgd_df = pd.read_excel(fileloc1, 'S-center')
Tbkgd_df = pd.read_excel(fileloc1, 'T-center')
#Ybkgd_df = pd.read_excel(fileloc1, 'Y-center')
Sbkgd_df = Sbkgd_df.drop([20,21])
Sbkgd_df = Sbkgd_df.drop(['Position:',0,-7,7], axis=1)
Tbkgd_df = Tbkgd_df.drop([20,21])
Tbkgd_df = Tbkgd_df.drop(['Position:',0,-7,7], axis=1)
Sbkgd_array = Sbkgd_df.to_numpy()
Tbkgd_array = Tbkgd_df.to_numpy()
#average S and T backgrounds
STbkgd_array = np.mean([Sbkgd_array,Tbkgd_array],axis=0)/100

#Convert observed frequencies to information content scores [matrix] per kinase
#IC(a,i)=P(a,i)*log((P(a,i)/Pref(a)),2)
IC_array = list(range(len(kinases)))
for b in IC_array:
    IC_array[b] = freq_array_mod[b]*np.log2(freq_array_mod[b]/STbkgd_array)

kinase_IC_dict = {k:v for k,v in zip(kinases,IC_array)}






