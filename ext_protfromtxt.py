# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 13:47:43 2020

@author: leokt
"""
#find klf6
import pandas as pd
import time

start_time = time.time()
fileloc = "J:\Depot - dDAVP-time course - Kirby\MaxQuant output files\\msms.txt"
input_file = open(fileloc, 'r')
klf6_list = []

with input_file as inF:
    for line in inF:
        if 'Klf6' in line:
            klf6_list.append(line)

print("My program took", time.time() - start_time, "to run")
klf6_df = pd.DataFrame([sub.split("\t") for sub in klf6_list])


#test with aqp2
input_file = open(fileloc, 'r')
aqp2_list = []
start_time = time.time()
with input_file as inF:
    for line in inF:
        if 'Aqp2' in line:
            aqp2_list.append(line)

print("My program took", time.time() - start_time, "to run")
aqp2_df = pd.DataFrame([sub.split("\t") for sub in aqp2_list])

#Clip1 and Slc14a2
import pandas as pd
import time

start_time = time.time()
fileloc = "J:\Depot - dDAVP-time course - Kirby\MaxQuant output files\\msms.txt"
input_file = open(fileloc, 'r')
clip1_list = []
slc14a2_list = []


with input_file as inF:
    for line in inF:
        if 'Clip1' in line:
            clip1_list.append(line)
        if 'Slc14a2' in line:
            slc14a2_list.append(line)

print("My program took", time.time() - start_time, "seconds to run")
clip1_df = pd.DataFrame([sub.split("\t") for sub in clip1_list])
slc14a2_df = pd.DataFrame([sub.split("\t") for sub in slc14a2_list])




