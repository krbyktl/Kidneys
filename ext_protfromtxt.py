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

from openpyxl import load_workbook
new_loc = "J:\Depot - dDAVP-time course - Kirby\Analysis\\Klf6_Msms.xlsx"
book = load_workbook(new_loc)
writer = pd.ExcelWriter(new_loc, engine = 'openpyxl')
writer.book = book
klf6_df.to_excel(writer, sheet_name = 'Klf6_msms')


#test
input_file = open(fileloc, 'r')
aqp2_list = []
start_time = time.time()
with input_file as inF:
    for line in inF:
        if 'Aqp2' in line:
            aqp2_list.append(line)

print("My program took", time.time() - start_time, "to run")
aqp2_df = pd.DataFrame([sub.split("\t") for sub in aqp2_list])


