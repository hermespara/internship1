#!/usr/bin/python3.6
import csv
import pandas as pd

final_path = "/home/hermesparaqindes/final_results.csv"
final_data = open(final_path, "r")
all_data = list(csv.reader(final_data, delimiter=','))

L1 = []
L2 = []

for row in all_data:
    tiss = row[3].strip()
    samID = row[2].strip()
    sraID = row[-1].strip()
    #print(tiss, samID, sraID)
    L3 = [samID, sraID,tiss]
    L1.append(L3)
#print(L1)

skin_dict  = {d[0]: d[1:] for d in L1}

new_dict = {}
for k,v in skin_dict.items():
    sID = "-".join(k.split("-",2)[:2])
    if sID not in new_dict:
        new_dict[sID]=[v]
    else:
        new_dict[sID].append(v)
print(len(new_dict.keys()))
for key,val in new_dict.items():
    if len(val) ==3:
        del(val[2])
el1 = new_dict.pop('GTEX-ZQG8')
print(len(new_dict.keys()))
el2 = new_dict.pop('GTEX-ZV7C')
print(len(new_dict.keys()))
new_dict.pop('SAMPID')
for i,k in new_dict.items():
    for n in k:
        if "Skin" in n:
            print("--p /home/hermesparaqindes/GTEx/All_data_analysis/IRF_finder/"+n[0]+"/IRFinder-IR-nondir.txt", end=' ', flush=True)
for num in range(1,72):
    print("--pID skin_" +str(num), end=' ')
for i,k in new_dict.items():
    for n in k:
        if "Blood" in n:
            print("--c /home/hermesparaqindes/GTEx/All_data_analysis/IRF_finder/"+n[0]+"/IRFinder-IR-nondir.txt", end=' ', flush=True)
