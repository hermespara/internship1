#!/usr/bin/python3.6
import csv

final_path = "/home/hermesparaqindes/final_results.csv"

final_data = open(final_path, "r")
all_data = list(csv.reader(final_data, delimiter=','))

#print(all_data[10])
b = []
#if all_data[10]
for line in all_data:
    if line[10] == '1': #or line[10] == '2':
        a = (line[11], line[1], line[3], line[10])
        b.append(a)
#print(b)
ind = []
list_ts = []
sam = []
for element in b[:20]:
    ind.append(element[1])
    print(element [0] + element[2] + '_'+element [0], element[1])
    tss = element[2].strip()
    list_ts.append(tss)
    sm = element[0].strip()
    sam.append(sm)

#print(set(ind))
print(*list_ts, sep = " ")

SAMPLES = []
# sra_id :  path to the file that contains all the sra run numbers listed in a column.
sra_id = '/home/hermesparaqindes/Bureau/internship1/sra_code.txt'

# open the sra_id and append to the SAMPLES list all the sra run ID.
with open(sra_id, "r") as code_file:
    for line in code_file:
        sra = line.strip()
        SAMPLES.append(sra)
# print the list of sra run ID
#print(SAMPLES)

set(sam).intersection(SAMPLES)
