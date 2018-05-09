#!/usr/bin/python3.6
import csv
all_data_path = "/home/hermesparaqindes/Bureau/internship1/complete.txt"
all_data = open(all_data_path, "r")
all_data_csv = csv.reader(all_data, delimiter='\t')


for row in all_data_csv:
    analysis = row[4]
    blood = row[2]
    if ' RNA' in analysis and 'lymphocytes' in row[3] or ' Skin' in blood or ' Brain' in blood:
        #print(row)
        print(row[0], '\t', row[1], '\t', row[2], '\t', row[3], '\t' ,row[4], '\t' ,row[5], '\t' ,row[6], '\t' ,row[7], '\t' , row[8], '\t', row[9])
