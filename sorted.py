#!/usr/bin/python3.6
import csv
import operator
reader = csv.reader(open("/home/hermesparaqindes/Bureau/internship1/sra_for_not_sun_exposed", newline=''), delimiter = '\t')
sorted_list = sorted(reader, key = operator.itemgetter(0))
print(sorted_list)

with open('sorted_for_not_sun_exposed.csv', 'w') as f:
	fieldnames = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
	writer = csv.writer(f)
	for line in sorted_list:
		writer.writerow(line)
'''
with open("/home/hermesparaqindes/Bureau/dbGaP-13871/files/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt/data_with_sra_code") as csvfile:
    spamreader = csv.DictReader(csvfile, delimiter="\t")
    sortedlist = sorted(spamreader, key=lambda row:(row[0]), reverse=False)


with open('sorted.csv', 'w') as f:
    #fieldnames = ['column_1', 'column_2', column_3]
    writer = csv.DictWriter(f) #fieldnames=fieldnames)
    writer.writeheader()
    for row in sortedlist:
        writer.writerow(row)
'''
