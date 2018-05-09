#!/usr/bin/python3.6
import csv
import re

class Prelevement:
    def __init__(self, row, header):
        self.__dict__ = dict(zip(header, row))
        #self.the_id = the_id
    def __str__(self):
        return str(self.__dict__)

path = "/home/hermesparaqindes/Bureau/ncbi/dbGaP-13871/files/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt/data"
open_file = open(path, "r")
data = list(csv.reader(open_file, delimiter ='\t'))
#data = list(csv.reader(open("/home/hermesparaqindes/Bureau/dbGaP-13871/files/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt/try")), delimiter = '\t')
instance = [Prelevement(i, data[0]) for i in data[1:]]
#instance1 = {a[1]:Prelevement(a,data[0], "prelevements_{}".format(i+1)) for i, a in enumerate(data[1:])}
#print(data[1])

#for i in instance:
    #print(i.SAMPID, i.SMTS, i.SMTSD, i.SMNABTCHT, i.SMTSISCH)
    #print(i)
class Individu:
    def __init__(self, row, header):
        self.__dict__ = dict(zip(header,row))

    def __str__(self):
        return str(self.__dict__)

ind_path = "/home/hermesparaqindes/Bureau/ncbi/dbGaP-13871/files/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt/data"
open_ind = open(ind_path, "r")
ind = list(csv.reader(open_ind, delimiter = '\t'))
ind_insta = [Individu(i, ind[0]) for i in ind[1:]]
d = []

for i in ind_insta:
    for a in instance:
        reg = r'^([\w]+-[\w]+)'
        short_id = re.match(reg,a.SAMPID).group()
        #print(short_id)
        if short_id in i.SUBJID:

            print(i.SUBJID, '\t', a.SAMPID,'\t', a.SMTS,'\t', a.SMTSD,'\t', a.SMNABTCHT,'\t', a.SMTSISCH,'\t', i.AGE,'\t', i.COHORT,'\t', i.RACE, '\t', i.SEX)
        
'''
list_short = []
for a in instance:
    short_id = a.SAMPID[:10]
    list_short.append(short_id)

print(len(set(list_short)))
'''
