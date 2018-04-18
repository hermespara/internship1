#!/usr/bin/python3.6
import csv

class Full_data:
    def __init__(self,row, header):
        self.__dict__ = dict(zip(header, row))

    def __str__(self):
        return str(self.__dict__)


full_data_path = "/home/hermesparaqindes/Bureau/dbGaP-13871/files/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt/not_sun_exposed"
full_data_open = open(full_data_path, "r")
full_data = list(csv.reader(full_data_open, delimiter ='\t'))
full_data_instance = [Full_data(i, full_data[0]) for i in full_data[1:]]





class Sra:
    def __init__(self, row, header):
        self.__dict__ = dict(zip(header, row))

    def __str__(self):
        return str(self.__dict__)


sra_data_path = "/home/hermesparaqindes/Bureau/dbGaP-13871/files/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt/SraRunTable.txt"
sra_data_open = open(sra_data_path, "r")
sra_data = list(csv.reader(sra_data_open, delimiter = '\t'))
sra_data_instance = [Sra(i, sra_data[0]) for i in sra_data[1:]]
shor_list = []
for sra in sra_data_instance:
    #print(sra.Sample_Name)
    for data in full_data_instance:
        shor = data.SAMPID.strip()
        if shor in sra.Sample_Name:
            #print(sra.Sample_Name, data.SAMPID, data.SUBJID, sra.Run)
            print(data.SUBJID, '\t', data.SAMPID,'\t', data.SMTS,'\t', data.SMTSD,'\t', data.SMNABTCHT,'\t', data.SMTSISCH,'\t', data.AGE,'\t', data.COHORT,'\t', data.RACE, '\t', data.SEX, '\t', sra.Run)
        '''
        if shor not in sra.Sample_Name:
            shor_list.append(shor)
print(len(set(shor_list)))
        '''
