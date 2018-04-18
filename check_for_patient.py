#!/usr/bin/python3.6
import csv
import re

lymphocytes_path = "/home/hermesparaqindes/Bureau/dbGaP-13871/files/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt/data_with_lympho"
lymphocytes_data = open(lymphocytes_path, "r")
all_lymphocytes_data = list(csv.reader(lymphocytes_data, delimiter='\t'))

sample_id = []
type_tissue = []
l = []
for row in all_lymphocytes_data:
    sample_id.append(row[1])
    type_tissue.append(row[3])

reg = r'^([\w]+-[\w]+)'

dicto_for_id_tissue = dict(zip(sample_id, type_tissue))
new_dict_id_tissue = {}
for ind_id, tis in dicto_for_id_tissue.items():
    #short_ID = ind_id[:12]
    #print(ind_id)
    short_ID = "-".join(ind_id.split("-",2)[:2])
    #print(short_ID)
    if short_ID not in new_dict_id_tissue:
        new_dict_id_tissue[short_ID]= [tis]
    else:
        new_dict_id_tissue[short_ID].append(tis)

#print(new_dict_id_tissue)
#for line in all_lymphocytes_data:
#    print(line)
important_patient = ['lymphocytes', 'Skin']
important_patient1 = ' lymphocytes'
important_patient2 = ' Not Sun'
important_patient3 = 'Skin'
#print(new_dict)

potential_patient = []
for new, new_ in new_dict_id_tissue.items():
        #print(new_)
        if important_patient2 in str(new_) and important_patient1 in str(new_): #or important_patient3 in str(new_):
            #print(new, new_)
            potential_patient.append(new)
        #if any(x in important_patient for x in str(new_)) == False:
            #print(new, new_)
        #else:
            #pass
            #print(new)
#print(potential_patient)

#print(len(potential_patient))

for ligne in all_lymphocytes_data:
    for field in ligne:
        for potential in potential_patient:
            if potential in field:
                print(ligne[0], '\t', ligne[1], '\t', ligne[2], '\t', ligne[3], '\t' ,ligne[4], '\t' ,ligne[5], '\t' ,ligne[6], '\t' ,ligne[7], '\t' , ligne[8], '\t', ligne[9])
'''
class Patient:
    def __init__(self, row, header):
        self.__dict__ = dict(zip(header,row))

    def __str__(self):
        return str(self.__dict__)

lymphocytes_path = "/home/hermesparaqindes/Bureau/dbGaP-13871/files/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt/data_with_lympho"
lymphocytes_data = open(lymphocytes_path, "r")
all_lymphocytes_data = list(csv.reader(lymphocytes_data, delimiter='\t'))
print(all_lymphocytes_data[0])
patient_instance = [Patient(i, all_lymphocytes_data[0]) for i in all_lymphocytes_data[1:]]

for patient in patient_instance:
    #print(patient.SUBJID)
    if 'lymphocytes' in patient.SMTSD:
        print(patient.SUBJID, patient.SMTS)
'''
