#!/usr/bin/python3.6
import pandas as pd

def try_parser(path):
    #read_csv file with pandas
    data = pd.read_csv(path, sep='\t', index_col=False, skiprows=10)
    #create a frame
    frame = pd.DataFrame(data)
    #dbGaP variable
    dbGaP_sample_id = frame[frame.columns[0]]
    #sample_ID
    sample_ID = frame[frame.columns[1]]
    #SMTS general tissue
    smts = frame[frame.columns[12]]
    #SMNABTCHT with type_of_analysis
    smnabtch = frame[frame.columns[20]]
    #smtsd with specific of tissues
    smtsd = frame[frame.columns[14]]
    #smtsisch with time of prelevation
    smtsisch = frame[frame.columns[17]]

    #create a list with sample_ID
    list_of_sampleID = []
    for num in sample_ID:
        list_of_sampleID.append(num)
    # list of general tissue
    general_tissue = []
    for ts in smts:
        general_tissue.append(ts)
    #list of the type of analysis
    list_of_analysis = []
    for analysis in smnabtch:
        list_of_analysis.append(analysis)
    #list of specific tissue
    specific_type_tissue = []
    for specific_tissue in smtsd:
        specific_type_tissue.append(specific_tissue)
    # list of timeof prelevation
    time_of_prelevation = []
    for time in smtsisch:
        time_of_prelevation.append(time)

    #print(time_of_prelevation)

    # dictionary with sample_ID and general type of tissue
    dictionary_with_tissue = dict(zip(list_of_sampleID, general_tissue))
    # dictionary with the type of analysis
    dictionary_with_analysis = dict(zip(list_of_sampleID, list_of_analysis))
    #dictionary with the specific tissue
    dictionary_with_specific_tissue = dict(zip(list_of_sampleID, specific_type_tissue))
    #dictionary with time
    dictionary_with_time = dict(zip(list_of_sampleID, time_of_prelevation))
    '''
    dict_id_analysis = {i:list(j) for i in dictionary_with_tissue.keys() for j in zip(dictionary_with_tissue.values(),dictionary_with_analysis.values())}
    dict_tissue = {i:list(j) for i in dict_id_analysis.keys() for j in zip(dict_id_analysis.values(),dictionary_with_specific_tissue.values())}
    print(dict_tissue)
    '''
    list_dict = [dictionary_with_time, dictionary_with_tissue, dictionary_with_analysis, dictionary_with_specific_tissue]
    output_of_dict = {}
    for d in list_dict:
        for k,v in d.items():
            output_of_dict.setdefault(k, []).append(v)

# output contains {1: [2, 3, 3], 2: [2, 1]}

    output_of_dict=[{k:v} for k,v in output_of_dict.items()]
    #print(output_of_dict)
    d = {}
    for i in output_of_dict:
        for ind, val in i.items():
            short = ind[:10]
            d.setdefault(short, {})[ind] = val
    print(d)

    list_of_interest = []
    for n in d:
        fichier = open("donnees.txt", "w")
        #for i in potential_patient:
        fichier.write(n + '\n')
        #print(n)
        #for  y in d[n]:
            #print(y)
            #print(y, ':', d[n][y])
            #for i in d[n][y]:
                #list_of_interest.append(i)
            #print(d[n][y])
            #fichier = open("donnees.txt", "w")
            #for i in potential_patient:
            #fichier.write(y)
    #print(list_of_interest)


    important_tissue = ['Brain', 'Blood', 'Skin']
    #for s_id, details in output_of_dict:

try_parser("/home/hermesparaqindes/Bureau/dbGaP-13871/files/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt/try")
