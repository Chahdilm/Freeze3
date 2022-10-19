import pandas as pd
import json

import datetime
from time import perf_counter
from script.path_variable import *


##############################################################################################################################

print("###############\nSTART 2_stepB1.py\n###############\n")


#################################################################

#################################################################
## INPUT TO BUILD THE DF

df_namefile_CC = pd.read_csv(PATH_OUTPUT_SOLVED_GENE,sep='\t')
df_namefile_CC = df_namefile_CC.drop(columns=['filename', 'nb_hpo'])

# input solved case
namefile_CC = list(set(df_namefile_CC['phenopacket']))
namefile_CC = [i+'.json.json' for i in namefile_CC]
print("number of CO available : ",len(namefile_CC))

# get df genes for cases
df_cases_gene_excels = pd.read_csv(PATH_OUTPUT_JSON_GENE,sep='\t')

print('\n####### STEP B1')

tstart = perf_counter()

# a set var to store orpha:script and phenopacket name this will be use to create the graph
case_case = list()
i = 0
while i <len(namefile_CC):
    # open one phenopact result per loop 
    with open(PATH_OUTPUT_RSLT_ALGO_PHENO+str(namefile_CC[i])) as file_phenopacket_result:
        one_phenopacket_result = json.load(file_phenopacket_result)

    # select only the result from a specific algo (these is 8 algo)
    for key_method in one_phenopacket_result['methods']:
        # go inside the dico which has the key name methods
        if key_method['method'] == "Resnik (symmetric)":

            for key_results in key_method['results']:
                # go inside the dico which has the key name results to get ORPHA:script info for each phenopacket
                # split part from 'P0000037.json.json' to 'P0000037,
                split_jsonfile_1 = namefile_CC[i].split('.')
                split_jsonfile_1 = split_jsonfile_1[0]

                # case
                # split part from 'PP:P0007304.json' to 'P0007304,
                split_jsonfile_2 = key_results['itemid'].split(':')
                split_jsonfile_2 = split_jsonfile_2[1].split('.')
                # ['P0007304', 'json'] that why we select the first element
                split_jsonfile_2 = split_jsonfile_2[0]
 

                case_case.append((split_jsonfile_1,split_jsonfile_2,float(key_results['score']),key_results['rank']))

    i=i+1

# convert the set of tuples into a dataframe because it's easier to handle 
df_case_case = pd.DataFrame(case_case, columns=["phenopacket",'cases','score','rank'])
print("B1\tTIME : BUILD the df : ",datetime.datetime.utcnow())


### get genes for cases/phenopackets 
df_add_gene_case = pd.merge(df_case_case, df_namefile_CC , on='phenopacket',how='outer')
df_add_gene_case.columns = [ 'phenopacket','case','score','rank','gene_P']

df_cases_gene_excels .columns = [ 'case','gene']
df_stepB1 = pd.merge(df_add_gene_case, df_cases_gene_excels , on='case',how='outer')
df_stepB1.columns = [ 'phenopacket','case','score','rank','gene_P','gene_C',]

df_stepB1 =df_stepB1.dropna(subset=['score'])
df_stepB1 =df_stepB1.dropna(subset=['phenopacket'])

print("B1\tTIME : ADD genes : ",datetime.datetime.utcnow())


df_stepB1.to_csv(PATH_OUTPUT+"stepB1.tsv", sep='\t',index=False)
print("B1\tEXPORT DF B1 all ",len(df_stepB1),'\t',datetime.datetime.utcnow())

# df_case_case['phenopacket'].to_csv(PATH_OUTPUT+"nb_case_B1.tsv", sep='\t',index=False)
print("B1\tEXPORT DF OUTPUT : contain all stepB1 interactions\nB1\tTIME : ",datetime.datetime.utcnow())


print("B1\tTIME END : ",datetime.datetime.utcnow())

print('####### END B1\n')
