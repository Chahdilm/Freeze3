import pandas as pd
import json
import os
from time import perf_counter
from SolveRD.script.path_variable import *

import logging
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()


logger.info("START\tB1\n ")
#################################################################
## INPUT TO BUILD THE DF

df_namefile_CC = pd.read_csv(PATH_OUTPUT_SOLVED_GENE,sep='\t')
df_namefile_CC = df_namefile_CC.drop(columns=['filename', 'nb_hpo'])

# input solved case
namefile_CC = list(set(df_namefile_CC['phenopacket']))
namefile_CC = [i+'.json.json' for i in namefile_CC]
logger.info("B1\tnumber of CO available : {}".format(len(namefile_CC)))

resultjsonjson = os.listdir(PATH_OUTPUT_RSLT_ALGO_PHENO)
match_jsonresult = set(namefile_CC).intersection(set(resultjsonjson))
match_jsonresult = list(match_jsonresult)

# get df genes for cases
df_cases_gene_excels = pd.read_csv(PATH_OUTPUT_JSON_GENE,sep='\t')


tstart = perf_counter()
# a set var to store orpha:script and phenopacket name this will be use to create the graph
case_case = list()
i = 0
while i <len(match_jsonresult):
    # open one phenopact result per loop 
    with open(PATH_OUTPUT_RSLT_ALGO_PHENO+str(match_jsonresult[i])) as file_phenopacket_result:
        one_phenopacket_result = json.load(file_phenopacket_result)

    # select only the result from a specific algo (these is 8 algo)
    for key_method in one_phenopacket_result['methods']:
        # go inside the dico which has the key name methods
        if key_method['method'] == "Resnik (symmetric)":

            for key_results in key_method['results']:
                # go inside the dico which has the key name results to get ORPHA:script info for each phenopacket
                # split part from 'P0000037.json.json' to 'P0000037,
                split_jsonfile_1 = match_jsonresult[i].split('.')
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
logger.info("B1\t : BUILD the df ")


### get genes for cases/phenopackets 
df_add_gene_case = pd.merge(df_case_case, df_namefile_CC , on='phenopacket',how='outer')
df_add_gene_case.columns = [ 'phenopacket','case','score','rank','gene_P']

df_cases_gene_excels .columns = [ 'case','gene']
df_stepB1 = pd.merge(df_add_gene_case, df_cases_gene_excels , on='case',how='outer')
df_stepB1.columns = [ 'phenopacket','case','score','rank','gene_P','gene_C',]

df_stepB1 =df_stepB1.dropna(subset=['score'])
df_stepB1 =df_stepB1.dropna(subset=['phenopacket'])

logger.info("B1\tADD genes")


df_stepB1.to_csv(PATH_OUTPUT+"stepB1.tsv", sep='\t',index=False)
logger.info("B1\tEXPORT DF B1 all {}".format(len(df_stepB1)))



logger.info("END\tB1\n ")
