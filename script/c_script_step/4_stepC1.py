
import pandas as pd
import os
import json
import numpy as np

import datetime
from time import perf_counter
from SolveRD.script.path_variable import *

import logging
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()
logger.info("START\tC1\n ")

def read_osfiles(path):
    """  store files names in a list"""
    os_files = os.listdir(path)
    len(os_files)
    # filter and keep only files names which have ".json.json" in their name because it containt the json results from the solve-rd algo
    list_os_files=[]
    for filename in os_files:
        if ".json.json" in filename:
            list_os_files.append(filename)
    
    return list_os_files

##############################################################################################################################



#### ORPHA gene
# df_pd4_pd6_with_child contains the product 4 and 6 which have HPO, Orphacode and genes orphacode information
df_pd4_pd6_with_child = pd.read_csv(PATH_OUTPUT_PRODUCT_1_6_child, sep='\t')

#################################################################
## INPUT TO BUILD THE DF
#### CASES gene
#  df_cases_gene_excels contains for each case its gene 
df_cases_gene_excels = pd.read_csv(PATH_OUTPUT_SOLVED_GENE, sep='\t')
df_cases_gene_excels = df_cases_gene_excels.drop(columns=['filename', 'nb_hpo'])

# get all phenopacket json name
namefile_CC = read_osfiles(PATH_OUTPUT_RSLT_ALGO_ORPHA)
logger.info("C1\tnumber of CO available : {} ".format(len(namefile_CC)))


#################################################################

# step A1 is needed to build C1 because for C1 cases have to be linked to the same gene as orphacode
stepA1_all = pd.read_csv(PATH_OUTPUT+r"stepA1.tsv", sep='\t')

# input of C1 are the case from A1
listephenoJS= list(set(list(stepA1_all['phenopacket'])))

logger.info("C1\tnb cases in the netwok based on A1 : {}".format(len(listephenoJS)))

#################################################################


case_ORDO_fromA1 = set()


# COI == case of interest 
for COI in listephenoJS:
    # here  for each case from A1 look for their orphacode and search for

    # extract the case of interest COI from the dataframe A1
    stepA1 = stepA1_all[stepA1_all["phenopacket"]==COI]

    # extract all orphacode and case related to the COI
    # put three list one parent only one child only and the last one exclude orpha which dont have gene but it parent or child have one
    A1_ORPHA_parent =list(stepA1['ORPHAcode'].drop_duplicates().dropna())
    A1_ORPHA_child =list(stepA1['ORPHAcode_child'].drop_duplicates().dropna())
    # filter each orpha their must have the same gene as the phenopacket
    A1_ORPHA_valid = stepA1[ (stepA1['gene_P'] == stepA1['gene_child'])  ]['ORPHAcode_child'].tolist()
    A1_ORPHA_valid = A1_ORPHA_valid + stepA1[(stepA1['gene_P'] == stepA1['symbol'] ) ]['ORPHAcode'].tolist()


    # for each case result of the algo solve-rd jsonjson 
    for case_from_jsonfile in namefile_CC:
        # open one cases per loop 
        with open(PATH_OUTPUT_RSLT_ALGO_ORPHA+str(case_from_jsonfile)) as file_phenopacket_result:
            one_phenopacket_result = json.load(file_phenopacket_result)

        # select only the result from a specific algo (these is 8 algo)
        for key_method in one_phenopacket_result['methods']:
            # go inside the dico which has the key name methods
            if key_method['method'] == "Resnik (symmetric)":

                for key_results in key_method['results']:
                    # the orpha stored in  key_results['itemid'] variable must be in the list orpha_valid because it conaint all orpha for the case on A1
                    if key_results['itemid'] in A1_ORPHA_valid:
                        # if it' s an orpha parent
                        if key_results['itemid'] in A1_ORPHA_parent:
                            # get the gene based on the key_results['itemid']  identified as orpha parent
                            gene_orphaP = stepA1[stepA1['ORPHAcode']==key_results['itemid']]['symbol']
                            gene_orphaP = gene_orphaP.iloc[0]

                            case_ORDO_fromA1.add((COI,stepA1['gene_P'].values[0],
                                                  key_results['itemid'],gene_orphaP,np.nan,np.nan,
                                                  stepA1['score'].values[0],case_from_jsonfile.replace('.json.json',''),float(key_results['score']),key_results['rank']))
                        # if it' s an orpha child no elif becasue orpha can be an orpha parent AND an orpha child
                        if key_results['itemid'] in A1_ORPHA_child:
                            # get the parent based on the key_results['itemid']  identified as orpha child
                            it_parent = stepA1[ (stepA1['ORPHAcode_child'] == key_results['itemid'] )]['ORPHAcode'].values[0]
                            # get its gene
                            gene_it_parent= stepA1[stepA1['ORPHAcode_child']==key_results['itemid']]['symbol']
                            gene_it_parent = gene_it_parent.iloc[0]
                            # get the key_results['itemid'] gene based on A1
                            gene_orphaC = stepA1[stepA1['ORPHAcode_child']==key_results['itemid']]['gene_child']
                            gene_orphaC = gene_orphaC.iloc[0]

                            case_ORDO_fromA1.add((COI,stepA1['gene_P'].values[0],
                                                  it_parent,gene_it_parent,key_results['itemid'],gene_orphaC,
                                                  stepA1['score'].values[0],case_from_jsonfile.replace('.json.json',''),float(key_results['score']),key_results['rank']))


# convert the set of tuples into a dataframe
df_case_ORDO_fromA1 = pd.DataFrame(case_ORDO_fromA1, columns=["phenopacket",'gene_P','ORPHAcode','symbol','ORPHAcode_child','gene_child','score_A1',"case",'score_C1','rank'])

df_cases_gene_excels.columns = ['case', 'gene_C']
df_stepC1 = pd.merge(df_case_ORDO_fromA1, df_cases_gene_excels , on='case',how='outer')


# remove information which have nothing in the score col which mean that algo solve-rd didn't run
df_stepC1 = df_stepC1.dropna(subset=['score_A1','score_C1'])

df_stepC1 = df_stepC1.drop_duplicates()

df_stepC1.to_csv(PATH_OUTPUT+"stepC1.tsv", sep='\t',index=False)
logger.info("C1\tEXPORT ALL C1 \t{}".format(len(df_stepC1)))




logger.info("END\tC1\n ")

