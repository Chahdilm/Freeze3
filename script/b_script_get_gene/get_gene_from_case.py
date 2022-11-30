
from SolveRD.script.path_variable import *
import json
import os
import pandas as pd

import logging
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()

logger.info("START\tget_gene_from_case.py\n ")

def get_gpap_solved(input, path):
    i = 0
    list_GPAP = set()
    while i < len(input):
        # open one phenopact per loop
        with open(path + '/' + str(input[i])) as file_phenopacket_result:

            one_result = json.load(file_phenopacket_result)

            one_phenopacket_result = one_result['phenopacket']
            id_phenopacket = one_phenopacket_result['id']

            if 'resolutionStatus' in one_result.keys():
                list_GPAP.add((input[i], id_phenopacket,one_result['resolutionStatus']))

            if 'interpretation' in one_result.keys():
                one_resolution_result = one_result['interpretation']
                # get list status SOLVED/UNSOLED/UNKNOW
                list_GPAP.add((input[i], id_phenopacket, one_resolution_result['resolutionStatus']))
        i = i + 1

    df_GPAP = pd.DataFrame(list_GPAP, columns=['filename', 'phenopacket', "status"])
    return df_GPAP


##################################################################################################################################################################################


def get_gene_case_solved(input, path):
    list_case_solved_with_gene = set()
    for oneinput in input:
        """ from each phenopacket files (not the phenopacket resultats json one ) we open them and extract genes info 
        """
        # open one phenopact per loop
        with open(path + '/' + str(oneinput)) as file_phenopacket_result:
            one_result = json.load(file_phenopacket_result)

            one_phenopacket_result = one_result['phenopacket']
            id_phenopacket = one_phenopacket_result['id']

            # condition on the keys because sometimes this keys is not available on the phenopacket
            if 'genes' in one_phenopacket_result.keys():
                if one_phenopacket_result['genes']:
                    # this condition is needed because sometime the gene dico exist but value are empty
                    for onegene in one_phenopacket_result['genes']:
                        # onegene == dict because one_phenopacket_result['genes'] is a dict into dict
                        list_case_solved_with_gene.add((oneinput, onegene['symbol'], id_phenopacket))

    df_case_solved_with_gene = pd.DataFrame(list_case_solved_with_gene, columns=['filename', "status", "phenopacket"])
    return df_case_solved_with_gene


##################################################################################################################################################################################


def exclude_parent(df_input,df_input_parent):
    # excluding parent
    list_parent = list(df_input_parent['PARENT_ID'].drop_duplicates())
    df_parent = df_input[~(df_input['phenopacket'].isin(list_parent))]
    return df_parent



##################################################################################################################################################################################
# HPO

def get_hpo_from_case(input, path,list_parent):
    list_case_5_HPO = list()
    list_case_1_HPO = list()
    for oneinput in input:
        """ from each phenopacket files (not the phenopacket resultats json one ) we open them and extract genes info 
        """
        oneinput_splited = oneinput.split('.')
        if oneinput_splited[0] in list_parent:
            # it s a parent
            pass
        else:
            # open one phenopact per loop
            with open(path + '/' + str(oneinput)) as file_phenopacket_result:

                one_result = json.load(file_phenopacket_result)

                one_phenopacket_result = one_result['phenopacket']
                id_phenopacket = one_phenopacket_result['id']

                # condition on the keys because sometimes this keys is not available on the phenopacket
                if 'phenotypicFeatures' in one_phenopacket_result.keys():
                    list_hpo = one_phenopacket_result['phenotypicFeatures']
                    hpo_temp = []
                    # fill hpo temp with hpo without netaged
                    for onehpo in list_hpo:
                        try:
                            if onehpo['negated']:
                                pass
                        except :
                            hpo_temp.append(onehpo)

                    if len(hpo_temp) >= 5:
                        # condition on the keys because sometimes this keys is not available on the phenopacket
                        list_case_5_HPO.append((oneinput, id_phenopacket, len(hpo_temp)))
                    if len(hpo_temp) >= 1:
                        # condition on the keys because sometimes this keys is not available on the phenopacket
                        list_case_1_HPO.append((oneinput, id_phenopacket, len(hpo_temp)))


    df_case_1_HPO = pd.DataFrame(list_case_1_HPO, columns=['filename', "phenopacket", 'nb_hpo'])
    df_case_5_HPO = pd.DataFrame(list_case_5_HPO, columns=['filename', "phenopacket", 'nb_hpo'])

    return df_case_1_HPO, df_case_5_HPO



def get_5HPO_gene_after_curation(input, path):
    all_interactions = set()

    for oneinput in input:
        """ from each phenopacket files (not the phenopacket resultats json one ) we open them and extract genes info 
        """
        # open one phenopact per loop
        with open(path + '/' + str(oneinput)) as file_phenopacket_result:
            one_result = json.load(file_phenopacket_result)

            # HPO_list = one_result['phenotypes']
            id_phenopacket = one_result['id']
            gene_list = one_result['genes']

            for oneg in gene_list:
                all_interactions.add((id_phenopacket, oneg['symbol']))
                #print(id_phenopacket, oneg['symbol'])

    df_interaction = pd.DataFrame(all_interactions, columns=["phenopacket", 'gene'])
    return df_interaction


# get list of the phenopacket folder
json_files = os.listdir(PATH_INPUT_PRODUCT_PHENOPACKET_BRUT)
len(json_files)
all_pheno_list = []
all_pheno_list_f = []
for onefile in json_files:
    # exclusion of the result folder
    if onefile != "result":
        all_pheno_list.append(onefile)

        all_pheno_list_f.append(onefile.split('.')[0])

logger.info("get_gene_from_case.py\tNumber of phenopacket available : {}".format(len(set(all_pheno_list))))
try:
    # Remove element from list phenopacket
    all_pheno_list.remove('results')
    logger.info("get_gene_from_case.py\tElement\tRemoved ")
except :
    logger.info("get_gene_from_case.py\tElement\tAlreadyRemoved ")



df_parent_interaction = pd.read_csv(PATH_OUPUT_5HPO_NOPARENT_DF,sep='\t')

# GPAP len 22264 no doublons
df_case_GPAP = get_gpap_solved(all_pheno_list,PATH_INPUT_PRODUCT_PHENOPACKET_BRUT)


# solved with gene orphanet case solved definition
input_solved_GPAP = list(df_case_GPAP[df_case_GPAP['status'] == "SOLVED"]['filename']) # 1410 no doublons
df_case_solved_with_gene = get_gene_case_solved(input_solved_GPAP,PATH_INPUT_PRODUCT_PHENOPACKET_BRUT) # 687 no doublons


# exclude parent no doublons 668
df_case_solved_with_gene_no_parent = exclude_parent(df_case_solved_with_gene,df_parent_interaction)


# HPO solved no doublons 1HPO 618    5HPO 474
input_solved_no_parent = list(df_case_solved_with_gene_no_parent['filename'].drop_duplicates())
df_final_1_hpo,df_final_5_hpo = get_hpo_from_case(input_solved_no_parent,PATH_INPUT_PRODUCT_PHENOPACKET_BRUT,df_parent_interaction['PARENT_ID'].drop_duplicates().tolist())




logger.info('get_gene_from_case.py\tNb phenopacket tot :\t{}'.format(len(all_pheno_list)))
logger.info('get_gene_from_case.py\tNb phenopacket SOLVED GPAP :\t{}'.format(len(df_case_GPAP[df_case_GPAP['status'] == "SOLVED"].drop_duplicates())))
logger.info('get_gene_from_case.py\tNb phenopacket UNSOLVED GPAP  :\t{}'.format(len(df_case_GPAP[df_case_GPAP['status'] == "UNSOLVED"].drop_duplicates())))
logger.info('get_gene_from_case.py\tNb phenopacket UNKNOW GPAP :\t{}'.format(len(df_case_GPAP[df_case_GPAP['status'] == "UNKNOWN"].drop_duplicates())))
logger.info('get_gene_from_case.py\tNb phenopacket UNKNOW + UNSOLVED GPAP :\t{}'.format((len(all_pheno_list)-len(df_case_GPAP[df_case_GPAP['status'] == "SOLVED"].drop_duplicates()))))
logger.info('get_gene_from_case.py\tNb phenopacket SOLVED Orphanet based on SOLVED GPAP :\t{}'.format(len(df_case_solved_with_gene['phenopacket'].drop_duplicates())))
logger.info('get_gene_from_case.py\tExcluding parent :\t{}'.format(len(df_case_solved_with_gene_no_parent['phenopacket'].drop_duplicates())))
logger.info('get_gene_from_case.py\t1 HPO :\t{}'.format(len(set(df_final_1_hpo['phenopacket']))))
logger.info('get_gene_from_case.py\t5 HPO :\t{}'.format(len(set(df_final_5_hpo['phenopacket']))))




json_files = os.listdir(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION)
try:
    # Remove element from list phenopacket
    json_files.remove('results')
    logger.info("get_gene_from_case.py\tElement\tRemoved ")
except :
    logger.info("get_gene_from_case.py\tElement\tAlreadyRemoved ")

df_5HPO_gene_afterC = get_5HPO_gene_after_curation(json_files,PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION)


df_final_5_hpo_afterC =df_final_5_hpo[(df_final_5_hpo['phenopacket'].isin(df_5HPO_gene_afterC['phenopacket'].tolist())) ]

df_final_5_hpo_gene_afterC = pd.merge(df_final_5_hpo_afterC, df_5HPO_gene_afterC, on='phenopacket',how='outer')
df_final_5_hpo_gene_afterC = df_final_5_hpo_gene_afterC.dropna()




logger.info('get_gene_from_case.py\tAFTER CURATION Nb phenopacket 5HPO gene on json :\t{}'.format(len(set(df_5HPO_gene_afterC['phenopacket']))))
logger.info('get_gene_from_case.py\tAFTER CURATION Nb phenopacket 5HPO gene on json SOLVED case :\t{}'.format(len(set(df_final_5_hpo_gene_afterC['phenopacket']))))



logger.info('get_gene_from_case.py\tExport df solved gene AND df json gene')
df_5HPO_gene_afterC.to_csv(PATH_OUTPUT_JSON_GENE,sep='\t',index=False)
df_final_5_hpo_gene_afterC.to_csv(PATH_OUTPUT_SOLVED_GENE,sep='\t',index=False)


logger.info("END\tget_gene_from_case.py\n ")

