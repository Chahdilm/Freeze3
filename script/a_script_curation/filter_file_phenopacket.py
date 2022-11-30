
import pandas as pd
import os
import json
import subprocess

from SolveRD.script.path_variable import *
import logging


logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()

logger.info("START\tfilter_file_phenopacket.py\n ")

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
            # phenopacket defined as parent aren't open
            pass
        else:

            with open(path + '/' + str(oneinput)) as file_phenopacket_result:

                one_result = json.load(file_phenopacket_result)
                # open json file

                one_phenopacket_result = one_result['phenopacket']
                id_phenopacket = one_phenopacket_result['id']

                # phenotypicFeatures key store hpo
                if 'phenotypicFeatures' in one_phenopacket_result.keys():
                    list_hpo = one_phenopacket_result['phenotypicFeatures']
                    hpo_temp = []
                    #  hpo_temp don't store  netaged HPO
                    #  list_hpo contains all HPO then thanks to hpo_temp list negated hpo are remove
                    for onehpo in list_hpo:
                        try:
                            if onehpo['negated']:
                                pass
                        except :
                            hpo_temp.append(onehpo)
                    # end loop contain all hpo for the phenopacket
                    if len(hpo_temp) >= 5:
                        # if the number of element in hpo_temp is hight than 5 then the phenopacket contains more than 5HPO
                        list_case_5_HPO.append((oneinput, id_phenopacket, len(hpo_temp)))
                    if len(hpo_temp) >= 1:
                        # same process for 1 HPO
                        list_case_1_HPO.append((oneinput, id_phenopacket, len(hpo_temp)))

    # convert the list of tuple into a df
    df_case_1_HPO = pd.DataFrame(list_case_1_HPO, columns=['filename', "phenopacket", 'nb_hpo'])
    df_case_5_HPO = pd.DataFrame(list_case_5_HPO, columns=['filename', "phenopacket", 'nb_hpo'])

    return df_case_1_HPO, df_case_5_HPO




################################################################################################################################################################################
################################################################################################################################################################################
# PARENT SECTION
# get all parents phenopacket thanks to ped file
parent_interaction = set()
ped_files = os.listdir(PATH_INPUT_PED)

for pedf in ped_files:
    with open(PATH_INPUT_PED+pedf) as file_p:
        # for each ped file name  we open the ped file
        file_p_list = file_p.readlines()
        for oneline_ped in file_p_list:
            # some ped file usr ' ' seperator and other '\t'.
            if '\t' in oneline_ped:
                oneline_ped_t = oneline_ped.split('\t')
            else:
                oneline_ped_t = oneline_ped.split(' ')

            individualID = oneline_ped_t[1]  # individualID
            paternalID = oneline_ped_t[2] # pere
            maternalID = oneline_ped_t[3]  # mere

            parent_interaction.add((oneline_ped_t[0],individualID,paternalID,'father'))
            parent_interaction.add((oneline_ped_t[0],individualID,maternalID,'mother'))
        # NEXT PED FILE

df_parent_interaction = pd.DataFrame(parent_interaction, columns=["FAM_ID",'IDV_ID','PARENT_ID','TYPE'])

# filtration on df parent
# keep only id which contain PO because a phenopacket id always containt PO example P0011778
df_parent_interaction =df_parent_interaction[(df_parent_interaction['PARENT_ID'].str.contains("P0")) ]
df_parent_interaction =df_parent_interaction[(df_parent_interaction['IDV_ID'].str.contains("P0")) ]
# remove all id which lenght is not 8 because phenopacket id are store into 8 caracters only.
df_parent_interaction=df_parent_interaction[df_parent_interaction.PARENT_ID.apply(lambda x: len(str(x))==8)]
df_parent_interaction=df_parent_interaction[df_parent_interaction.IDV_ID.apply(lambda x: len(str(x))==8)]




##################################################################################################################################################################################
##################################################################################################################################################################################
#store all phenopacket names on json_files liste
json_files = os.listdir(PATH_INPUT_PRODUCT_PHENOPACKET_BRUT)
len(json_files)
all_pheno_list = []
all_pheno_list_f = []
for onefile in json_files:
    # exclusion of the result folder
    if onefile != "result":
        all_pheno_list.append(onefile)
        # remove the date on the file name
        all_pheno_list_f.append(onefile.split('.')[0])
logger.info("filter_file_phenopacket.py\tNumber of phenopacket available : \t{}".format(len(set(all_pheno_list))))


try:
    # Remove element from list phenopacket
    all_pheno_list.remove('results')
    logger.info("filter_file_phenopacket.py\tElement\tRemoved")

except :
    logger.info("filter_file_phenopacket.py\tElement\tAlreadyRemoved")




##################################################################################################################################################################################
# EXCLUDING PHENOPACKET PARENT
# tot phenopacket Exlcuding parent filename
tot_no_parent=[]
for onep in all_pheno_list_f:
    if onep not in df_parent_interaction['PARENT_ID'].tolist():
        tot_no_parent.append(onep)



# HPO tot phenopacket
df_all_pheno_1_hpo_no_parent,df_all_pheno_5_hpo_no_parent = get_hpo_from_case(all_pheno_list,PATH_INPUT_PRODUCT_PHENOPACKET_BRUT,df_parent_interaction['PARENT_ID'].drop_duplicates().tolist())
logger.info("filter_file_phenopacket.py\tNb phenopacket tot Excluding parent :\t{}".format(len(tot_no_parent)))
logger.info("filter_file_phenopacket.py\tNb phenopacket tot 5 HPO :\t{}".format(len(set(df_all_pheno_5_hpo_no_parent['phenopacket']))))
logger.info("filter_file_phenopacket.py\tNb phenopacket tot 1 HPO :\t{}".format(len(set(df_all_pheno_1_hpo_no_parent['phenopacket']))))






####################################################################################################
# # get a valid list of all phenopacket filtered no parent and with 5 HPO or more
list_valid = df_all_pheno_5_hpo_no_parent['filename'].tolist()



for file_valid in list_valid:
    # # window command
    # command = 'copy  ' + PATH_INPUT_PRODUCT_PHENOPACKET_BRUT+str(file_valid) +'  '+ PATH_OUPUT_5HPO_NOPARENT+str(file_valid)
    # # Linux command
    command = 'cp  ' + PATH_INPUT_PRODUCT_PHENOPACKET_BRUT+str(file_valid) +'  '+ PATH_OUPUT_5HPO_NOPARENT+str(file_valid)
    subprocess.run(command,shell=True)
    logger.info("filter_file_phenopacket.py\tcopy phenopacket\t{}".format(file_valid))


logger.info("filter_file_phenopacket.py\tDONE copy phenopacket keep only no parent 5HPO ")


# # export tsv parent
df_parent_interaction.to_csv(PATH_OUPUT_5HPO_NOPARENT_DF,sep='\t',index=False)
logger.info("filter_file_phenopacket.py\tExport DF parent ")

logger.info("END\t1_filter_file_phenopacket\n ")



