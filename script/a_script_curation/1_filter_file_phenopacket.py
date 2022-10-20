
import pandas as pd
import os
import json

from script.path_variable import *

import logging
logging.basicConfig(filename='example.log', filemode='w', level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()

logger.debug("{} - Reading {} rows in {} after skipping {} rows...".format())


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
            # print('it s a parent\t',oneinput,oneinput_splited[0])
            pass
        else:
            # open one phenopact per loop
            with open(path + '\\' + str(oneinput)) as file_phenopacket_result:

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





################################################################################################################################################################################
################################################################################################################################################################################

parent_interaction = set()
ped_files = os.listdir(PATH_INPUT_PED)
for pedf in ped_files:
    with open(PATH_INPUT_PED+pedf) as file_p:
        file_p_list = file_p.readlines()
        for oneline_ped in file_p_list:
            if '\t' in oneline_ped:
                oneline_ped_t = oneline_ped.split('\t')
            else:
                oneline_ped_t = oneline_ped.split(' ')

            individualID = oneline_ped_t[1]  # individualID
            paternalID = oneline_ped_t[2] # pere
            maternalID = oneline_ped_t[3]  # mere

            parent_interaction.add((oneline_ped_t[0],individualID,paternalID,'father'))
            parent_interaction.add((oneline_ped_t[0],individualID,maternalID,'mother'))
        #print('\n NEXT FILE')

df_parent_interaction = pd.DataFrame(parent_interaction, columns=["FAM_ID",'IDV_ID','PARENT_ID','TYPE'])
df_parent_interaction2 = pd.DataFrame(parent_interaction, columns=["FAM_ID",'IDV_ID','PARENT_ID','TYPE'])

df_parent_interaction =df_parent_interaction[(df_parent_interaction['PARENT_ID'].str.contains("P0")) ]
df_parent_interaction =df_parent_interaction[(df_parent_interaction['IDV_ID'].str.contains("P0")) ]
df_parent_interaction=df_parent_interaction[df_parent_interaction.PARENT_ID.apply(lambda x: len(str(x))==8)]
df_parent_interaction=df_parent_interaction[df_parent_interaction.IDV_ID.apply(lambda x: len(str(x))==8)]

# df_parent_interaction[df_parent_interaction['PARENT_ID']=="P0012322-3"]
len(set(df_parent_interaction['PARENT_ID']))
# len(set(df_parent_interaction['FAM_ID']))

# df_parent_interaction[df_parent_interaction['IDV_ID']=="P0012322"]
# df_parent_interaction[df_parent_interaction['IDV_ID']=="P0012298"]


# df_parent_interaction[df_parent_interaction['PARENT_ID']=="P0011305"]
# df_parent_interaction2[df_parent_interaction2['PARENT_ID']=="P0011305"]


##################################################################################################################################################################################
##################################################################################################################################################################################

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

print("1_filter_file_phenopacket.py\tNumber of phenopacket available : ",len(set(all_pheno_list)))
try:
    # Remove element from list phenopacket
    all_pheno_list.remove('results')
    print("1_filter_file_phenopacket.py\tElement\tRemoved ")
except :
    print("1_filter_file_phenopacket.py\tElement\tAlreadyRemoved ")



##################################################################################################################################################################################

# tot phenopacket Exlcuding parent filename and no filename
tot_no_parent=[]
for onep in all_pheno_list_f:
    if onep not in df_parent_interaction['PARENT_ID'].tolist():
        tot_no_parent.append(onep)



# HPO tot phenopacket
# df_all_pheno_1_hpo,df_all_pheno_5_hpo = get_hpo_from_case(all_pheno_list,PATH_INPUT_PRODUCT_PHENOPACKET_BRUT)
df_all_pheno_1_hpo_no_parent,df_all_pheno_5_hpo_no_parent = get_hpo_from_case(all_pheno_list,PATH_INPUT_PRODUCT_PHENOPACKET_BRUT,df_parent_interaction['PARENT_ID'].drop_duplicates().tolist())

print(
'1_filter_file_phenopacket.py\tNb phenopacket tot Excluding parent :\t',len(tot_no_parent),'\n',
'1_filter_file_phenopacket.py\tNb phenopacket tot 5 HPO :\t',len(set(df_all_pheno_5_hpo_no_parent['phenopacket'])),'\n',
'1_filter_file_phenopacket.py\tNb phenopacket tot 1 HPO :\t',len(set(df_all_pheno_1_hpo_no_parent['phenopacket'])),'\n',
)

####################################################################################################




list_valid = df_all_pheno_5_hpo_no_parent['filename'].tolist()
len(list_valid)
for file_valid in list_valid:
    command = 'copy  ' + PATH_INPUT_PRODUCT_PHENOPACKET_BRUT+str(file_valid) +'  '+ PATH_OUPUT_5HPO_NOPARENT+str(file_valid)
    os.system(command)
print("1_filter_file_phenopacket.py\tcopy phenopacket keep only no parent 5HPO ")


# export tsv parent
df_parent_interaction.to_csv(PATH_OUPUT_5HPO_NOPARENT_DF,sep='\t',index=False)
print("1_filter_file_phenopacket.py\tExport DF parent ")

# vari = "P0003865"
#
# df_all_pheno_5_hpo_no_parent[df_all_pheno_5_hpo_no_parent['phenopacket']==vari]
# df_parent_interaction[df_parent_interaction['PARENT_ID']==vari]
# df_parent_interaction2[df_parent_interaction2['PARENT_ID']==vari]
#
#
#
# json_files_oh = os.listdir(r"C:\Users\mchahdil\Documents\Freeze3\input\RunSolveRD\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_RunSolveRD")
# json_files_mc = os.listdir(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION)
#
# json2=[]
# for onefile in json_files_mc:
#     # exclusion of the result folder
#     of = onefile.split('.')
#     of2 = of[0]+'.json'
#     json2.append(of2)
# set(json2).intersection(set(json_files_oh))
# set(json2).difference(set(json_files_oh)) # present mc mais pas oh len = 46
# tofile = set(json_files_oh).difference(set(json2)) # present oh mais pas mc len = 1202
#
#
# f = open(r"C:\Users\mchahdil\Documents\phenopacket_difference.txt",'w')
# for one in list(tofile):
#     f.write(str(one))
#     f.write('\n')
#
# f.close()
#
# dfdf = get_gene_case_solved(json_files,r"C:\Users\mchahdil\Documents\Freeze3\output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete")


