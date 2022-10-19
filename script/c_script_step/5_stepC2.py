import pandas as pd
import json


import datetime
from time import perf_counter

from script.path_variable import *





##############################################################################################################################
print("###############\nSTART 5_stepC2.py\n###############\n")
save_start_time =datetime.datetime.utcnow()

print('\n####### STEP C2')

#################################################################
## INPUT TO BUILD THE DF


#  df_excels contains for each case its gene
df_excels = pd.read_csv(PATH_OUTPUT_SOLVED_GENE, sep='\t')
df_excels = df_excels.drop(columns=['filename', 'nb_hpo'])
df_excels.columns = [ 'phenopacket','gene_P']


# df_pd4_pd6_with_child contains the product 4 and 6 which have HPO, Orphacode and genes orphacode information
df_pd4_pd6_with_child = pd.read_csv(PATH_OUTPUT_PRODUCT_1_6_child, sep='\t')
 
# import C1
stepC1 = pd.read_csv(PATH_OUTPUT+r"stepC1.tsv", sep='\t')
case_C1 =list( set(stepC1['case'].tolist()))
#get all cases from C1 -> it's the input

print("C2\tnumber of CO available  : ",len(case_C1))


#################################################################


case_ORDO = list()

# time start
tstart = perf_counter()
for COI in case_C1:
    # here for each C1 case we search the 50 closest orpha thanks to the algorithm

    with open(PATH_OUTPUT_RSLT_ALGO_ORPHA + str(COI) + '.json.json') as file_phenopacket_result:
        one_phenopacket_result = json.load(file_phenopacket_result)

    # select only the result from a specific algo (these is 8 algo)
    for key_method in one_phenopacket_result['methods']:
        # go inside the dico which has the key name methods
        if key_method['method'] == "Resnik (symmetric)":

            for key_results in key_method['results']:
                # go inside the dico which has the key name results to get ORPHA:script info for each phenopacket

                # split part from 'P0000037.json.json' to 'P0000037,
                split_jsonfile = COI.split('.')
                split_jsonfile = split_jsonfile[0]


                case_ORDO.append((split_jsonfile, key_results['itemid'], float(key_results['score']), key_results['rank']))


# print("C2\tnumber of interactions :",  len(case_ORDO))
df_case_ORDO = pd.DataFrame(case_ORDO, columns=["case", 'ORPHAcode_C2', 'score_C2', 'rank_C2'])

# print('C2\tmerge all df step C2 and df_pd4_pd6_with_child')
# select some cols and rename them
minidf_pd4_pd6 = df_pd4_pd6_with_child[["ORPHAcode", "symbol", "ORPHAcode_child", "gene_child"]]
minidf_pd4_pd6.columns = ["ORPHAcode_C2", "gene_orpha_C2", "ORPHAcode_child_C2", 'gene_child_C2']

# add genes information related to orphacode
df_case_ORDO = pd.merge(df_case_ORDO, minidf_pd4_pd6, on='ORPHAcode_C2', how='outer')




df_C2 = df_case_ORDO.dropna(subset=['gene_orpha_C2','ORPHAcode_child_C2','gene_child_C2'], how='all')
df_excels.columns = [ 'case','gene_C']
df_C2 = pd.merge(df_C2, df_excels, on='case', how='outer')
df_C2 = df_C2.dropna(subset=['case'], how='all')
df_C2 = df_C2.dropna(subset=['score_C2'], how='all')

df_C2 = df_C2.drop_duplicates()
df_C2 = df_C2.drop(columns=['gene_C'])



# df_C2.to_csv(PATH_OUTPUT+"stepC2_only.tsv", sep='\t',index=False)
print("C2\tONLY C2 \tNB :  ", len(df_C2),"\tDONE")

df_C2_all = pd.merge(stepC1, df_C2, on='case', how='outer')


df_C2_all = df_C2_all.drop_duplicates()
df_C2_all = df_C2_all.dropna(subset=['case'])
df_C2_all = df_C2_all.dropna(subset=['phenopacket'])

print("C2\tEXPORT EACH STEP C2 all(merged with C1)\tNB :  ", len(df_C2_all),"\tDONE")
print("C2\tONLY C2 \tNB :  ", len(df_C2),"\tDONE")

df_C2_all.to_csv(PATH_OUTPUT+"stepC2.tsv", sep='\t',index=False)

# time end
tstop = perf_counter()    
print("C2\tTIME build all df : ",(tstop - tstart))


print("C2\tTIME END : ",datetime.datetime.utcnow())
print('####### END C2\n')



