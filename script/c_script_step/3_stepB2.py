import pandas as pd

import json

import datetime
from time import perf_counter

from script.path_variable import *


##############################################################################################################################

print("###############\nSTART 3_stepB2.py\n###############\n")
save_start_time =datetime.datetime.utcnow()


#################################################################
## INPUT TO BUILD THE DF
# df_pd4_pd6_with_child contains the product 4 and 6 which have HPO, Orphacode and genes orphacode information
df_pd4_pd6_with_child = pd.read_csv(PATH_OUTPUT_PRODUCT_1_6_child, sep='\t')



# gene_cases containt for each case it related genes
gene_cases = pd.read_csv(PATH_OUTPUT_SOLVED_GENE,sep='\t')
gene_cases = gene_cases.drop(columns=['filename', 'nb_hpo'])

#################################################################
print('LOAD B1 STEP')
# step B1 is needed to build B2
df_stepB1 = pd.read_csv(PATH_OUTPUT+r"stepB1.tsv", sep='\t')


print('\n####### STEP B2')

##########################################################
# BUILD DF ORDO-CASES INTERACTIONS  FROM CASE CASE DATA
##########################################################


print("B2\tTIME START build B2 : ",datetime.datetime.utcnow())


tstart = perf_counter()


list_cases = list(df_stepB1["case"])
print("B2\tnumber of interaction CO available  based on B1: ", len(list_cases))

list_cases_json =list(set( [str(i) + ".json.json" for i in list_cases]))
print("B2\tnumber of case available  based on B1: ", len(list_cases_json))

# a set var to store orpha:codes and phenopackets  this will be use to create the graph
case_ORDO = set()


for case_B1 in list_cases_json:

    # open one case_B1 per loop
    with open(PATH_OUTPUT_RSLT_ALGO_ORPHA + str(case_B1)) as file_phenopacket_result:
        one_phenopacket_result = json.load(file_phenopacket_result)

    the_case = case_B1.replace(".json.json", '')

    # select only the result from a specific algo (there is 8 algo)
    for key_method in one_phenopacket_result['methods']:
        # go inside the dico which has the key name methods
        if key_method['method'] == "Resnik (symmetric)":
            for key_results in key_method['results']:
                case_ORDO.add(( the_case,  key_results['itemid'], float(key_results['score']),key_results['rank']))

# convert the set of tuples into a dataframe
df_case_ORDO = pd.DataFrame(case_ORDO, columns=[ 'case',  'ORPHAcode', 'score_B2', 'rank_B2'])
#df_case_ORDO.to_csv(PATH_OUTPUT+"stepB2_only.tsv", sep='\t',index=False)

print("B2\tonly B2 DF OUTPUT size:  : {}".format(len(df_case_ORDO)))
tstop = perf_counter()
print("B2\tTIME THREADED BATCH : {}".format(tstop - tstart))

# add gene orpha in the B2 df
merge_allpd = pd.merge(df_case_ORDO, df_pd4_pd6_with_child, on='ORPHAcode',how='outer')





# add B1 info on B2
merge_allpdcasepheno=  pd.merge(merge_allpd, df_stepB1, on='case',how='outer')
merge_allpdcasepheno.columns = ['case','ORPHAcode','score_B2','rank_B2','symbol','ORPHAcode_child','gene_child','phenopacket','score_B1','rank_B1','gene_P','gene_C']


# merge_allpdcasepheno.columns = ['phenopacket','case','score_B1','ORPHAcode','score_B2','rank_B2','symbol','ORPHAcode_child','gene_child','gene_C','gene_P']



merge_allpdcasepheno=  merge_allpdcasepheno.dropna(subset=['score_B2'])
merge_allpdcasepheno=  merge_allpdcasepheno.drop_duplicates()

print("B2\tTIME : ADD genes and B1 col : ",datetime.datetime.utcnow())

# filtration with orpha remove orpha child when orpha parent have gene
load=""
all_interractions=set()
# remove row which containt orpha gene and orpha child gene
# keep row containing orpha gene and NO orpha gene child or orpha gene child and NO orpha gene :

dict_stepB2_df = merge_allpdcasepheno.to_dict('index')

for value in dict_stepB2_df.values():
    onepheno = value['phenopacket']
    onecase = value['case']

    oneORPHA = value['ORPHAcode']
    oneORPHAgene= value['symbol']

    onescoreB1 = value['score_B1']
    onescoreB2 = value['score_B2']

    oneORPHAchild = value['ORPHAcode_child']
    oneORPHAchildgene = value['gene_child']

    onegenecase = value['gene_C']
    onegenePHENO = value['gene_P']

    onerank = value['rank_B2']
    onerankB1= value['rank_B1']

    # orpha nan
    if pd.isna(oneORPHAgene):
        # orpha child nan
        if pd.isna(oneORPHAchildgene):
            #print('both empty')
            pass
        elif not pd.isna(oneORPHAchildgene):
          #print('only orpha empty')
            load = 'yes'

    # orpha child nan
    if pd.isna(oneORPHAchildgene):
            # orpha  nan
        if pd.isna(oneORPHAgene):
          # print('both empty')
            pass
        elif not pd.isna(oneORPHAgene):
        #    print('only orpha child empty')
            load = 'yes'

    if not pd.isna(oneORPHAgene):
        if not pd.isna(oneORPHAchildgene):
        #       print('both empty')
        #     oneORPHAchildgene =np.nan
            pass


    if load == 'yes':

        all_interractions.add((
            onepheno,
            onegenePHENO,
            onecase,
            onegenecase,
            onescoreB1,
            oneORPHA,

            onescoreB2,

            oneORPHAgene,
            oneORPHAchild,
            oneORPHAchildgene,
            onerank,
            ))
    load = ""




merge_allpd_2 = pd.DataFrame(all_interractions, columns=["phenopacket","gene_P","case","gene_C",'score_B1',"ORPHAcode",'score_B2','symbol','ORPHAcode_child','gene_child','rank'])
print("B2\tTIME : FILTRATION ORPHA:script child gene : ",datetime.datetime.utcnow())



merge_allpd_2 = merge_allpd_2.drop_duplicates()
merge_allpd_2 = merge_allpd_2.dropna(subset=['score_B1','score_B2'])






print("B2\t EXPORT DF OUTPUT : contain all stepB2 interactions\t",len(merge_allpd_2),"\nB2\tTIME : ",datetime.datetime.utcnow())

merge_allpd_2.to_csv(PATH_OUTPUT+"stepB2.tsv", sep='\t',index=False)
print("B2\t EXPORT DF OUTPUT : contain all stepB2 interactions\t",len(merge_allpd_2),"\nB2\tTIME : ",datetime.datetime.utcnow())


print("B2\tTIME build all df : ",(tstop - tstart))


print("###############\nEND 3_stepB2.py\n###############\n")
