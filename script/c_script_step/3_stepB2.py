import pandas as pd
import os
import json

from SolveRD.script.path_variable import *
import logging
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()

logger.info("START\tB2\n ")


#################################################################
## INPUT TO BUILD THE DF
# df_pd4_pd6_with_child contains the product 4 and 6 which have HPO, Orphacode and genes orphacode information
df_pd4_pd6_with_child = pd.read_csv(PATH_OUTPUT_PRODUCT_1_6_child, sep='\t')



# gene_cases containt for each case it related genes
gene_cases = pd.read_csv(PATH_OUTPUT_SOLVED_GENE,sep='\t')
gene_cases = gene_cases.drop(columns=['filename', 'nb_hpo'])

#################################################################
logger.info('B2\tLOAD B1 STEP')
# step B1 is needed to build B2
df_stepB1 = pd.read_csv(PATH_OUTPUT+r"stepB1.tsv", sep='\t')


##########################################################
# BUILD DF ORDO-CASES INTERACTIONS  FROM CASE CASE DATA
##########################################################



list_cases = list(df_stepB1["case"])
logger.info("B2\tnumber of interaction CO available  based on B1:{} ".format(len(list_cases)))

list_cases_json =list(set( [str(i) + ".json.json" for i in list_cases]))
logger.info("B2\tnumber of case available  based on B1: {}".format(len(list_cases_json)))

resultjsonjson = os.listdir(PATH_OUTPUT_RSLT_ALGO_ORPHA)
match_jsonresult = set(list_cases_json).intersection(set(resultjsonjson))
match_jsonresult = list(match_jsonresult)


# a set var to store orpha:codes and phenopackets  this will be use to create the graph
case_ORDO = set()


for case_B1 in match_jsonresult:

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

logger.info("B2\tonly B2 DF OUTPUT size:  : {}".format(len(df_case_ORDO)))


# add gene orpha in the B2 df
merge_allpd = pd.merge(df_case_ORDO, df_pd4_pd6_with_child, on='ORPHAcode',how='outer')





# add B1 info on B2
merge_allpdcasepheno=  pd.merge(merge_allpd, df_stepB1, on='case',how='outer')
merge_allpdcasepheno.columns = ['case','ORPHAcode','score_B2','rank_B2','symbol','ORPHAcode_child','gene_child','phenopacket','score_B1','rank_B1','gene_P','gene_C']


# merge_allpdcasepheno.columns = ['phenopacket','case','score_B1','ORPHAcode','score_B2','rank_B2','symbol','ORPHAcode_child','gene_child','gene_C','gene_P']

merge_allpdcasepheno=  merge_allpdcasepheno.dropna(subset=['score_B2'])
merge_allpdcasepheno=  merge_allpdcasepheno.drop_duplicates()

logger.info("B2\tTIME : ADD genes and B1 col ")

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
            #both empty
            pass
        elif not pd.isna(oneORPHAchildgene):
          #only orpha empty
            load = 'yes'

    # orpha child nan
    if pd.isna(oneORPHAchildgene):
            # orpha  nan
        if pd.isna(oneORPHAgene):
          # both empty
            pass
        elif not pd.isna(oneORPHAgene):
        #    only orpha child empty
            load = 'yes'

    if not pd.isna(oneORPHAgene):
        if not pd.isna(oneORPHAchildgene):
        #      both empty
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
logger.info("B2\tFILTRATION ORPHA:script child gene ")



merge_allpd_2 = merge_allpd_2.drop_duplicates()
merge_allpd_2 = merge_allpd_2.dropna(subset=['score_B1','score_B2'])





merge_allpd_2.to_csv(PATH_OUTPUT+"stepB2.tsv", sep='\t',index=False)
logger.info("B2\t EXPORT DF OUTPUT : contain all stepB2 interactions\t{}".format(len(merge_allpd_2)))


logger.info("END\tB2\n ")
