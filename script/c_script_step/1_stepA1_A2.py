import pandas as pd
import json
import os

from SolveRD.script.path_variable import *
import logging
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()

logger.info("START\tA1_A2\n ")

#################################################################

## INPUT TO BUILD THE DF

df_excels = pd.read_csv(PATH_OUTPUT_SOLVED_GENE, sep='\t')
df_excels = df_excels.drop(columns=['filename', 'nb_hpo'])

# input solved case
namefile_CO = list(set(df_excels['phenopacket']))
namefile_CO = [i+'.json.json' for i in namefile_CO]
logger.info("A1_A2\tnumber of CO available : {}".format(len(namefile_CO)))

logger.info('A2')

# intersection between list from solved case with gene and results from algo
resultjsonjson = os.listdir(PATH_OUTPUT_RSLT_ALGO_ORPHA)
match_jsonresult = set(namefile_CO).intersection(set(resultjsonjson))
match_jsonresult = list(match_jsonresult)
# a set var to store orpha:script and phenopacket name this will be use to create the graph
case_ORPHA = list()
i = 0

while i <len(match_jsonresult):
    # open one phenopact result per loop 
    with open(PATH_OUTPUT_RSLT_ALGO_ORPHA+str(match_jsonresult[i]),'r',encoding = 'utf8') as file_phenopacket_result:
        one_phenopacket_result = json.load(file_phenopacket_result,encoding = 'utf8')

    # select only the result from a specific algo (these is 8 algo)
    for key_method in one_phenopacket_result['methods']:
        # go inside the dico which has the key name methods
        #if key_method['method'] == "Cosine Weighted Similarity":
        if key_method['method'] == "Resnik (symmetric)":
            # key_method contains all information of the phenopacket
            for key_results in key_method['results']:

                # split part from 'P0000037.json.json' to 'P0000037,
                split_jsonfile = match_jsonresult[i].split('.')
                split_jsonfile = split_jsonfile[0]
                # ,key_results['itemid'] contain orpha:script, ,key_results['score'] contain score and so on
                case_ORPHA.append((split_jsonfile,key_results['itemid'],float(key_results['score']),key_results['rank']))
    i=i+1

# convert the set of tuples into a dataframe 
df_case_ORPHA = pd.DataFrame(case_ORPHA, columns=["phenopacket",'ORPHAcode','score','rank'])
logger.info("A2\tBUILD the df ")

# add genes related to case
df_add_gene_case = pd.merge(df_case_ORPHA, df_excels, on='phenopacket',how='outer')
df_add_gene_case.columns = [ 'phenopacket','ORPHAcode','score','rank','P gene']

# add genes related to orpha
df_pd_with_child = pd.read_csv(PATH_OUTPUT_PRODUCT_1_6_child, sep='\t')

# merge df which gene case and df with gene orpha:script
df_stepA = pd.merge(df_add_gene_case, df_pd_with_child, on='ORPHAcode',how='outer')
df_stepA = df_stepA.dropna(subset=['score'])
df_stepA = df_stepA[(~df_stepA['symbol'].isnull()) | (~df_stepA['gene_child'].isnull())]

df_stepA.columns = ['phenopacket','ORPHAcode','score','rank','gene_P','symbol','ORPHAcode_child','gene_child']


logger.info("A2\tADD genes ")

df_stepA.to_csv(PATH_OUTPUT+"stepA2.tsv", sep='\t',index=False)

logger.info("A2\tEXPORT DF OUTPUT : contain all stepA2 interactions\t{}".format(len(df_stepA)))




logger.info("A2\tEND ")


########################################################## 
# A1  A1  A1  A1  A1  A1  A1  A1  A1  A1  A1  A1  A1  A1
##########################################################

logger.info('A1')

# we rebuild the dataframe base on pd6 and pd_child
# we remove line where orphacode parent AND child have genes
# we keep line where orphaparent have gene AND NOT child OR orphacode chidl have gene AND NOT parents 
load=""
all_interractions=set()

dict_stepA_df = df_stepA.to_dict('index')
for value in dict_stepA_df.values():
    onepheno = value['phenopacket']
    oneORPHA = value['ORPHAcode']
    oneORPHAgene= str(value['symbol'])

    onescore = value['score']
    onePgene = str(value['gene_P'])
    onerank = str(value['rank'])



    oneORPHAchild = value['ORPHAcode_child']
    oneORPHAchildgene = value['gene_child']

    if onePgene == oneORPHAgene:
        load = 'yes'
    elif onePgene == oneORPHAchildgene:
        load = "yes"
    
    if load == 'yes':
        
        all_interractions.add((
            onepheno,
            onescore ,
            onePgene,
            oneORPHA,
            oneORPHAgene,
            oneORPHAchild ,
            oneORPHAchildgene,
            onerank
            ))
    load = ""

stepA1_all = pd.DataFrame(all_interractions, columns=["phenopacket",'score','gene_P','ORPHAcode','symbol','ORPHAcode_child','gene_child','rank'])


stepA1_all = stepA1_all.dropna(subset=['symbol','gene_P'])
stepA1_all.to_csv(PATH_OUTPUT+r"stepA1.tsv", sep='\t',index=False)
# stepA1_all['phenopacket'].to_csv(PATH_OUTPUT+"nb_case_A1.tsv", sep='\t',index=False)

logger.info("A1\tEXPORT DF OUTPUT : contain all stepA1 interactions\t{}".format(len(stepA1_all)))

logger.info("A1\tEND ")

logger.info("END\tA1_A2\n ")

