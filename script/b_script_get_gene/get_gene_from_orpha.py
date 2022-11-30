from SolveRD.script.path_variable import *

import pandas as pd
import logging

# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')

logger = logging.getLogger()

logger.info("START\tget_gene_from_orpha.py\n ")


list_pd1 = []
orphacode_id = disordergroup = disordertype = ""
with open(PATH_INPUT_PRODUCT_PD1 ,'r',encoding="ISO-8859-1") as f:
    ligne = f.readline()
    while ligne != "":
        if "<JDBOR" in ligne:
            logger.info("get_gene_from_orpha.py\tPD_1:\t{}".format(ligne))
        elif "<DisorderList" in ligne:
            logger.info("get_gene_from_orpha.py\tPD_1:\t{}".format(ligne))
        elif '<OrphaCode>' in ligne:
            new_ligne = ligne.replace('<OrphaCode>' ,'ORPHA:').strip()
            new_ligne = new_ligne.replace('</OrphaCode>' ,'').strip()

            orphacode_id = new_ligne
        elif  '<DisorderType' in ligne:
            dec_ligne =f.readline()
            new_ligne = dec_ligne.replace('<Name lang="en">', '').strip()
            new_ligne = new_ligne.replace('</Name>', '').strip()
            disordertype = new_ligne
        elif '<DisorderGroup' in ligne:
            dec_ligne = f.readline()
            new_ligne = dec_ligne.replace('<Name lang="en">', '').strip()
            new_ligne = new_ligne.replace('</Name>', '').strip()
            disordergroup = new_ligne

        list_pd1.append((orphacode_id,disordertype, disordergroup))
        ligne = f.readline()
logger.info('get_gene_from_orpha.py\tPD1:\tread xml product 1 HPO which containt all orphas with HPO ')

list_pd6 = []
orphacode_id = disordergroup = disordertype = ""
symbol = 'symbol'
with open(PATH_INPUT_PRODUCT_PD6, 'r',encoding="ISO-8859-1") as f:
    ligne = f.readline()

    while ligne != "":
        if "<JDBOR" in ligne:
            logger.info("get_gene_from_orpha.py\tPD_6:\t{}".format(ligne))
        elif "<DisorderList" in ligne:
            logger.info("get_gene_from_orpha.py\tPD_6:\t{}".format(ligne))
        elif '<OrphaCode>' in ligne:
            new_ligne = ligne.replace('<OrphaCode>', 'ORPHA:').strip()
            new_ligne = new_ligne.replace('</OrphaCode>', '').strip()
            orphacode_id = new_ligne
        elif '<DisorderType' in ligne:
            dec_ligne = f.readline()
            new_ligne = dec_ligne.replace('<Name lang="en">', '').strip()
            new_ligne = new_ligne.replace('</Name>', '').strip()
            disordertype = new_ligne
        elif '<DisorderGroup' in ligne:
            dec_ligne = f.readline()
            new_ligne = dec_ligne.replace('<Name lang="en">', '').strip()
            new_ligne = new_ligne.replace('</Name>', '').strip()
            disordergroup = new_ligne
        elif '<Symbol>' in ligne:

            new_ligne = ligne.replace('<Symbol>', '').strip()
            new_ligne = new_ligne.replace('</Symbol>', '').strip()
            symbol = new_ligne

        list_pd6.append((orphacode_id,disordertype, disordergroup,symbol))
        symbol = 'symbol'
        ligne = f.readline()

logger.info('get_gene_from_orpha.py\tPD6:\tread xml product 6 genes which containt all orphas with genes ')



get_pd_child =list()

orphacode_id_parent = orphacode_id_child = ""
symbol = 'symbol'
only_for_first_el = 0
with open(PATH_INPUT_PRODUCT_PD_ORPHACHILD, 'r',encoding="ISO-8859-1") as f:
    ligne = f.readline()
    while ligne != "":
        if "<JDBOR" in ligne:
            logger.info("get_gene_from_orpha.py\tPD_CHILD:\t{}".format(ligne))
        elif "<DisorderList" in ligne:
            logger.info("get_gene_from_orpha.py\tPD_CHILD:\t{}".format(ligne))
        elif ('<OrphaCode>' in ligne) :
            dec_ligne = f.readline()
            dec_ligne = f.readline()
            if '<ClassificationNodeList' in dec_ligne:
                new_ligne = ligne.replace('<OrphaCode>', 'ORPHA:').strip()
                new_ligne = new_ligne.replace('</OrphaCode>', '').strip()
                orphacode_id_parent=new_ligne
            else:
                new_ligne = ligne.replace('<OrphaCode>', 'ORPHA:').strip()
                new_ligne = new_ligne.replace('</OrphaCode>', '').strip()
                orphacode_id_child = new_ligne
        elif '<Symbol>' in ligne:

            new_ligne = ligne.replace('<Symbol>', '').strip()
            new_ligne = new_ligne.replace('</Symbol>', '').strip()
            symbol = new_ligne

        get_pd_child.append((orphacode_id_parent,orphacode_id_child, symbol))
        symbol = 'symbol'
        ligne = f.readline()

logger.info('get_gene_from_orpha.py\tPD_CHILD:\tread xml special product child which containt all orphas with child ')

# build and filter df for child
df_pd_child = pd.DataFrame(get_pd_child, columns=["ORPHACode",'ORPHACode_child',"gene_child"])
df_pd_child = df_pd_child.drop_duplicates()
df_pd_child_filtered = df_pd_child[df_pd_child['gene_child'] != 'symbol']


#####################
# build and filter df for orpha genes
df_pd6 = pd.DataFrame(list_pd6, columns=['ORPHACode', 'DisorderType_6', 'DisorderGroup_6', 'symbol'])
df_pd6=df_pd6.drop_duplicates()
df_pd6_filtered = df_pd6[df_pd6['DisorderType_6'] != 'Category']

#####################

# build and filter df for orpha hpo
df_pd1 = pd.DataFrame(list_pd1, columns=['ORPHACode', 'DisorderType_1', 'DisorderGroup_1'])
df_pd1=df_pd1.drop_duplicates()
df_pd1 = df_pd1.drop_duplicates(subset='ORPHACode', keep="last")
df_pd1_filtered = df_pd1[df_pd1['DisorderGroup_1'] != 'Group of disorders']
#####################

# MERGE genes and hpo (pd 6 and 1)
df_pd1_pd6 = pd.merge(df_pd6_filtered, df_pd1_filtered, on='ORPHACode',how='outer')
df_pd1_pd6 = df_pd1_pd6.drop('DisorderType_6', 1)
df_pd1_pd6 = df_pd1_pd6.drop('DisorderGroup_6', 1)

#####################
# MERGE ALL genes hpo and child
df_pd1_pd6_child = pd.merge(df_pd1_pd6, df_pd_child_filtered, on='ORPHACode',how='outer')


# loops to keep only parent which have gene but child not OR parent with no gene but child yes
dict_df_pd1_pd6_child = df_pd1_pd6_child.to_dict('index')
all_interractions = list()
for value in dict_df_pd1_pd6_child.values():
    oneORPHAchildgene = value['gene_child']
    oneORPHA = value['ORPHACode']
    oneORPHAgene = value['symbol']
    oneORPHAchild = value['ORPHACode_child']

    if pd.isna(oneORPHAchildgene):
        if pd.isna(oneORPHAgene):
            pass
        else:
            load = "yes"  # gene child empty but not its gene parent

    if pd.isna(oneORPHAgene):
        if pd.isna(oneORPHAchildgene):
            #both empty
            pass
        else:
            # orpha child filled
            load = "yes"  # gene parent not empty but chidl yes

    if (pd.isna(oneORPHAgene) and pd.isna(oneORPHAchildgene)):
        # both empty
        pass
    else:
        load = "yes"

    if load == "yes":
        all_interractions.append((
            oneORPHA,
            oneORPHAgene,
            oneORPHAchild,
            oneORPHAchildgene,
        ))

    load = ""


# build and filter from the merged one now filtered for child/parent genes
df_product = pd.DataFrame(all_interractions, columns=["ORPHAcode", 'symbol', 'ORPHAcode_child', 'gene_child'])
df_product = df_product.drop_duplicates()
df_product = df_product[df_product["symbol"] != "symbol"]

df_pd6=df_pd6.drop_duplicates()




logger.info('get_gene_from_orpha.py\texport the df PATH_INPUT_PRODUCT_1_6_child')
df_product.to_csv(PATH_OUTPUT_PRODUCT_1_6_child, sep='\t',index=False)


logger.info("END\tget_gene_from_orpha.py\n ")

