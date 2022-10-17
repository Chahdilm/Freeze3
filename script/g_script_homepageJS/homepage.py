import pandas as pd
import time
import datetime

import json

from script.path_variable import *

##############################################################################################################################
def build_datajsondict(listofdict):
    """  build a json format available for cytoscape js """
    datadict = {}
    listlistofdict = []
    for onedict in listofdict:
        # loop for each dict of the list
        datadict = {
            "data": onedict
        }
        # add this new dict into a new list
        listlistofdict.append(datadict)
    return listlistofdict


##############################################################################################################################


##############################################################################################################################
start = time.process_time()



df_case_gene = pd.read_csv(PATH_OUTPUT_SOLVED_GENE,sep = '\t')
input_cases = list(set(df_case_gene['phenopacket']))


# loops on the json phenopacket
case_case = list()
i = 0
while i < len(input_cases):
    # open one phenopact result per loop
    with open(PATH_OUTPUT_RSLT_ALGO_PHENO + str(
            input_cases[i]+".json.json")) as file_phenopacket_result:
        one_phenopacket_result = json.load(file_phenopacket_result)

    # select only the result from a specific algo (these is 8 algo)
    for key_method in one_phenopacket_result['methods']:
        # go inside the dico which has the key name methods
        # if key_method['method'] == "Cosine Weighted Similarity":
        if key_method['method'] == "Resnik (symmetric)":

            for key_results in key_method['results']:
                # go inside the dico which has the key name results to get ORPHA:code info for each phenopacket
                # print("name of the ORPHA:code : \t",key_results['itemid'])

                # phenopacket is list_home_network_cases
                # split part from 'P0000037.json.json' to 'P0000037,
                split_jsonfile_1 = input_cases[i].split('.')
                split_jsonfile_1 = split_jsonfile_1[0]

                # case
                # split part from 'PP:P0007304.json' to 'P0007304,
                split_jsonfile_2 = key_results['itemid'].split(':')
                split_jsonfile_2 = split_jsonfile_2[1].split('.')
                # ['P0007304', 'json'] that why we select the first element
                split_jsonfile_2 = split_jsonfile_2[0]

                # in the beggining put only phenopacket id and orpha:code
                # THEN   store score and rank (but score and rank not save in file yet)
                case_case.append((split_jsonfile_1, split_jsonfile_2, float(key_results['score']),key_results['rank']))

    i = i + 1
# convert the set of tuples into a dataframe
df_case_case = pd.DataFrame(case_case, columns=["phenopacket", 'case', 'score','rank'])
print("HOMEPAGE\tTIME : BUILD the df : ", datetime.datetime.utcnow())


list_phenopacket = set(df_case_case['phenopacket'] )


keep_best_interaction = set()
for onephenopacket in list(list_phenopacket):
    df_single_row  = df_case_case[(df_case_case['phenopacket'] == onephenopacket)   ]
    df_onerow = df_single_row.iloc[:1]
    keep_best_interaction.add((df_onerow['phenopacket'].values[0],df_onerow['case'].values[0],df_onerow['score'].values[0]))


df_keep_best_interaction = pd.DataFrame(keep_best_interaction, columns=["node1", 'node2', 'score'])

df_keep_best_interaction.to_csv(PATH_OUTPUT_CYTOSCAPE+r"\cytoscape_homepage_ALL.tsv",sep='\t')

###############################################################################################################################################################################




## JSON VERSION for the cytoscape js app
# nodetype tsv
# type phenopacket UNSOLVED or SOLVED usless here all of them are solved !!
# nodetype_pheno = pd.read_csv(PATH_INPUT_NODE_TYPE_PHENO,sep='\t')

# dataframe of cases related to ERNs
nodetype_ERN = pd.read_csv(PATH_OUPUT_NODE_TYPE_ERN,sep='\t')

all_nodes = list(df_keep_best_interaction['node1']) + list(df_keep_best_interaction['node2'])
print('START BUILDING JSON')

nodedict = {}
nodelist = []

for onenode in set(all_nodes):
    # for each node thus each orphacodes or cases or genes of the case of interest

    # condition to detect if the node is a case

    # useless because all of them are solved
    # array(['P0002635', 'UNSOLVED'], dtype=object)
    # nodeinfo = nodetype_pheno[nodetype_pheno['cases'] == onenode].values[0]  # it's a list
    if onenode in set(nodetype_ERN['element']):
        # array(['P0002635', 'NMD'], dtype=object)
        nodeERN = nodetype_ERN[nodetype_ERN['element'] == onenode].values[0]  # it's a list
        ERNinfo = nodeERN[1]
    else:
        ERNinfo = "no ERN"

    nodedict = {
        "id": onenode,  # example P0002635
        "group": "SOLVED",
        "type_variant": 'no',  # need to put this key type_variant because it for gene but not use for case
        "ERN": ERNinfo  # example NMD
    }
    # if it's not a case condition to detect if the node is an orphacode


    # add a dict into a list (which is a list of dict)
    nodelist.append(nodedict)
    # reset the dict  for the next node
    nodedict = {}

# same process for edges
edgedict = {}
edgelist = []
dict_df_keep_best_interaction = df_keep_best_interaction.to_dict('index')

# df is the cytoscape tsv
for value in dict_df_keep_best_interaction.values():
    # for each row we extract node1 score node2 because it's for edges
    node1 = value['node1']
    node2 = value['node2']
    score = str(value['score'])

    edgedict = {
        "id": node1 + node2,
        "source": node1,
        "target": node2,
        "score": score
    }

    edgelist.append(edgedict)
    edgedict = {}

print('edges+nodes in a dict')
# build a new list of dict to match the json cytoscape js format
global_node = build_datajsondict(nodelist)
global_edge = build_datajsondict(edgelist)

# create a new dict which will containt the list of dict
globaldict = {}

# globaldict is a dict of list which have 2 keys nodes and egdes
# globaldict respect the json format for cytoscape js
globaldict["nodes"] = global_node
globaldict["edges"] = global_edge

print("EXPORT JSON")
with open(PATH_OUTPUT_JSON + r"\homepage_ALL.json",'w') as f:
    f.write(json.dumps(globaldict, indent=2))

