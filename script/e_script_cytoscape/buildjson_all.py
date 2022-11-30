import pandas as pd
import json
import os
import datetime
from time import perf_counter
from multiprocessing import Process

from SolveRD.script.path_variable import *

import logging
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()


def read_osfiles(path):
    """  store files names in a list"""
    os_files = os.listdir(path)
    len(os_files)

    return list(os_files)


##############################################################################################################################



##############################################################################################################################

def WP_type_interaction(thenode, df_of_interest):
    """ define node type interaction resnik (come from solved-rd algo) wikipath (come from the related DB) or both (possible with node high degree)"""
    ###### wikipathway ######
    WK_type_var = ""
    # TEST type interaction whay kind of interaction the node can have
    WP_interaction = df_of_interest[(df_of_interest['node1'] == thenode) | (df_of_interest['node2'] == thenode)][
        'type_interaction']
    WP_interaction_list = list(set(WP_interaction))
    if (('wikipathways' in WP_interaction_list) and ('resnik' in WP_interaction_list)):
        WK_type_var = 'both'
    elif (('resnik' in WP_interaction_list) and ('wikipathways' not in WP_interaction_list)):
        WK_type_var = 'resnik'
    elif (('wikipathways' in WP_interaction_list) and ('resnik' not in WP_interaction_list)):
        WK_type_var = 'wikipathways'
    else:
        # print("/tFUNCTION WP_type_interaction \t PROBLEME ?")
        pass
    # print("/tFUNCTION WP_type_interaction \t",thenode," :\t",WK_type_var)
    return WK_type_var


##############################################################################################################################
def build_datajsondict(listofdict):
    """  build a json format available for cytoscape js """
    # input a list of several dict
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


def set_autocomplete(path, case, T_step, list_nodes):
    # for each nodes (all cases/genes/orphacode) we store it into a txt file for the autocomplete
    with open(path + case + "_" + T_step + ".txt", 'w') as f:
        for onenode in set(list_nodes):
            # we use ; to split because it's esier to handle in the js files
            f.write(str(onenode) + ";")
        # print(type_step,"\tEND  AUTOCOMPLETE\t",removetsv)


def get_variant(set_step, df_variantimported, phenointerest):
    """  function which detect if a case have variant AND if it's pathogenetic and what type of pathogenicity"""
    listofvar_accepted = [
        'Pathogenic',
        'Likely pathogenic',
        'Pathogenic/Likely pathogenic',
        'Pathogenic/Likely pathogenic, risk factor',
        'Pathogenic//Likely pathogenic',
        'Benign/Likely benign//Pathogenic//Pathogenic',
        'Benign//Pathogenic',
        'Conflicting interpretations of pathogenicity',
        'Conflicting interpretations of pathogenicity, other',
        'Conflicting interpretations of pathogenicity//Uncertain significance',
        'Conflicting interpretations of pathogenicity, association, risk factor'
    ]
    info_var = set()
    # variant_phenopacket = list(df_variantimported.iloc[:, :1])
    variant_phenopacket = df_variantimported['phenotips_id']

    # condition to see if the case of interest is into the list set_step extract throught the cytoscape file (
    # it mean that the case of interest is available in the network)
    if phenointerest in set_step:
        # print("\tFUNCTION get_variant \t ",phenointerest,' is in df step')
        if str(phenointerest) in set(variant_phenopacket):
            # print("\tFUNCTION get_variant \t ",phenointerest," have variant")
            # tinyDF containt only the dataframe of one case
            tinyDF = df_variantimported[df_variantimported['phenotips_id'] == phenointerest]
            # for each case we add the type of variant and its gene
            dict_tinyDF = tinyDF.to_dict('index')
            for value in dict_tinyDF.values():
                onegenePHENO = value['gene']
                onevariant = value['ClinVar_ClinicalSignificance']

                if onevariant in listofvar_accepted:
                    if onevariant == "Pathogenic":
                        info_var.add((onevariant, onegenePHENO))
                    if onevariant == "Likely pathogenic":
                        info_var.add((onevariant, onegenePHENO))
                    if onevariant == "Pathogenic/Likely pathogenic":
                        info_var.add((onevariant, onegenePHENO))

                    if onevariant == "Pathogenic/Likely pathogenic, risk factor":
                        info_var.add(("Pathogenic/Likely pathogenic", onegenePHENO))
                    if onevariant == "Pathogenic//Likely pathogenic":
                        info_var.add(("Pathogenic/Likely pathogenic", onegenePHENO))
                    if onevariant == "Benign/Likely benign//Pathogenic//Pathogenic":
                        info_var.add(("Benign/Likely benign/Pathogenic/Pathogenic", onegenePHENO))
                    if onevariant == "Benign//Pathogenic":
                        info_var.add((onevariant, onegenePHENO))

                    if onevariant == "Conflicting interpretations of pathogenicity":
                        info_var.add((onevariant, onegenePHENO))
                    if onevariant == "Conflicting interpretations of pathogenicity, other":
                        info_var.add(("Conflicting interpretations of pathogenicity", onegenePHENO))
                    if onevariant == "Conflicting interpretations of pathogenicity//Uncertain significance":
                        info_var.add(("Conflicting interpretations of pathogenicity", onegenePHENO))
                    if onevariant == "Conflicting interpretations of pathogenicity, association, risk factor":
                        info_var.add(("Conflicting interpretations of pathogenicity", onegenePHENO))
                        # gene_variant =set( df_variantimported[ df_variantimported['phenotips_id'] == phenointerest]['gene'].to_list())
                    # type_variant =set( df_variantimported[ df_variantimported['phenotips_id'] == phenointerest]['ClinVar_ClinicalSignificance'].to_list())
        else:
            # print("\tFUNCTION get_variant \t ",phenointerest," don't have variant")
            info_var = []
    else:
        # print("\tFUNCTION get_variant \t ",phenointerest,' is not in df step')
        info_var = []
    # return a set of tuple for the case of interest genes variant /type of variant
    return info_var



def set_variant( list_node, the_df_step_variant, case):

    # detect if there is variant into the case of interest
    info_variant = get_variant(set(list_node), the_df_step_variant, case)
    # output of info_variant is a list of tuple which will help building a dataframe
    df_I_V = pd.DataFrame(info_variant, columns=["type variant", 'gene'])
    # print(type_step,"\tBUILD df variant",df_info_variant.values.tolist())
    # print(type_step,"\tEXPORT VARIANT in tsv format")

    # df_I_V.to_csv(path + T_step + '_variant_' + str(case) + '.tsv', sep='\t')
    # print(type_step,"\tEND  VARIANT\t",case)
    return df_I_V


##############################################################################################################################################################################
##############################################################################################################################################################################

def build_variant_AC_json(type_step, df_step_variant, nodetype_pheno, nodetype_ERN, nodetype_orpha, df_gene,
                          cytoscapeWiki_folder):
    """ the main function : from cytoscape df and wikipath cytoscape df :
            -> get variant build txt output variant
            -> autocomplete txt output
            -> create json format for cytoscape js two dict nodes and egdes
            -> convert dict into json """
    # store all names files of one folder

    ###############################################################################################################################################################################
    gene_all_orphanet = df_gene['symbol'].tolist()

    i = 0
    # cytoscapeWiki_folder = ['cytoscape_P0000044_A2.tsv']
    for cytoscapefile in cytoscapeWiki_folder:
        if type_step in cytoscapefile:
            i = i + 1
            # cytoscape file is like cytoscape_P0000044_A2.tsv once slited : ['cytoscape', 'P0000044', 'A2.tsv']
            splited = cytoscapefile.split('_')
            removetsv = str(splited[1])
            # removetsv = "P0000044"

            # print(type_step,"\tTHE CASE : ", removetsv)

            #####################
            # df_cyto_or_WK_step = pd.read_csv(PATH_OUTPUT_CYTOSCAPE+str(cytoscapefile),sep='\t',usecols=[1,2,3])
            # read the cytoscape tsv file with WIKIPATHWAY !!!
            df_cyto_or_WK_step = pd.read_csv(
                PATH_OUTPUT_CYTOSCAPE +"cytoscape_" + str(removetsv)+'_'+str(type_step)+'.tsv',
                sep='\t', usecols=[1, 2, 3, 4])
            #####################

            #####################
            # store all nodes (case,gene,orphacode) in one liste
            all_nodes = list(df_cyto_or_WK_step['node1']) + list(df_cyto_or_WK_step['node2'])
            #####################

            #####################
            # print(type_step,"\tSTART variant : ", removetsv)
            # export vairant df AND return a list of all variant for the case (removetsv) and its df
            df_info_variant = set_variant( all_nodes, df_step_variant, removetsv)
            # a list of all variant for the case (removetsv)
            all_variant = list(df_info_variant['gene'])
            #print(all_variant)
            # print(type_step,"\tall variant for the COI : ",str(removetsv),"\t",len(all_variant),all_variant)

            #####################

            #####################
            # print(type_step,"\tSTART AUTOCOMPLETE : ", removetsv)
            set_autocomplete(PATH_OUTPUT_AUTOCOMPLETE, removetsv, type_step, all_nodes)
            #####################

            #####################
            # load TSV file for the orphachild part
            df_step_nocyto = pd.read_csv(PATH_OUTPUT_TSV + str(removetsv) + "_" + type_step + ".tsv", sep='\t')
            df_for_orphaCHILD = pd.DataFrame(df_step_nocyto, columns=["ORPHAcode", 'ORPHAcode_child'])

            if type_step == "C2":
                list_orpha_C2_onenode = list(df_step_nocyto["ORPHAcode"]) + list(df_step_nocyto["ORPHAcode_C2"])
                list_child_C2_onenode = list(df_step_nocyto["ORPHAcode_child"]) + list(
                    df_step_nocyto["ORPHAcode_child_C2"])

                # create a tuple for the df
                list_for_df = list(zip(list_orpha_C2_onenode, list_child_C2_onenode))
                # here we create a new df which have the same name col as the other step
                df_for_orphaCHILD = pd.DataFrame(list_for_df, columns=["ORPHAcode", 'ORPHAcode_child'])

            df_for_orphaCHILD = df_for_orphaCHILD.loc[df_for_orphaCHILD["ORPHAcode_child"] != "ORPHA:nan"]
            df_for_orphaCHILD = df_for_orphaCHILD.dropna()
            #####################

            # print(type_step,"\tSTART BUILDING JSON : ", removetsv)
            # print(type_step,"\tNODE")

            nodedict = {}
            nodelist = []

            for onenode in set(all_nodes):
                ###### wikipathway ######
                WK_var = WP_type_interaction(onenode, df_cyto_or_WK_step)

                # for each node thus each orphacodes or cases or genes of the case of interest

                # condition to detect if the node is a case
                if onenode in set(nodetype_pheno['case']):
                    # print(type_step,"\tINSIDE THE CASE type : ", onenode)

                    nodeinfo = nodetype_pheno[nodetype_pheno['case'] == onenode].values[
                        0]  # it's a list ['P0002635', 'UNSOLVED']
                    if onenode in set(nodetype_ERN['element']):
                        nodeERN = nodetype_ERN[nodetype_ERN['element'] == onenode].values[
                            0]  # it's a list ['P0002635', 'NMD']
                        ERNinfo = nodeERN[1]
                    else:
                        ERNinfo = "no"

                    nodedict = {
                        "id": str(nodeinfo[0]),  # P0002635
                        "group": nodeinfo[1],  # UNSOLVED
                        "type_variant": 'no',
                        # need to put this key type_variant because it for gene but not use for case
                        "ERN": ERNinfo,  # NMD
                        "id_gene": 'no',
                        "couleur_orpha": 'no',
                        'type_orpha': 'no',
                        'type_interaction': WK_var
                    }

                    # if it's not a case condition to detect if the node is an orphacode
                elif onenode in set(nodetype_orpha['ORDO']):
                    # print(type_step,"\tINSIDE THE ORPHA:script : ", onenode)

                    orphainfo = list(
                        nodetype_orpha[nodetype_orpha['ORDO'] == onenode].values[
                            0])  # it's a list ['ORPHA:98', 'Disease']

                    # if it s a child
                    if (onenode in list(df_for_orphaCHILD['ORPHAcode_child'])):
                        orphachild = onenode  # store the name node because it s a child

                    elif (onenode in list(df_for_orphaCHILD['ORPHAcode'])) and (
                            onenode in list(df_for_orphaCHILD['ORPHAcode_child'])):
                        orphachild = onenode  # store the name it s a child even is he s available on parent col

                    elif onenode in list(df_for_orphaCHILD['ORPHAcode']):
                        orphachild = "no"  # it s not a child
                    else:
                        orphachild = "no"  # it s not a child

                    nodedict = {
                        "id": orphainfo[0],  # ORPHA:353327
                        "group": "ORPHA:code",  # this can be change if we want de specify disease type
                        "type_variant": 'no',
                        # need to put this key type_variant because it for gene but not use for case or orphacode
                        "ERN": "no",  # need to put this key ERN because it for case but not use for orphacode
                        "id_gene": 'no',
                        "couleur_orpha": orphachild,
                        'type_orpha': orphainfo[1],
                        'type_interaction': WK_var

                    }


                # if it's not a case and not an orphacode thus it's a gene
                elif onenode in set(gene_all_orphanet):

                    # print(type_step,"\tINSIDE THE GENE : ", onenode)

                    # all_gene = df_gene[(df_gene['symbol'] == onenode) | (df_gene["ref"] == onenode)]['id'].values[0]
                    all_gene = df_gene[(df_gene['symbol'] == onenode)].values[0]

                    # gene_number = str(all_gene)

                    # if the gene have variants
                    if onenode in all_variant:
                        nodedict = {
                            "id": onenode,
                            "group": "variant",
                            "type_variant": str(
                                df_info_variant[df_info_variant['gene'] == onenode]['type variant'].values[0]),
                            "ERN": "no",
                            "id_gene": str(all_gene[1]),
                            "couleur_orpha": 'no',
                            'type_orpha': 'no',
                            'type_interaction': WK_var

                        }
                    # else the gene don't have variant
                    else:
                        nodedict = {
                            "id": onenode,
                            "group": "gene",
                            "type_variant": 'no',
                            "ERN": "no",
                            "id_gene": str(all_gene[1]),
                            "couleur_orpha": 'no',
                            'type_orpha': 'no',
                            'type_interaction': WK_var
                        }

                elif WK_var == "wikipathways":
                    nodedict = {
                        "id": onenode,
                        "group": "wikipathways",
                        "type_variant": 'no',
                        "ERN": "no",
                        "id_gene": 'no',
                        "couleur_orpha": 'no',
                        'type_orpha': 'no',
                        'type_interaction': WK_var
                    }

                else:
                    nodedict = {
                        "id": onenode,
                        "group": "other",
                        "type_variant": 'no',
                        "ERN": "no",
                        "id_gene": 'no',
                        "couleur_orpha": 'no',
                        'type_orpha': 'no',
                        'type_interaction': 'no'
                    }

                # add a dict into a list (which is a list of dict)
                nodelist.append(nodedict)
                # reset the dict  for the next node
                nodedict = {}

            # same process for edges
            edgedict = {}
            edgelist = []

            # print(type_step,"\tEDGES")
            dict_df_cyto_or_WK_step = df_cyto_or_WK_step.to_dict('index')
            # df is the cytoscape tsv
            for value in dict_df_cyto_or_WK_step.values():
                # for each row we extract node1 score node2 because it's for edges
                node1 = value['node1']
                node2 = value['node2']
                score = str(value['score'])

                edgedict = {
                    "id": node1 + node2,  # need an unique id thus merge both node to make sure id is unique
                    "source": node1,
                    "target": node2,
                    "score": score
                }
                # add a dict into a list (which is a list of dict)
                edgelist.append(edgedict)
                # reset the dict  for the next row
                edgedict = {}

            # print(type_step,"\tedges+nodes in a dict")
            # build a new list of dict to match the json cytoscape js format
            global_node = build_datajsondict(nodelist)
            global_edge = build_datajsondict(edgelist)

            # create a new dict which will containt the list of dict
            globaldict = {}

            globaldict["nodes"] = global_node
            globaldict["edges"] = global_edge

            # globaldict is a dict of list which have 2 keys nodes and egdes
            # globaldict respect the json format for cytoscape js
            # print(type_step,"\tEXPORT JSON",removetsv)
            with open(PATH_OUTPUT_JSON + type_step + "_" + removetsv + ".json",
                      'w') as f:
                f.write(json.dumps(globaldict, indent=2))

            logger.info('buildjsonall.py\t{}\tTIME END  JSON: {}\t {}FOR THE CASE'.format(type_step,datetime.datetime.utcnow(),removetsv))
    logger.info('buildjsonall.py\t{}\tINPUT cytsocape files size : \t{}'.format(type_step, i))

##############################################################################################################################################################################


##############################################################################################################################################################################
##############################################################################################################################################################################
if __name__ == '__main__':

    # time start
    tstart = perf_counter()

    tab_step = ["A1", "A2", 'B1', 'B2', 'C1', 'C2']
    # tab_step = ["A2"]

    variant_step = ["StepA1_variants_28.1.22.xlsx", "StepA2_variants_28.1.22.xlsx", 'StepB1_variants_28.1.22.xlsx',
                    'StepB2_variants_25.2.22.xlsx', 'StepC1_variants_25.2.22.xlsx', 'StepC2_variants_25.2.22.xlsx']
    # tab_step = ["A1"]
    # variant_step = ["StepA2_variants_28.1.22.xlsx"]

    # type phenopacket UNSOLVED or SOLVED
    NT_pheno = pd.read_csv(PATH_OUPUT_NODE_TYPE_PHENO, sep='\t')
    # dataframe of orphacode related to type of subtype or type of disorder
    NT_orpha = pd.read_csv(PATH_OUPUT_NODE_TYPE_ORPHA, sep='\t')
    # dataframe of cases related to ERNs
    NT_ERN = pd.read_csv(PATH_OUPUT_NODE_TYPE_ERN, sep='\t')
    ## GENE ID all gene from orphanet (for the gene links in the app)
    NT_gene = pd.read_csv(PATH_OUPUT_NODE_TYPE_GENE, sep='\t')

    # wikipathway folder dont have a cytoscape_homepage_ALL.tsv file but keep this code to prevent
    cyto_F = read_osfiles(PATH_OUTPUT_CYTOSCAPE)

    try:
        # Create target Directory
        cyto_F.remove("cytoscape_homepage_ALL.tsv")
        logger.info('buildjsonall.py\tRemove cytoscape_homepage_ALL.tsv')

    except:
        logger.info('buildjsonall.py\tcytoscape_homepage_ALL.tsv Already removed')



    # test_df = build_variant_AC_json(tab_step[0], pd.read_excel(PATH_INPUT_VARIANT + variant_step[1], usecols=[2, 8, 37]), NT_pheno, NT_ERN,NT_orpha, NT_gene, cyto_F)



    list_proc = []
    for i in range(len(tab_step)):
        logger.info('buildjsonall.py\t{}\t{}\t'.format(tab_step[i],variant_step[i]))

        P_step = Process(target=build_variant_AC_json, args=(
        tab_step[i], pd.read_excel(PATH_INPUT_VARIANT + variant_step[i], usecols=[2, 8, 37]), NT_pheno, NT_ERN,
        NT_orpha, NT_gene, cyto_F))
        list_proc.append(P_step)


    for onejob in list_proc:
        onejob.start()

    for onejob in list_proc:
        onejob.join()

        logger.info('buildjsonall.py\tIs Process alive ?\t{}\t'.format(onejob.is_alive()))
        logger.info('buildjsonall.py\tTIME END:\tprocess id: {}'.format(datetime.datetime.utcnow(),os.getpid()))


    tstop = perf_counter()

    logger.info('buildjsonall.py\tIME build all df :{}'.format((tstop - tstart)))
    logger.info('buildjsonall.py\tEND')
