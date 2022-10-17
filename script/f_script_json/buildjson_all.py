import pandas as pd
import json 
import os

import datetime
from time import perf_counter

from multiprocessing import Queue, Process, Pool

from script.path_variable import *
 
def read_osfiles(path):
    """  store files names in a list"""
    os_files = os.listdir(path)
    len(os_files)

    return list(os_files)


##############################################################################################################################


def get_variant(set_step, df_variantimported,phenointerest):
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
    variant_phenopacket = list( df_variantimported.iloc[:,:1])
    # condition to see if the case of interest is into the list set_step extract throught the cytoscape file (
    # it mean that the case of interest is available in the network)
    if phenointerest in set_step:
        # print("/tFUNCTION get_variant \t ",phenointerest,' is in df step')
        if str(phenointerest) in set(variant_phenopacket):
            # print("\tFUNCTION get_variant \t ",phenointerest," have variant")
            # tinyDF containt only the dataframe of one case
            tinyDF =  df_variantimported[ df_variantimported['phenotips_id'] == phenointerest]
            # for each case we add the type of variant and its gene
            dict_tinyDF = tinyDF.to_dict('index')
            for value in dict_tinyDF.values():
                onegenePHENO = value['gene']
                onevariant = value['ClinVar_ClinicalSignificance']

                if onevariant in listofvar_accepted:
                    if onevariant == "Pathogenic":
                        info_var.add((onevariant,onegenePHENO)) 
                    if onevariant == "Likely pathogenic":
                        info_var.add((onevariant,onegenePHENO)) 
                    if onevariant == "Pathogenic/Likely pathogenic":
                        info_var.add((onevariant,onegenePHENO)) 

                        
                    if onevariant == "Pathogenic/Likely pathogenic, risk factor":
                        info_var.add(("Pathogenic/Likely pathogenic",onegenePHENO)) 
                    if onevariant == "Pathogenic//Likely pathogenic":
                        info_var.add(("Pathogenic/Likely pathogenic",onegenePHENO))
                    if onevariant == "Benign/Likely benign//Pathogenic//Pathogenic":
                        info_var.add(("Benign/Likely benign/Pathogenic/Pathogenic",onegenePHENO)) 
                    if onevariant == "Benign//Pathogenic":
                        info_var.add((onevariant,onegenePHENO)) 

                        
                    if onevariant == "Conflicting interpretations of pathogenicity":
                        info_var.add((onevariant,onegenePHENO)) 
                    if onevariant == "Conflicting interpretations of pathogenicity, other":
                        info_var.add(("Conflicting interpretations of pathogenicity",onegenePHENO)) 
                    if onevariant == "Conflicting interpretations of pathogenicity//Uncertain significance":
                        info_var.add(("Conflicting interpretations of pathogenicity",onegenePHENO)) 
                    if onevariant == "Conflicting interpretations of pathogenicity, association, risk factor":
                        info_var.add(("Conflicting interpretations of pathogenicity",onegenePHENO))                                                                         
                    # gene_variant =set( df_variantimported[ df_variantimported['phenotips_id'] == phenointerest]['gene'].to_list())
                    # type_variant =set( df_variantimported[ df_variantimported['phenotips_id'] == phenointerest]['ClinVar_ClinicalSignificance'].to_list())
        else :
            # print("\tFUNCTION get_variant \t ",phenointerest," don't have variant")
            info_var = []      
    else :
        # print("\tFUNCTION get_variant \t ",phenointerest,' is not in df step')
        info_var = []
    # return a set of tuple for the case of interest genes variant /type of variant
    return info_var

##############################################################################################################################

def WP_type_interaction(thenode,df_of_interest):
    """ define node type interaction resnik (come from solved-rd algo) wikipath (come from the related DB) or both (possible with node high degree)"""
    ###### wikipathway ######
    WK_type_var = ""
    # TEST type interaction whay kind of interaction the node can have
    WP_interaction = df_of_interest[(df_of_interest['node1'] == thenode )| (df_of_interest['node2'] == thenode)]['type_interaction']
    WP_interaction_list= list(set(WP_interaction))
    if (('wikipathway' in WP_interaction_list) and  ('resnik' in WP_interaction_list)):
        WK_type_var = 'both'
    elif (('resnik' in WP_interaction_list) and  ('wikipathway' not in WP_interaction_list)):
        WK_type_var = 'resnik'
    elif (('wikipathway' in WP_interaction_list )and  ('resnik' not in WP_interaction_list)):
        WK_type_var = 'wikipathway'
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
            "data" : onedict
        }
        # add this new dict into a new list
        listlistofdict.append(datadict)
    return listlistofdict




##############################################################################################################################################################################
##############################################################################################################################################################################

def build_variant_AC_json(type_step, df_step_variant):
    """ the main function : from cytoscape df and wikipath cytoscape df :
            -> get variant build txt output variant
            -> autocomplete txt output
            -> create json format for cytoscape js two dict nodes and egdes
            -> convert dict into json """
    # store all names files of one folder

    ###############################################################################################################################################################################
    ## JSON
    # nodetype tsv

    # type phenopacket UNSOLVED or SOLVED
    nodetype_pheno = pd.read_csv(PATH_OUPUT_NODE_TYPE_PHENO,sep='\t')
    # dataframe of orphacode related to type of subtype or type of disorder
    nodetype_orpha = pd.read_csv(PATH_OUPUT_NODE_TYPE_ORPHA, sep='\t')
    # dataframe of cases related to ERNs
    nodetype_ERN = pd.read_csv(PATH_OUPUT_NODE_TYPE_ERN,sep='\t')

    ## GENE ID all gene from orphanet (for the gene links in the app)
    df_gene = pd.read_csv(PATH_OUPUT_NODE_TYPE_GENE,sep='\t')
    gene_all_orphanet = list(df_gene['symbol'])

    ###############################################################################################################################################################################

    cytoscape_folder = read_osfiles(PATH_OUTPUT_CYTOSCAPE_WIKIPATHWAY)
    print(type_step, "\tINPUT cytsocape files size : \t", len(cytoscape_folder))

    no_home_file_STEP_cytoscape_folder = []

    for cyto_file in cytoscape_folder:
        if type_step in cyto_file:
            no_home_file_STEP_cytoscape_folder.append(cyto_file)


    for cytoscapefile in no_home_file_STEP_cytoscape_folder:
        splited = cytoscapefile.split('_')
        # cytoscape file is like wikipathway_cytoscape_P0000044_A2.tsv once slited : ['wikipathway','cytoscape', 'P0000044', 'A2.tsv']
        removetsv = str(splited[2])
        # removetsv = "P0000044"

        # print(type_step,"\tTHE CASE : ", removetsv)

        # read the cytoscape tsv file
        # df_cyto_or_WK_step = pd.read_csv(PATH_OUTPUT_CYTOSCAPE+str(cytoscapefile),sep='\t',usecols=[1,2,3])
        # read the cytoscape tsv file with WIKIPATHWAY !!!
        df_cyto_or_WK_step = pd.read_csv(
            PATH_OUTPUT_CYTOSCAPE_WIKIPATHWAY + r"cytoscape_wikipathway_" + str(removetsv) + "_" + str(type_step) + ".tsv",
            sep='\t', usecols=[0, 1, 2, 3])

        # store all nodes (case,gene,orphacode) in one liste
        all_nodes = list(df_cyto_or_WK_step['node1']) + list(df_cyto_or_WK_step['node2'])
        # print(all_nodes)

        # print(type_step,"\tSTART variant : ", removetsv)
        # detect if there is variant into the case of interest
        info_variant = get_variant(set(all_nodes), df_step_variant, removetsv)
        # output of info_variant is a list of tuple which will help building a dataframe
        df_info_variant = pd.DataFrame(info_variant, columns=["type variant", 'gene'])
        # print(type_step,"\tBUILD df variant",df_info_variant.values.tolist())
        # print(type_step,"\tEXPORT VARIANT in tsv format")
        df_info_variant.to_csv(PATH_OUTPUT_VARIANT + type_step + '_variant_' + str(
                removetsv) + '.tsv', sep='\t')
        all_variant = list(df_info_variant['gene'])
        # print(type_step,"\tall variant for the COI : ",str(removetsv),"\t",len(all_variant),all_variant)
        # print(type_step,"\tEND  VARIANT\t",removetsv)

        # print(type_step,"\tSTART AUTOCOMPLETE : ", removetsv)
        # for each nodes (all cases/genes/orphacode) we store it into a txt file for the autocomplete
        with open(PATH_OUTPUT_AUTOCOMPLETE + removetsv + "_" + type_step + ".txt",'w') as f:
            for onenode in set(all_nodes):
                # we use ; to split because it's esier to handle in the js files
                f.write(str(onenode) + ";")
                # print(type_step,"\tEND  AUTOCOMPLETE\t",removetsv)

        df_step_nocyto = pd.read_csv(PATH_OUTPUT_TSV +str(removetsv) + "_" + type_step + ".tsv", sep='\t')
        df_for_orphaCHILD = pd.DataFrame(df_step_nocyto, columns=["ORPHAcode", 'ORPHAcode_child'])

        if type_step == "C2":
            list_orpha_C2_onenode = list(df_step_nocyto["ORPHAcode"]) + list(df_step_nocyto["ORPHAcode_C2"])
            list_child_C2_onenode = list(df_step_nocyto["ORPHAcode_child"]) + list(df_step_nocyto["ORPHAcode_child_C2"])

            # create a tuple for the df
            list_for_df = list(zip(list_orpha_C2_onenode, list_child_C2_onenode))
            # here we create a new df which have the same name col as the other step
            df_for_orphaCHILD = pd.DataFrame(list_for_df, columns=["ORPHAcode", 'ORPHAcode_child'])

        df_for_orphaCHILD = df_for_orphaCHILD.loc[df_for_orphaCHILD["ORPHAcode_child"] != "ORPHA:nan"]
        df_for_orphaCHILD = df_for_orphaCHILD.dropna()

        # minidf_parentORchild
        # this list will be usefull to build the json format
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
                    ERNinfo = "no ERN"

                nodedict = {
                    "id": str(nodeinfo[0]),  # P0002635
                    "group": nodeinfo[1],  # UNSOLVED
                    "type_variant": 'no',  # need to put this key type_variant because it for gene but not use for case
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
                    nodetype_orpha[nodetype_orpha['ORDO'] == onenode].values[0])  # it's a list ['ORPHA:98', 'Disease']

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
                    "group": "ORPHA:script",  # a modifier plus tard bails de type et tout
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
                all_gene = df_gene[(df_gene['symbol'] == onenode) ].values[0]

                # gene_number = str(all_gene)

                # if the gene have variants
                if onenode in all_variant:
                    nodedict = {
                        "id": onenode,
                        "group": "variant",
                        "type_variant": str(
                            df_info_variant[df_info_variant['gene'] == onenode]['type variant'].values[0]),
                        "ERN": "no",
                        "id_gene": str(all_gene[3]),
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
                        "id_gene": str(all_gene[3]),
                        "couleur_orpha": 'no',
                        'type_orpha': 'no',
                        'type_interaction': WK_var
                    }

            elif WK_var == "wikipathway":
                nodedict = {
                    "id": onenode,
                    "group": "wikipathway",
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

        # print("sometime variant are found but not put in the json because gene is not in the network")

        print(type_step, "\tTIME END  JSON: ", datetime.datetime.utcnow(), "FOR THE CASE\t", removetsv)


##############################################################################################################################################################################



##############################################################################################################################################################################
##############################################################################################################################################################################
if __name__ == '__main__':

     # time start
    tstart = perf_counter()


    # build_variant_AC_json("A2",pd.read_excel(PATH_INPUT_VARIANT+"StepA2_variants_28.1.22.xlsx", usecols=[2,8,37]))



    tab_step = ["A1","A2",'B1','B2','C1','C2']
    variant_step =  ["StepA1_variants_28.1.22.xlsx","StepA2_variants_28.1.22.xlsx",'StepB1_variants_28.1.22.xlsx','StepB2_variants_25.2.22.xlsx','StepC1_variants_25.2.22.xlsx','StepC2_variants_25.2.22.xlsx']
     # tab_step = ["A1"]
     # variant_step = ["StepA2_variants_28.1.22.xlsx"]

     list_proc = []
    for i  in range(len(tab_step)):
        print(tab_step[i],'\t',variant_step[i])
        P_step = Process(target=build_variant_AC_json, args=(tab_step[i],pd.read_excel(PATH_INPUT_VARIANT+variant_step[i], usecols=[2,8,37])))
        list_proc.append(P_step)

    for onejob in list_proc:
        onejob.start()

    for onejob in list_proc:
        onejob.join()
        print(f'Is Process alive ?: {onejob.is_alive()}')
        print("TIME END: ", datetime.datetime.utcnow(),'\tprocess id: ', os.getpid())

    tstop = perf_counter()
    print("\n\tTIME build all df : ",(tstop - tstart))


    print("###############\nEND buildjson_all.py\n###############\n")
