PATH= r"C:\Users\mchahdil\Documents\Freeze3\\"
PATH_R = r'C:\Program Files\R\R-4.2.1\bin\Rscript'


############################################################################################################################
PATH_SCRIPT = PATH +  r"script\\"

PATH_INPUT = r"C:\Users\mchahdil\Documents\Freeze3\\input\\"
PATH_OUTPUT = r"C:\Users\mchahdil\Documents\Freeze3\\output_5HPO\\"


PATH_INPUT_PRODUCT_PD_ORPHACHILD = PATH_INPUT +"product_ORPHApackets_childs.xml"
PATH_INPUT_PRODUCT_PD6 = PATH_INPUT+"RunSolveRD\en_product6_sep2022.xml"
PATH_INPUT_PRODUCT_PD1 = PATH_INPUT+"RunSolveRD\en_product1_sep2022.xml"

PATH_INPUT_WK = PATH_INPUT+"linkset_for_py4cy_wikipathway//"

######################################

PATH_OUPUT_5HPO_NOPARENT = PATH_OUTPUT + "curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete\\"
PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION = PATH_OUTPUT + "curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\\"
PATH_OUPUT_5HPO_NOPARENT_DF = PATH_OUTPUT + "curation_input\parent_phenopacket.tsv"

######################################

PATH_OUTPUT_RSLT_ALGO = PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION+r"results\results_noduplicates\\"
PATH_OUTPUT_RSLT_ALGO_ORPHA = PATH_OUTPUT_RSLT_ALGO+ r"resultsORDO\\"
PATH_OUTPUT_RSLT_ALGO_PHENO = PATH_OUTPUT_RSLT_ALGO+ r"resultsPhenopackets\\"

######################################

PATH_INPUT_VARIANT = PATH_INPUT+"variant//"
PATH_INPUT_PRODUCT_PHENOPACKET_BRUT = PATH_INPUT + r"RunSolveRD\Freezes1_2_3_noduplicates_rawdata\\"
# ped_master file for the parents phenopackets
PATH_INPUT_PED = PATH_INPUT + r"RunSolveRD\ped_master\\"

# a changer
PATH_INPUT_PRODUCT_PHENOPACKET = PATH_INPUT + r"RunSolveRD\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_RunSolveRD\\"
# PATH_INPUT_PRODCUT_SOLVED_GENE =  PATH_INPUT + r"RunSolveRD\Solved_with_gene_noparents_with_5phenotypes_aftercuration.xlsx"
############################################################################################################################
############################################################################################################################

PATH_OUTPUT_TSV = PATH_OUTPUT+'tsv//'
PATH_OUTPUT_CYTOSCAPE = PATH_OUTPUT+'cytoscape//'
PATH_OUTPUT_CYTOSCAPE_WIKIPATHWAY = PATH_OUTPUT+'cytoscape_wikipathways//'
PATH_OUTPUT_JSON = PATH_OUTPUT+'json//'
PATH_OUTPUT_VARIANT = PATH_OUTPUT+'variant//'
PATH_OUTPUT_AUTOCOMPLETE = PATH_OUTPUT+'autocomplete//'


# output from get_gene_from_case.py
PATH_OUTPUT_JSON_GENE = PATH_OUTPUT + "gene_info/gene_from_json_case.tsv" # use only for B1 22300 to 22499 associations
PATH_OUTPUT_SOLVED_GENE = PATH_OUTPUT + "gene_info/gene_from_solved_case.tsv"

# output from 4_get_gene_from_orpha.py
PATH_OUTPUT_PRODUCT_1_6_child = PATH_OUTPUT + "gene_info/pd_1_6_child.tsv"


############################################################################################################################
PATH_INPUT_NODE_TYPE_ERN = PATH_INPUT + "2022_10_14_cohort_ERN.csv"

PATH_OUPUT_NODE_TYPE_ERN =PATH_OUTPUT +  r"node_type\nodetype_ERN.tsv" # dont have ERN file update its thejamboree one
PATH_OUPUT_NODE_TYPE_GENE =PATH_OUTPUT +  r"node_type\nodetype_gene_id.tsv"
PATH_OUPUT_NODE_TYPE_ORPHA =PATH_OUTPUT +  r"node_type\nodetype_DisorderType.tsv"
PATH_OUPUT_NODE_TYPE_PHENO =PATH_OUTPUT +  r"node_type\nodetype_PhenoType.tsv"
######################################

