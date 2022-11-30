import pandas as pd


from threading import Thread
from SolveRD.script.path_variable import *
import os

import logging
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()


logger.info("START\t1_cytoscape_all.py\n ")
# class for the multithreading
class T(Thread):
    def __init__(self, df_step, letter_step):
        super().__init__()  # super init to initiate the constructeur of the parent class
        self.df_step = df_step
        self.letter_step = letter_step
        self.df_final = None  # output this will contain the output cytoscape df

    def run(self):
        """ this is the function that will be exe on multithrading"""
        self.df_final = step_loop(self.df_step, self.letter_step)

def get_WK_node(path_input,path_nodetype):
    BIG_gmt_wkname = BIG_gmt_gene = []
    with open(path_input, 'r') as f:
        gmt_files = f.readlines()

        for onegmtline in gmt_files:
            gmt_files_splited = onegmtline.split('\t')
            # ['GDNF signaling%WikiPathways_20221010%WP5143%Homo sapiens',
            #  'http://www.wikipathways.org/instance/WP5143_r120629',
            #  '6482',
            #  '3576',
            #  '2331\n']
            gmt_wk = gmt_files_splited[0] # extract wk info
            gmt_wkname = gmt_wk.split('%')[0] # extract wk name
            # pslit by % and keep name only
            # ['Sphingolipid metabolism: integrated pathway',
            #  'WikiPathways_20221010',
            #  'WP4726',
            #  'Homo sapiens']
            gmt_gene = gmt_files_splited[2:len(gmt_files_splited)] # extract all gene links to it
            gmt_gene = [int(oneg.strip()) for oneg in gmt_gene] # remove the /n
            BIG_gmt_gene= BIG_gmt_gene + gmt_gene # it s a list of list
            list_gtm_wkname = [gmt_wkname]* len(gmt_gene) # create liste of duplicate wkname depending on the len of the gmt_gene list
            BIG_gmt_wkname = BIG_gmt_wkname + list_gtm_wkname


        df_gtm = pd.DataFrame({'wk_name' : BIG_gmt_wkname,
                             'id' :BIG_gmt_gene})


    df_NT_gene = pd.read_csv(path_nodetype,sep='\t')
    df_NT_gene.columns = ['symbol', 'id_orphanet','id']

    df_gmt_rslt_o  = pd.merge(df_NT_gene,df_gtm, on='id',how='outer' )
    df_gmt_rslt_oNA = df_gmt_rslt_o.dropna()
    # df_gmt_rslt_oNA['id_orphanet'] = df_gmt_rslt_oNA['id_orphanet'].astype(int)


    return df_gmt_rslt_oNA


def get_wk_interaction(dfwk,gene,wk_interaction):
    list_wk = dfwk[dfwk['symbol'] == gene]['wk_name']
    for one_path in list_wk:
        wk_interaction.add((gene, one_path,"wikipathways",'w'))

    return wk_interaction


def stepA(df_stepOI,df_WK):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""
    wk_symbol = df_WK['symbol'].tolist()

    # convert input df into a dict each values represent a row
    df_transformed_to_dict = df_stepOI.to_dict('index')
    all_interractions = set()
    for value in df_transformed_to_dict.values():
        onePHENO = value['phenopacket']
        onegenePHENO = value['gene_P']

        oneORPHA = value['ORPHAcode']
        onegeneORPHA = value['symbol']

        oneORPHAchild = value['ORPHAcode_child']
        onegeneORPHAchild = value['gene_child']

        onescore = value['score']

        all_interractions.add((onePHENO, oneORPHA,'resnik', onescore))

        all_interractions.add((onePHENO, onegenePHENO,'resnik', 'g'))

        all_interractions.add((oneORPHA, onegeneORPHA,'resnik', 'g'))
        all_interractions.add((oneORPHA, oneORPHAchild, 'resnik','o'))
        all_interractions.add((oneORPHAchild, onegeneORPHAchild,'resnik', 'g'))

        if (onegeneORPHA in wk_symbol ):
            all_interractions = get_wk_interaction(df_WK, onegeneORPHA, all_interractions)
        elif (onegeneORPHAchild in wk_symbol) :
            all_interractions = get_wk_interaction(df_WK, onegeneORPHAchild, all_interractions)
        elif  (onegenePHENO in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegenePHENO, all_interractions)


    return all_interractions


def stepB1(df_stepOI,df_WK):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""
    wk_symbol = df_WK['symbol'].tolist()

    # convert input df into a dict each values represent a row
    df_transformed_to_dict = df_stepOI.to_dict('index')

    all_interractions = set()
    for value in df_transformed_to_dict.values():
        onePHENO = value['phenopacket']
        onegenePHENO = value['gene_P']

        onescore = value['score']

        onecase = value['case']
        onegenecase = value['gene_C']

        all_interractions.add((onePHENO, onegenePHENO,'resnik', 'g'))
        all_interractions.add((onecase, onegenecase,'resnik', 'g'))
        all_interractions.add((onePHENO, onecase,'resnik', onescore))

        if (onegenePHENO in wk_symbol ):
            all_interractions = get_wk_interaction(df_WK, onegenePHENO, all_interractions)
        elif (onegenecase in wk_symbol) :
            all_interractions = get_wk_interaction(df_WK, onegenecase, all_interractions)


    return all_interractions


def stepB2(df_stepOI,df_WK):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""
    wk_symbol = df_WK['symbol'].tolist()

    # convert input df into a dict each values represent a row
    df_transformed_to_dict = df_stepOI.to_dict('index')

    all_interractions = set()
    for value in df_transformed_to_dict.values():

            onePHENO = value['phenopacket']
            onegenePHENO = value['gene_P']

            oneORPHA = value['ORPHAcode']
            onegeneORPHA = value['symbol']

            oneORPHAchild = value['ORPHAcode_child']
            onegeneORPHAchild = value['gene_child']

            onescoreB1 = value['score_B1']
            onescoreB2 = value['score_B2']

            onecase = value['case']
            onegenecase = value['gene_C']

            all_interractions.add((onePHENO, onegenePHENO,'resnik', 'g'))

            all_interractions.add((onecase, onegenecase,'resnik', 'g'))
            all_interractions.add((oneORPHA, onegeneORPHA,'resnik', 'g'))
            all_interractions.add((oneORPHA, oneORPHAchild,'resnik', 'o'))
            all_interractions.add((oneORPHAchild, onegeneORPHAchild,'resnik', 'g'))

            all_interractions.add((onecase, oneORPHA, onescoreB2))
            all_interractions.add((onePHENO, onecase, onescoreB1))

            if (onegenePHENO in wk_symbol):
                all_interractions = get_wk_interaction(df_WK, onegenePHENO, all_interractions)
            elif (onegenecase in wk_symbol):
                all_interractions = get_wk_interaction(df_WK, onegenecase, all_interractions)
            elif (onegeneORPHA in wk_symbol):
                all_interractions = get_wk_interaction(df_WK, onegeneORPHA, all_interractions)
            elif (onegeneORPHAchild in wk_symbol):
                all_interractions = get_wk_interaction(df_WK, onegeneORPHAchild, all_interractions)

    return all_interractions



def stepC1(df_stepOI,df_WK):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""
    wk_symbol = df_WK['symbol'].tolist()

    # convert input df into a dict each values represent a row
    df_transformed_to_dict = df_stepOI.to_dict('index')

    all_interractions = set()
    for value in df_transformed_to_dict.values():
        onePHENO = value['phenopacket']
        onegenePHENO = value['gene_P']

        oneORPHA = value['ORPHAcode']
        onegeneORPHA = value['symbol']

        onescoreA1 = value['score_A1']
        onescoreC1 = value['score_C1']

        oneORPHAchild = value['ORPHAcode_child']
        onegeneORPHAchild = value['gene_child']

        onecase = value['case']
        onegenecase = value['gene_C']

        all_interractions.add((onePHENO, onegenePHENO,'resnik', 'g'))
        all_interractions.add((oneORPHA, onegeneORPHA,'resnik', 'g'))
        all_interractions.add((oneORPHAchild, onegeneORPHAchild,'resnik', 'g'))
        all_interractions.add((onecase, onegenecase,'resnik', 'g'))
        all_interractions.add((oneORPHA, oneORPHAchild,'resnik','o'))

        all_interractions.add((onePHENO, oneORPHA,'resnik', onescoreA1))
        all_interractions.add((onecase, oneORPHA,'resnik', onescoreC1))

        if (onegenePHENO in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegenePHENO, all_interractions)
        elif (onegeneORPHA in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegeneORPHA, all_interractions)
        elif (onegeneORPHAchild in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegeneORPHAchild, all_interractions)
        elif (onegenecase in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegenecase, all_interractions)

    return all_interractions


def stepC2(df_stepOI,df_WK):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""
    wk_symbol = df_WK['symbol'].tolist()

    # convert input df into a dict each values represent a row
    df_transformed_to_dict = df_stepOI.to_dict('index')

    all_interractions = set()
    for value in df_transformed_to_dict.values():
        onePHENO = value['phenopacket']
        onegenePHENO = value['gene_P']

        onescoreA1 = value['score_A1']
        onescoreC1 = value['score_C1']

        onecase = value['case']
        onegenecase = value['gene_C']
        oneORPHA = value['ORPHAcode']
        onegeneORPHA = value['symbol']
        oneORPHAchild = value['ORPHAcode_child']
        onegeneORPHAchild = value['gene_child']

        oneORPHAcode_C2 = value['ORPHAcode_C2']
        oneGene_orpha_C2 = value['gene_orpha_C2']
        onescore_C2 = value['score_C2']
        oneGene_child_C2 = value['gene_child_C2']
        oneORPHAcode_child_C2 = value['ORPHAcode_child_C2']

        # gene interactions
        # issu de A1
        all_interractions.add((onePHENO, oneORPHA,'resnik', onescoreA1))
        # gene part
        all_interractions.add((onePHENO, onegenePHENO,'resnik', 'g'))
        all_interractions.add((oneORPHA, onegeneORPHA,'resnik', 'g'))
        all_interractions.add((oneORPHAchild, onegeneORPHAchild,'resnik', 'g'))
        # orpha parent and child interaction
        all_interractions.add((oneORPHA, oneORPHAchild,'resnik', 'o'))

        # issu de C1
        all_interractions.add((oneORPHA, onecase,'resnik', onescoreC1))
        # gene part
        all_interractions.add((onecase, onegenecase,'resnik', 'g'))

        # issu de C2
        all_interractions.add((onecase, oneORPHAcode_C2,'resnik', onescore_C2))
        # gene part
        all_interractions.add((oneORPHAcode_C2, oneGene_orpha_C2,'resnik', 'g'))
        all_interractions.add((oneORPHAcode_child_C2, oneGene_child_C2,'resnik', 'g'))
        # orpha parent and child interaction
        all_interractions.add((oneORPHAcode_C2, oneORPHAcode_child_C2,'resnik', 'o'))


        if (onegenePHENO in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegenePHENO, all_interractions)
        elif (onegeneORPHA in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegeneORPHA, all_interractions)
        elif (onegeneORPHAchild in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegeneORPHAchild, all_interractions)
        elif (onegenecase in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, onegenecase, all_interractions)
        elif (oneGene_orpha_C2 in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, oneGene_orpha_C2, all_interractions)
        elif (oneGene_child_C2 in wk_symbol):
            all_interractions = get_wk_interaction(df_WK, oneGene_child_C2, all_interractions)


    return all_interractions


def step_loop(df_step, letter_step):
    """ Main function call the step[letter_step] function get its tuple output build the df which is the cytoscape output """
    df_WK_gene = get_WK_node(PATH_INPUT_WK, PATH_OUPUT_NODE_TYPE_GENE)

    if letter_step == "A2":
        all_cyto_interractions = stepA(df_step,df_WK_gene)
    elif letter_step == "A1":
        all_cyto_interractions = stepA(df_step,df_WK_gene)
    elif letter_step == "B1":
        all_cyto_interractions = stepB1(df_step,df_WK_gene)
    elif letter_step == "B2":
        all_cyto_interractions = stepB2(df_step,df_WK_gene)

    elif letter_step == 'C1':
        all_cyto_interractions = stepC1(df_step,df_WK_gene)
    elif letter_step == 'C2':
        all_cyto_interractions = stepC2(df_step,df_WK_gene)

    else:
        logger.info('1_cytoscape_all.py\tcan t build df_all_interractions')

    # convert the list of tuple into a df
    df_all_interractions = pd.DataFrame(all_cyto_interractions, columns=["node1", 'node2','type_interaction','score'])

    # FILTRATION
    # remove row depending on orpha:nan is was for orpha child remove this because it s a fake interaction
    df_all_interractions = df_all_interractions.loc[df_all_interractions["node2"] != "ORPHA:nan"]
    df_all_interractions = df_all_interractions.loc[df_all_interractions["node1"] != "ORPHA:nan"]
    # remove Nan because we don't need node named NaN
    df_all_interractions = df_all_interractions.dropna()

    return df_all_interractions


def run_multithreads(df_step, idx):
    """ create subset of the df build each row thanks to multithreading and concat each subset in a final one"""
    list_thread = []  # store all created threads.

    for i in range(0, len(df_step), 1000):  # clive the df into 1000 row subset row
        df_subset = df_step.iloc[i:i + 1000, :]  # subset of the df wich contain 1000 row
        threadclass = T(df_subset, idx)  # instantiate the classe
        threadclass.start()  # start/run the thread
        list_thread.append(threadclass)  # add the thread in a list

    list_result_thread = []  # contain result thread
    for threadclass in list_thread:
        threadclass.join()  # wait the end of the thread job !
        result = threadclass.df_final  # get value
        list_result_thread.append(result)  # contain all values of thread

    df_result_all_df = pd.concat(list_result_thread)  # concat and build the output df
    df_result_all_df = df_result_all_df.drop_duplicates()  # remove doublons

    return df_result_all_df



def step_info(step,path):
    """  store files names in a list"""
    os_files = os.listdir(path)

    list_os_files = []
    for filename in os_files:
        if step in filename:
            list_os_files.append(filename)

    return list_os_files


if __name__ == '__main__':

    tab_step = ["A1","A2",'B1','B2','C1','C2']
    # tab_step = ["A2"]
    # tab_step = ['C2','B2']

    # tab_step = ['A2']

    for onestep in tab_step:
        # fonction which contain all phenopacket available for the step of interest
        list_all_pheno = step_info(onestep, PATH_OUTPUT_TSV)

        for i in range(len(list_all_pheno)):

            # make the name of the input file example :P0000048_B2.tsv
            # build_fname = str( list_all_pheno[i] + "_" + onestep + ".tsv")
            build_fname = str( list_all_pheno[i] )

            # import input df
            df_step = pd.read_csv(PATH_OUTPUT_TSV + "/" + build_fname, sep='\t')
            # df_step = df_step.iloc[:20000,:]

            df_result_multithreads = run_multithreads(df_step, idx=onestep)


            # EXPORT DF cytoscape OUTPUT FOR THE CASE
            df_result_multithreads.to_csv(PATH_OUTPUT_CYTOSCAPE + "cytoscape_" + list_all_pheno[i] ,sep='\t')


            logger.info('1_cytoscape_all.py\tSTEP {}\tPHENOPACKET  {}\tDONE'.format(onestep,list_all_pheno[i]))


logger.info("END\t1_cytoscape_all.py\n ")

