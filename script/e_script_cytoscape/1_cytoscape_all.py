import pandas as pd
from multiprocessing import Pool, Process
import datetime
from time import perf_counter, sleep

from threading import Thread
from script.path_variable import *
import os

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


def stepA(df_stepOI):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""

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

        all_interractions.add((onePHENO, oneORPHA, onescore))

        all_interractions.add((onePHENO, onegenePHENO, 'g'))

        all_interractions.add((oneORPHA, onegeneORPHA, 'g'))
        all_interractions.add((oneORPHA, oneORPHAchild, 'o'))
        all_interractions.add((oneORPHAchild, onegeneORPHAchild, 'g'))

    return all_interractions


def stepB1(df_stepOI):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""

    # convert input df into a dict each values represent a row
    df_transformed_to_dict = df_stepOI.to_dict('index')

    all_interractions = set()
    for value in df_transformed_to_dict.values():
        onePHENO = value['phenopacket']
        onegenePHENO = value['gene_P']

        onescore = value['score']

        onecase = value['case']
        onegenecase = value['gene_C']

        all_interractions.add((onePHENO, onegenePHENO, 'g'))
        all_interractions.add((onecase, onegenecase, 'g'))
        all_interractions.add((onePHENO, onecase, onescore))

    return all_interractions


def stepB2(df_stepOI):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""

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

            all_interractions.add((onePHENO, onegenePHENO, 'g'))

            all_interractions.add((onecase, onegenecase, 'g'))
            all_interractions.add((oneORPHA, onegeneORPHA, 'g'))
            all_interractions.add((oneORPHA, oneORPHAchild, 'o'))
            all_interractions.add((oneORPHAchild, onegeneORPHAchild, 'g'))

            all_interractions.add((onecase, oneORPHA, onescoreB2))
            all_interractions.add((onePHENO, onecase, onescoreB1))

    return all_interractions



# def stepB2_only(df_stepOI):
#     """ get the input df convert into a dict and build the tuple for the df cytoscape output"""
#
#     # convert input df into a dict each values represent a row
#     df_transformed_to_dict = df_stepOI.to_dict('index')
#
#     all_interractions = set()
#     for value in df_transformed_to_dict.values():
#             onecase = value['case']
#             oneORPHA = value['ORPHAcode']
#             onegeneORPHA = value['symbol']
#             onescoreB2 = value['score_B2']
#             oneORPHAchild = value['ORPHAcode_child']
#             onegeneORPHAchild = value['gene_child']
#             onegenecase = value['gene_C']
#
#
#             all_interractions.add((onecase, onegenecase, 'g'))
#             all_interractions.add((oneORPHA, onegeneORPHA, 'g'))
#             all_interractions.add((oneORPHA, oneORPHAchild, 'o'))
#             all_interractions.add((oneORPHAchild, onegeneORPHAchild, 'g'))
#             all_interractions.add((onecase, oneORPHA, onescoreB2))
#
#
#     return all_interractions
#


def stepC1(df_stepOI):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""

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

        all_interractions.add((onePHENO, onegenePHENO, 'g'))
        all_interractions.add((oneORPHA, onegeneORPHA, 'g'))
        all_interractions.add((oneORPHAchild, onegeneORPHAchild, 'g'))
        all_interractions.add((onecase, onegenecase, 'g'))
        all_interractions.add((oneORPHA, oneORPHAchild, 'o'))

        all_interractions.add((onePHENO, oneORPHA, onescoreA1))
        all_interractions.add((onecase, oneORPHA, onescoreC1))

    return all_interractions


def stepC2(df_stepOI):
    """ get the input df convert into a dict and build the tuple for the df cytoscape output"""

    # convert input df into a dict each values represent a row
    df_transformed_to_dict = df_stepOI.to_dict('index')

    all_interractions = set()
    for value in df_transformed_to_dict.values():
        # print(value)
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
        all_interractions.add((onePHENO, oneORPHA, onescoreA1))
        # gene part
        all_interractions.add((onePHENO, onegenePHENO, 'g'))
        all_interractions.add((oneORPHA, onegeneORPHA, 'g'))
        all_interractions.add((oneORPHAchild, onegeneORPHAchild, 'g'))
        # orpha parent and child interaction
        all_interractions.add((oneORPHA, oneORPHAchild, 'o'))

        # issu de C1
        all_interractions.add((oneORPHA, onecase, onescoreC1))
        # gene part
        all_interractions.add((onecase, onegenecase, 'g'))

        # issu de C2
        all_interractions.add((onecase, oneORPHAcode_C2, onescore_C2))
        # gene part
        all_interractions.add((oneORPHAcode_C2, oneGene_orpha_C2, 'g'))
        all_interractions.add((oneORPHAcode_child_C2, oneGene_child_C2, 'g'))
        # orpha parent and child interaction
        all_interractions.add((oneORPHAcode_C2, oneORPHAcode_child_C2, 'o'))

    return all_interractions


# def stepC2_only(df_stepOI):
#     """ get the input df convert into a dict and build the tuple for the df cytoscape output"""
#
#     # convert input df into a dict each values represent a row
#     df_transformed_to_dict = df_stepOI.to_dict('index')
#
#     all_interractions = set()
#     for value in df_transformed_to_dict.values():
#
#
#         onecase = value['case']
#         onegenecase = value['gene_C']
#
#
#         oneORPHAcode_C2 = value['ORPHAcode_C2']
#         oneGene_orpha_C2 = value['gene_orpha_C2']
#         onescore_C2 = value['score_C2']
#         oneGene_child_C2 = value['gene_child_C2']
#         oneORPHAcode_child_C2 = value['ORPHAcode_child_C2']
#
#         # gene part
#         all_interractions.add((onecase, onegenecase, 'g'))
#
#         # issu de C2
#         all_interractions.add((onecase, oneORPHAcode_C2, onescore_C2))
#         # gene part
#         all_interractions.add((oneORPHAcode_C2, oneGene_orpha_C2, 'g'))
#         all_interractions.add((oneORPHAcode_child_C2, oneGene_child_C2, 'g'))
#         # orpha parent and child interaction
#         all_interractions.add((oneORPHAcode_C2, oneORPHAcode_child_C2, 'o'))
#
#     return all_interractions


def step_loop(df_step, letter_step):
    """ Main function call the step[letter_step] function get its tuple output build the df which is the cytoscape output """

    if letter_step == "A2":
        all_cyto_interractions = stepA(df_step)
    elif letter_step == "A1":
        all_cyto_interractions = stepA(df_step)
    elif letter_step == "B1":
        all_cyto_interractions = stepB1(df_step)
    elif letter_step == "B2":
        all_cyto_interractions = stepB2(df_step)
    # elif letter_step == "B2_only":
    #     all_cyto_interractions = stepB2_only(df_step)

    elif letter_step == 'C1':
        all_cyto_interractions = stepC1(df_step)
    elif letter_step == 'C2':
        all_cyto_interractions = stepC2(df_step)
    # elif letter_step == 'C2_only':
    #     all_cyto_interractions = stepC2_only(df_step)
    else:
        print('can t build df_all_interractions')

    # convert the list of tuple into a df
    df_all_interractions = pd.DataFrame(all_cyto_interractions, columns=["node1", 'node2', 'score'])

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


# def step_info(step,path):
#     """ get all cases files input for the step of interest """
#     list_all_pheno_step = []
#     if ((step == "A1") or (step =="C1)" or (step=="C2")or (step=="C2_only"))):
#         df_step_info = pd.read_csv( path + "nb_case_A1.tsv", sep='\t')
#         # list_all_pheno_step = set(df_step_info.iloc[:, 0])
#         list_all_pheno_step = set(df_step_info['phenopacket'])
#     elif ((step == "A2") or (step =="B1)" or (step=="B2")or (step=="B2_only"))):
#         df_step_info = pd.read_csv( path + "nb_case_B1.tsv", sep='\t')
#
#         # list_all_pheno_step = set(df_step_info.iloc[:, 0])
#         list_all_pheno_step = set(df_step_info['phenopacket'])
#
#     return list(list_all_pheno_step)


def step_info(step,path):
    """  store files names in a list"""
    os_files = os.listdir(path)

    list_os_files = []
    for filename in os_files:
        if step in filename:
            list_os_files.append(filename)

    return list_os_files


if __name__ == '__main__':
    print("###############\nSTART 1_cytoscape_all.py\n###############\n")

    list_p = []

    tab_step = ["A1","A2",'B1','B2','C1','C2']
    # # tab_step = ["A1","A2",'B1','B2_only','C1','C2_only']
    # # tab_step = ['B2_only','C2_only']
    # tab_step = ['C2','B2']

    # tab_step = ['A2']

    print("\nSTART : ", datetime.datetime.utcnow())
    for onestep in tab_step:
        # fonction which contain all phenopacket available for the step of interest
        list_all_pheno = step_info(onestep, PATH_OUTPUT_TSV)
        #print(list_all_pheno)

        for i in range(len(list_all_pheno)):
            # print("\nSTEP : ", onestep, " FOR THE CASE : ", list_all_pheno[i])

            # make the name of the input file example :P0000048_B2.tsv
            # build_fname = str( list_all_pheno[i] + "_" + onestep + ".tsv")
            build_fname = str( list_all_pheno[i] )

            # import input df
            df_step = pd.read_csv(PATH_OUTPUT_TSV + "\\" + build_fname, sep='\t')
            # df_step = df_step.iloc[:20000,:]
            msg_dfinput_size = "DF INPUT size:  {}".format(len(df_step))

            tstart = perf_counter()
            df_result_multithreads = run_multithreads(df_step, idx=onestep)
            # msg_df_size = "DF OUTPUT size:  : {}".format(len(df_result_multithreads))
            tstop = perf_counter()
            # msg_thr = "THREADED BATCH : {}".format(tstop - tstart)

            #print("EXPORT DF cytoscape OUTPUT : ", list_all_pheno[i], " FOR THE CASE : ", onestep)
            df_result_multithreads.to_csv(PATH_OUTPUT_CYTOSCAPE + "cytoscape_" + list_all_pheno[i] ,sep='\t')

            # print(msg_dfinput_size)
            # print(msg_df_size)
            # print(msg_thr)
            print("TIME: ", datetime.datetime.utcnow(), '\tSTEP ',onestep,'\tPHENOPACKET ', list_all_pheno[i],'\tDONE')


print("END : ",datetime.datetime.utcnow())
print("###############\nEND 1_cytoscape_all.py\n###############\n")

