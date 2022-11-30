
import pandas as pd
from SolveRD.script.path_variable import *

from multiprocessing import  Process


import logging
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logger = logging.getLogger()

logger.info("START\tbuild_minidf.py\n ")


def build_mini(df_step,step,path,list_phenopacket):
    for onepheno in list_phenopacket:
        mini_df = df_step[(df_step.phenopacket.values ==onepheno)]

        mini_df.to_csv(path+str(onepheno)+"_"+str(step)+".tsv", sep='\t',index=False)

def build_mini_df_each_pheno(df_step,step,path):
    # get all phenopacket
    list_phenopacket = df_step["phenopacket"].drop_duplicates().tolist()

    build_mini(df_step, step, path, list_phenopacket)

    logger.info("build_minidf.py\t{}\tEXPORT EACH CASE\tNB : {}".format( str(step),len(list_phenopacket)))


if __name__ == '__main__':

    # df_stepA1 = pd.read_csv(PATH_OUTPUT+r"stepA1.tsv", sep='\t')
    # df_stepA2 = pd.read_csv(PATH_OUTPUT+r"stepA2.tsv", sep='\t')
    # df_stepB1 = pd.read_csv(PATH_OUTPUT+r"stepB1.tsv", sep='\t')
    # df_stepB2 = pd.read_csv(PATH_OUTPUT+r"stepB2.tsv", sep='\t') #stepC2_only  stepB2
    # df_stepC1 = pd.read_csv(PATH_OUTPUT+r"stepC1.tsv", sep='\t')
    # df_stepC2 = pd.read_csv(PATH_OUTPUT+r"stepC2.tsv", sep='\t') # stepC2_only

    # build_mini(df_stepA2,'A2',PATH_OUTPUT_TSV)
    # build_mini(df_stepA1,'A1',PATH_OUTPUT_TSV)
    # build_mini(df_stepB1,'B1',PATH_OUTPUT_TSV)
    # build_mini(df_stepB2,'B2',PATH_OUTPUT_TSV)
    # build_mini(df_stepC1,'C1',PATH_OUTPUT_TSV)
    # build_mini(df_stepC2,'C2',PATH_OUTPUT_TSV)


    tab_step = ["A1","A2",'B1','B2','C1','C2']
    # tab_step = ["A1","A2",'B1']
    list_proc = []
    for i in range(len(tab_step)):
        logger.info("build_minidf.py\t {}".format(tab_step[i]))
        P_step = Process(target=build_mini_df_each_pheno, args=((pd.read_csv(PATH_OUTPUT+"step"+tab_step[i]+".tsv",sep='\t'),(tab_step[i]),(PATH_OUTPUT_TSV))))
        list_proc.append(P_step)



    for onejob in list_proc:
        onejob.start()

    for onejob in list_proc:
        onejob.join()

    logger.info("END\tbuild_minidf.py\n ")







