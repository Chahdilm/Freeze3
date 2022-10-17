

from time import perf_counter
import pandas as pd
from script.path_variable import *

# from multiprocessing import Queue, Process, Pool
import time
import datetime

def build_mini(df_step,step,path):
    phenopacket =  df_step.iloc[:, 0]
    # time start
    tstart = perf_counter()

    # this loop will export a df for each case of interest (named phenopacket here)
    for onepheno in phenopacket:

        mini_df = df_step[df_step.iloc[:, 0] == str(onepheno) ]
        mini_df.to_csv(path+str(onepheno)+"_"+str(step)+".tsv", sep='\t',index=False)

    # time end
    tstop = perf_counter()
    print(str(step)+"\tTIME build all df : ",(tstop - tstart))
    print(str(step)+"\tEXPORT EACH CASE\tNB :  ", len(phenopacket),"\tDONE")


if __name__ == '__main__':

    start = perf_counter()
    df_stepA1 = pd.read_csv(PATH_OUTPUT+r"stepA1.tsv", sep='\t')

    df_stepA2 = pd.read_csv(PATH_OUTPUT+r"stepA2.tsv", sep='\t')
    df_stepB1 = pd.read_csv(PATH_OUTPUT+r"stepB1.tsv", sep='\t')

    df_stepB2 = pd.read_csv(PATH_OUTPUT+r"stepB2.tsv", sep='\t') #stepC2_only  stepB2

    df_stepC1 = pd.read_csv(PATH_OUTPUT+r"stepC1.tsv", sep='\t')
    df_stepC2 = pd.read_csv(PATH_OUTPUT+r"stepC2.tsv", sep='\t') # stepC2_only




    build_mini(df_stepA2,'A2',PATH_OUTPUT_TSV)

    build_mini(df_stepA1,'A1',PATH_OUTPUT_TSV)

    build_mini(df_stepB1,'B1',PATH_OUTPUT_TSV)

    build_mini(df_stepB2,'B2',PATH_OUTPUT_TSV)

    build_mini(df_stepC1,'C1',PATH_OUTPUT_TSV)

    build_mini(df_stepC2,'C2',PATH_OUTPUT_TSV)

    print('END time B_folder: ', time.process_time() - start)

    # tab_step = ["A1","A2",'B1','B2','C1','C2']
    # tab_step = ["A1"]
    # list_proc = []
    # for i in range(len(tab_step)):
    #     print("build_minidf.py\t",tab_step[i],'\t')
    #     P_step = Process(target=build_mini, args=((pd.read_csv(PATH_OUTPUT+r"step"+tab_step[i]+".tsv",sep='\t'),(tab_step[i]),(PATH_OUTPUT_TSV))))
    #     list_proc.append(P_step)
    #
    #
    # for onejob in list_proc:
    #     onejob.start()
    #
    # for onejob in list_proc:
    #     onejob.join()
    #     print(f'Is Process alive ?: {onejob.is_alive()}')
    #     print("TIME END: ", datetime.datetime.utcnow())







