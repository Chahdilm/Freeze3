import subprocess
import sys
from multiprocessing import Process
import time

import datetime
import os

from script.path_variable import *

def script_step(name_file, path):
    result = subprocess.run([sys.executable, path + name_file], capture_output=True, text=True)
    print("stdout:", result.stdout)
    print("stderr:", result.stderr)

def create_dir(folder_name):
    try:
        # Create target Directory
        os.mkdir(folder_name)
        print("Directory " , folder_name ,  " Created ")
    except FileExistsError:
        print("Directory " , folder_name ,  " already exists")


if __name__ == '__main__':
    start = time.process_time()

    create_dir(PATH + "output_5HPO")

    create_dir(PATH_OUTPUT + "curation_input")
    create_dir(PATH_OUTPUT + "curation_input/Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete")
    create_dir(
        PATH_OUTPUT + "curation_input/Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete")

    create_dir(PATH_OUTPUT + "gene_info")

    create_dir(PATH_OUTPUT + "tsv")

    create_dir(PATH_OUTPUT + "cytoscape")
    create_dir(PATH_OUTPUT + "cytoscape_wikipathways")
    create_dir(PATH_OUTPUT + "node_type")

    create_dir(PATH_OUTPUT + "json")
    create_dir(PATH_OUTPUT + "variant")
    create_dir(PATH_OUTPUT + "autocomplete")



    os_files = os.listdir(PATH_SCRIPT)
    files_curation = os.listdir(PATH_SCRIPT+'/'+"a_script_curation") # folder a_script_curation
    files_curation.remove('__init__.py')

    files_gene = os.listdir(PATH_SCRIPT+'/'+"b_script_get_gene") # folder b_script_get_gene
    files_gene.remove('__init__.py')

    files_all_step = os.listdir(PATH_SCRIPT+'/'+"c_script_step") # folder c_script_step
    files_all_step.remove('__init__.py')

    files_minitsv = os.listdir(PATH_SCRIPT+'/'+"d_script_minitsv") # folder d_script_minitsv
    files_minitsv.remove('__init__.py')

    files_cytoscape = os.listdir(PATH_SCRIPT+'/'+"e_script_cytoscape") # folder e_script_cytoscape
    files_cytoscape.remove('__init__.py')

    files_json = os.listdir(PATH_SCRIPT+'/'+"f_script_json") # folder f_script_json
    files_cytoscape.remove('__init__.py')

    files_hp = os.listdir(PATH_SCRIPT+'/'+"g_script_homepageJS") # g_script_homepageJS
    files_hp.remove('__init__.py')


    ##################################
    ####         CURATION      #######
    ##################################
    start = time.process_time()
    list_proc = []
    for file in files_curation:
        print(file, '\t', files_curation.index(file))
        P_step = Process(target=script_step, args=(file, str(PATH_SCRIPT + "a_script_curation" + '/')))
        list_proc.append(P_step)

    list_proc[0].start()    # 1_filter_file_phenopacket
    list_proc[0].join()     # 1_filter_file_phenopacket

    list_proc[1].start()    # 2_curation.py
    list_proc[1].join()     # 2_curation.py
    print("\nTIME END\t ", datetime.datetime.utcnow())
    print(f'Process 1_filter_file_phenopacket.py  is alive: {list_proc[1].is_alive()}')

    list_proc[2].start()     # 3_transform_phenopackets_cleaned_for_RunSolveRD.py
    list_proc[2].join()     # 3_transform_phenopackets_cleaned_for_RunSolveRD.py
    print("\nTIME END\t ", datetime.datetime.utcnow())
    print(f'Process 2_curation.py is alive: {list_proc[2].is_alive()}\n\n')

    list_proc[3].start()     # 4_RUN_runSolvedRD.py
    list_proc[3].join()     # 4_RUN_runSolvedRD.py
    print("\nTIME END\t ", datetime.datetime.utcnow())
    print(f'Process 3_transform_phenopackets_cleaned_for_RunSolveRD.py is alive: {list_proc[3].is_alive()}\n\n')

    print('END a_folder time : ', time.process_time() - start)

    ###################################
    #####         GET GENES     #######
    ###################################

    start = time.process_time()
    list_proc_input = []
    for file in files_gene:
        print(file,'\t', files_gene.index(file))
        P_step = Process(target=script_step, args = (file,str(PATH_SCRIPT+"b_script_get_gene"+'/')))
        list_proc_input.append(P_step)

    for one_proc in list_proc_input:
        one_proc.start()

    for one_proc in list_proc_input:
        one_proc.join()
        print(f'Process  is alive: {one_proc.is_alive()}\n')
        print("\nTIME END\t ", datetime.datetime.utcnow())

    print('END time B_folder: ', time.process_time() - start)
    ####################################


    ###################################
    #####         STEPS         #######
    ###################################

    start = time.process_time()
    list_proc = []
    for file in files_all_step:
        print(file, '\t', files_all_step.index(file))
        P_step = Process(target=script_step, args=(file, str(PATH_SCRIPT +"c_script_step" + '/')))
        list_proc.append(P_step)


    list_proc[0].start()    # A1,A2
    list_proc[1].start()    # B1
    list_proc[0].join()     # A1,A2
    list_proc[1].join()     # B1
    print("\nTIME END\tA1\tA2\tB1: ", datetime.datetime.utcnow())
    print(f'Process A1 A2 is alive: {list_proc[1].is_alive()}')
    print(f'Process B1 is alive: {list_proc[2].is_alive()}\n\n')

    list_proc[2].start()    # B2
    list_proc[3].start()    # C1
    list_proc[2].join()     # B2
    list_proc[3].join()     # C1
    print("\nTIME END\tB2\tC1: ", datetime.datetime.utcnow())
    print(f'Process B2 is alive: {list_proc[3].is_alive()}')
    print(f'Process C1 is alive: {list_proc[4].is_alive()}\n\n')

    list_proc[4].start()     # C2
    list_proc[4].join()     # C2
    print("\nTIME END  C2: ", datetime.datetime.utcnow())
    print(f'Process C2 is alive: {list_proc[5].is_alive()}\n\n')

    print('time : ', time.process_time() - start)
    ####################################


    ##################################
    ####        MINI TSV       ######
    ##################################
    start = time.process_time()
    list_proc = []

    print(files_minitsv[0],'\t', files_minitsv.index(files_minitsv[0]))
    P_step = Process(target=script_step, args = (files_minitsv[0],str(PATH_SCRIPT+"d_script_minitsv"+'/')))
    list_proc.append(P_step)


    for one_proc in list_proc:
        one_proc.start()

    for one_proc in list_proc:
        one_proc.join()
        print(f'Process  is alive: {one_proc.is_alive()}\n')
        print("\nTIME END\t ", datetime.datetime.utcnow())
    print('time : ', time.process_time() - start)

    ###################################
    #####        CYTOSCAPE      #######
    ###################################
    start = time.process_time()
    list_proc = []
    for file in files_cytoscape:
        print(file, '\t', files_cytoscape.index(file))
        P_step = Process(target=script_step, args=(file, str(PATH_SCRIPT + "e_script_cytoscape" + '/')))
        list_proc.append(P_step)

    list_proc[0].start()  # 1_cytoscape_all
    list_proc[2].start()  # 3_get_ALL_nodetype

    list_proc[0].join()  # 1_cytoscape_all
    list_proc[2].join()  # 3_get_ALL_nodetype
    print(f'Process 0_build_output_folder is alive: {list_proc[0].is_alive()}')
    print(f'Process 1_cytoscape_all  is alive: {list_proc[1].is_alive()}\n\n')
    print(f'Process 3_get_ALL_nodetype  is alive: {list_proc[3].is_alive()}\n\n')

    list_proc[1].start()  # 2_wikipathway_all_py4cy
    list_proc[1].join()  # 2_wikipathway_all_py4cy

    print(f'Process 2_wikipathway_all_py4cy  is alive: {list_proc[2].is_alive()}\n\n')
    print("\n\nTIME END\t: ", datetime.datetime.utcnow())


    ####################################
    ######        JSON      #######
    ####################################
    start = time.process_time()
    list_proc_input = []
    for file in files_json:
        print(file,'\t', files_json.index(file))
        P_step = Process(target=script_step, args = (file,str(PATH_SCRIPT+"f_script_json"+'/')))
        list_proc_input.append(P_step)

    for one_proc in list_proc_input:
        one_proc.start()

    for one_proc in list_proc_input:
        one_proc.join()
        print(f'Process  is alive: {one_proc.is_alive()}\n')

    print('time : ', time.process_time() - start)

    ###################################
    ##### MINI TSV  AND HOMEPAGE ######
    ###################################
    start = time.process_time()
    list_proc = []


    print(files_hp[0],'\t', files_hp.index(files_hp[0]))
    P_step = Process(target=script_step, args = (files_hp[0],str(PATH_SCRIPT+"g_script_homepageJS"+'/')))
    list_proc.append(P_step)

    for one_proc in list_proc:
        one_proc.start()

    for one_proc in list_proc:
        one_proc.join()
        print(f'Process  is alive: {one_proc.is_alive()}\n')
        print("\nTIME END\t ", datetime.datetime.utcnow())
    print('time : ', time.process_time() - start)