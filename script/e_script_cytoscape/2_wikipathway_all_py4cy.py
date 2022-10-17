# use conda python 3.9.7

import os
import pandas as pd
import py4cytoscape as p4c

import numpy as np  # pour les nan replace to 0 in the score col

import datetime
from time import perf_counter

from script.path_variable import *

##############################################################################################################################


def read_osfiles(path):
    """  store files names in a list"""
    os_files = os.listdir(path)
    len(os_files)

    return list(os_files)


##############################################################################################################################


def main_build_cytoWK(onestep):
    # store all names files of one folder
    cytoscape_folder = read_osfiles(PATH_OUTPUT_CYTOSCAPE)
    print(onestep, "\tINPUT cytsocape files size : \t", len(cytoscape_folder))

    no_home_file_STEP_cytoscape_folder = []

    for cyto_file in cytoscape_folder:
        if onestep in cyto_file:
            no_home_file_STEP_cytoscape_folder.append(cyto_file)

    # no_home_file_STEP_cytoscape_folder = ['cytoscape_P0000044_B2.tsv']

    for cytoscapefile in no_home_file_STEP_cytoscape_folder:

        if onestep in cytoscapefile:

            splited = cytoscapefile.split('_')
            removetsv = str(splited[1])
            # print(onestep,"\tCYTOSCAPE file ",cytoscapefile,'\t',removetsv,'\t')

            # removetsv = "P0012693"
            # transform cytoscape_P0011778_A2.tsv" into P0011778

            # open df cytoscape file (col : node1,node2,score)
            path_cytoscape_format_file = PATH_OUTPUT_CYTOSCAPE + str(cytoscapefile)

            # store the tsv cytoscape pd
            df_cys_format = pd.read_csv(path_cytoscape_format_file, sep='\t', usecols=[1, 2, 3])

            # build a specific pd compatible for cytoscape
            edge_data = {'source': df_cys_format["node1"],
                         'target': df_cys_format["node2"],
                         'interaction': "resnik",
                         'score': df_cys_format["score"]
                         }
            edges = pd.DataFrame(data=edge_data, columns=['source', 'target', 'interaction', 'score'])

            # print(onestep,"\t",removetsv,"\tBuild df for cytoscape\tcolumns=['source', 'target', 'interaction','score']")

            # load the edges df in cytoscape !
            p4c.networks.create_network_from_data_frames(None, edges,
                                                         title='network ' + str(removetsv) + ' ' + str(onestep),
                                                         collection='collection ')
            # print(onestep,"\tTIME",datetime.datetime.utcnow(),"\tLOAD this df into cytoscape ")

            # IMPORtANT !!! launch Cytargetlinker command throught py4cytoscape !!!!

            # wikipathways-20220511-hsa-REACTOME.xgmml  Gene-RD-Prov_v2.xgmml
            p4c.commands.commands_run(
                'cytargetlinker extend idAttribute="shared name" linkSetFiles="'+PATH_INPUT_WK+'"wikipathways-20220511-hsa-REACTOME.xgmml" network=current')
            # # export image
            # p4c.notebook_export_show_image()
            # load the table network from cytoscape
            # network_table_data = p4c.tables.get_table_columns(columns=["shared name","CTL.label","CTL.PathwayName"])
            # store all nodes from cytosape in a list
            # all_nodes = p4c.networks.get_all_nodes()
            # print(onestep,"\tTIME",datetime.datetime.utcnow(),"\tOPEN wikipathway + cytoscape commands")

            # export the network in a SIF format
            SIF_network = p4c.export_network('network', 'SIF')
            # print(onestep,"\tTIME",datetime.datetime.utcnow(),"\tEXPORT df + cytoscape commands")

            # p4c.delete_network()

            # put the SIF in a df pandas
            df_sif = pd.read_csv(SIF_network['file'], sep='\t', header=None)
            # print(onestep,"\tLOAD SIF EXPORTED from cytoscape ")

            # rename df
            df_sif.columns = ['node1', 'type_interaction', 'node2']
            # replace some values in one col
            new_col = df_sif["type_interaction"].replace('-', "wikipathway", regex=True).to_list()
            # remove the wrong col
            df_sif.pop('type_interaction')
            # add the corect one
            df_sif['type_interaction'] = new_col

            # merge the two df
            df_merge = pd.merge(df_sif, df_cys_format, how='left', left_on=['node1', 'node2'],
                                right_on=['node1', 'node2'])
            # print(onestep,"\tMERGE SIF df with cytoscape df")

            # remove score col wrong and replace it by a new one same way as before for df_sif
            new_col_score = df_merge["score"].replace(np.nan, 'wk', regex=True).to_list()
            df_merge.pop('score')
            df_merge['score'] = new_col_score

            # df_merge[df_merge['node1']=="ORPHA:444099"]
            # df_merge[df_merge['node2']=="ORPHA:444099"]

            # remove the sif file because it s useless to keep it on the pc
            os.remove(SIF_network['file'])
            # print(onestep,"\tDELETE SIF")

            # p4c.delete_network()
            df_merge.to_csv(
                PATH_OUTPUT_CYTOSCAPE_WIKIPATHWAY+"cytoscape_wikipathway_" + removetsv + "_" + onestep + ".tsv",
                sep='\t', index=False)
            print(onestep, "\tTIME END  JSON: ", datetime.datetime.utcnow(), "FOR THE CASE\t", removetsv)

        else:
            print(onestep, '\tno in the step', cytoscapefile)  # jamais exe celui la apart pour les homepage

        # print(onestep,"\tDELETE ALL NETWORK on Cytoscape")
        p4c.delete_all_networks()


if __name__ == '__main__':
    print("###############\nSTART 3_wikipathway_all.py\n###############\n")

    ## CYTOSCAPE
    # install pluging
    p4c.install_app("WikiPathways")

    # open cytoscape //
    p4c.cytoscape_ping()
    p4c.cytoscape_version_info()



    tstart = perf_counter()

    main_build_cytoWK("A1")
    main_build_cytoWK("A2")
    main_build_cytoWK("B1")
    main_build_cytoWK("B2")
    main_build_cytoWK("C1")
    main_build_cytoWK("C2")

    tstop = perf_counter()
    print("TIME build all df : ", (tstop - tstart))

    # nettoyer tout les réseaux de la session
    p4c.session.close_session(False)
    print("TIME END : ", datetime.datetime.utcnow())

    print("###############\nEND 3_wikipathway_all.py\n###############\n")

