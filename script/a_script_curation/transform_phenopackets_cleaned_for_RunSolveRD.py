# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:59:07 2021

@author: ohongnat,mchahdil
"""

import json
import os

from SolveRD.script.path_variable import *

import logging
# logging.basicConfig(level=logging.INFO,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')
logging.basicConfig(filename=PATH_INIT+'/project.log',  level=logging.DEBUG,format='%(asctime)-10s:%(levelname)-20s:%(name)s:%(message)-20s')

logger = logging.getLogger()

logger.info("START\t3_transform_phenopackets_cleaned_for_RunSolveRD\n ")

all_pheno = os.listdir(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION_tmp)

try:
    # Remove element from list phenopacket
    all_pheno.remove('results')
    logger.info("transform_phenopackets_cleaned_for_RunSolveRD.py\tElement\tRemoved ")
except :
    logger.info("transform_phenopackets_cleaned_for_RunSolveRD.py\tElement\tAlreadyRemoved ")

# Open all json files
for filename in all_pheno :
    with open(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION_tmp+filename) as jsonfile:
        data = json.load(jsonfile)
        id = data['id']
        phenotypes = data['phenotypes']
        # condition to avoid the runsolvedRD error : "Expect an array but found: {}"
        if not data['genes']:
            # print(id)
            data['genes']=[]
        genes = data['genes']

    # We create the dico & lists for the output
    data_output = {"id": "",
                   "phenotypes": [],
                   "genes": []}

    # We rename the keys id & label
    phenotypes_output = []
    for hpo in phenotypes:
        dict_json = {"type": {"id": "", "label": ""}}
        dict_json["type"]["id"] = hpo['type.id']
        dict_json["type"]["label"] = hpo['type.label']
        phenotypes_output.append(dict_json)

    # We put these into the output list
    data_output['phenotypes'] = phenotypes_output

    # We prepare the output
    data_output['id'] = id
    data_output['genes'] = genes

    # We export the outputs

    with open(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION+filename, mode='w') as jsonfileoutput:
        json.dump(data_output, jsonfileoutput)
        jsonfileoutput.close()
        jsonfile.close()

logger.info("transform_phenopackets_cleaned_for_RunSolveRD.py\tRe-structured all json files for algo")


# rename all files (keep only the id
for onefile in all_pheno:
    onefile_strip = onefile.split('.')
    os.rename(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION+onefile,PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION+onefile_strip[0]+str('.json'))

logger.info("transform_phenopackets_cleaned_for_RunSolveRD.py\tRename json file")
logger.info("END\t3_transform_phenopackets_cleaned_for_RunSolveRD\n ")

