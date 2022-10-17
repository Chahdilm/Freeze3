# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 11:54:36 2021

@author: ohongnat
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:59:07 2021

@author: ohongnat,mchahdil
"""

import json
import os

from script.path_variable import *

print("3_transform_phenopackets_cleaned_for_RunSolveRD.py")

all_pheno = os.listdir(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION)

try:
    # Remove element from list phenopacket
    all_pheno.remove('results')
    print("3_transform_phenopackets_cleaned_for_RunSolveRD.py\tElement\tRemoved ")
except :
    print("3_transform_phenopackets_cleaned_for_RunSolveRD.py\tElement\tAlreadyRemoved ")

# Open all json files
for filename in all_pheno :
    with open(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION+filename) as jsonfile:
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

print("3_transform_phenopackets_cleaned_for_RunSolveRD.py\tRename json file")
# rename all files (keep only the id
for onefile in all_pheno:
    onefile_strip = onefile.split('.')
    os.rename(PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION+onefile,PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION+onefile_strip[0]+str('.json'))