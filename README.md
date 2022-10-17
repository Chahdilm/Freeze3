### Auteur
Maroua CHAHDIL  

# Goal
[...] **WIP**


# Install : 
- Réalisation du projet **sous Mac OS** Catalina Version 10.15.7 (19H15)

  - R 4.2.1
use R Language for IntelliJ on pycharm (223.6646—223.6646.*)
  - Python 3.9

# Package needed  :
## Python
- **Pandas** : 
  - install librairie on python console : use conda or pip
```
pandas 
os 
json
subprocess
datetime
time
multiprocessing
threading
numpy
py4cytoscape
sys
```
- **R**   :
    - install librairie on R console :
```
install.packages("tidyverse")
install.packages("rlist")
install.packages("data.table")
install.packages("dplyr")
install.packages("XML")
install.packages("methods")
install.packages("xlsx")
install.packages("jsonlite")
install.packages("writexl")
install.packages("ontologyIndex")
```

## Description
- **FOLDER :a_script_curation** :
- From brut  phenopacket to filtered no parent with curation phenopacket 
- Output : curation_input folder
  -   **1_filter_file_phenopacket.py** : 
  -   Description : 
  -   Filter phenopacket remove parent all parent stored in output *output_5HPO\curation_input\parent_phenopacket.tsv*
  -   Filter phenopacket keep only phenopacket 5HPO stored in output folder *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete\*
  
  -   **2_curation.py** : 
  -    Description :
  -    Lauch R script *curation_excluding_redundant_terms_with_ontologyIndex.R* to execute **ontologyIndex** on phenopacket output stored in folder :  *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\*
 
  -   **3_transform_phenopackets_cleaned_for_RunSolveRD.py** : 
  -   Description :
  -   Set json structure for RunSolveRD algo
  -   Rename phenopacket file output stored in folder  *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\*
  
  -   **4_RUN_runSolvedRD.py** : 
  -    Bash command run algo **RunSolveRD** output stored in folder *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results*

  -   **curation_excluding_redundant_terms_with_ontologyIndex.R** : 
  -   R file run in script 2_curation.py 
  -   /!\ need to change path variable.

- Phenopackets used for la suite is on folder *Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete* 

- **FOLDER :b_script_get_gene** :
- Description : pas de chiffre devant car pas d'ordre need ici
  -   **get_gene_from_case.py** : 
  -   Get genes from phenopackets output save in folder *output_5HPO\gene_info\gene_from_json_case.tsv*
  -   Get genes from phenopackets solved only GPAP output save in folder *output_5HPO\gene_info\gene_from_solved_case.tsv*
  
  -   **get_gene_from_orpha.py** : 
  -   Get genes from orphanet website based on folder *input\RunSolveRD* *en_product1_sep2022.xml* file and *en_product6_sep2022.xml* file and *input\product_ORPHApackets_childs.xml* file output stored in *output_5HPO\gene_info\pd_1_6_child.tsv*
  
- **FOLDER :c_script_step** :
- Description : pas de chiffre devant car pas d'ordre need ici
  -   **1_stepA1_A2.py** : 
  -   Description
  
  -   **2_stepB1.py** : 
  -    Description
 
  -   **3_stepB2.py** : 
  -   Description
  
  -   **4_stepC1.py** : 
  -    Description

  -   **5_stepC2.R** : 
  -   Description

- **FOLDER :d_script_minitsv** :
- Description : pas de chiffre devant car pas d'ordre need ici
  -   **build_minidf.py** : 
  -   Description

- **FOLDER :e_script_cytoscape** :
- Description : pas de chiffre devant car pas d'ordre need ici
  -   **1_cytoscape_all.py** : 
  -   Description
 
  -   **2_wikipathway_all_py4cy.py** : 
  -   Description
    
  -   **3_get_ALL_nodetype.py** : 
  -   Description
 
 - **FOLDER :e_script_cytoscape** :
  - Description : pas de chiffre devant car pas d'ordre need ici
  -   **buildjson_all.py** : 
  -   Description
 
  - **FOLDER :g_script_homepageJS** :
  - Description : pas de chiffre devant car pas d'ordre need ici
  -   **homepage.py** : 
  -   Description


- Remarque:
    - blabla
    - blabla
    - blabla

section R modifier path on the R script [name]
```
folder = "C:\\Users\\mchahdil\\Documents\\Freeze3\\\\output_5HPO\\\\curation_input\\"
setwd(folder)

input_json = "Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete\\"
output_json = "Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\\"
```

mettre paragraphe runsolvedRD
  expliquer les output generer 

expliquer le déroulement de ce code
paragraphe input 

mettre paragraphe sur le temps de chaque script 

docker paragraphe 







