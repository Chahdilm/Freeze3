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
## Remarque:
expliquer phenopacket et cases 

## Description
#### **FOLDER :a_script_curation** :
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
- Description : exe depend on number order
  -   **1_stepA1_A2.py** : 
  -   input phenopacket result from RunSolveRD algo folder : *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results\results_noduplicates\resultsORDO*
  -   Check the similarity between a case related to rare diseases with the same gene output stored in folder *output_5HPO\stepA1.tsv* 
  -   Check the similarity between a case related to rare diseases with different genes output stored in folder *output_5HPO\stepA2.tsv*
  
  -   **2_stepB1.py** : 
  -   input phenopacket result from RunSolveRD algo folder : *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results\results_noduplicates\resultsPhenopackets* 
  -   Find the corresponding between cases output stored in folder *output_5HPO\stepB1.tsv*
  
  -   **3_stepB2.py** : 
  -   input phenopacket result from RunSolveRD algo folder : *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results\results_noduplicates\resultsORDO* and output B1 cases
  -   Finding the unsolved cases similar to the solved case in folder *output_5HPO\stepB2.tsv*
  
  
  -   **4_stepC1.py** : 
  -   input phenopacket result from RunSolveRD algo folder : [...] en cours
  -   Working out the unsolved cases similar to the RD related to the solved case *output_5HPO\stepC1.tsv*
  

  -   **5_stepC2.R** : 
  -   input phenopacket result from RunSolveRD algo folder : [...] en cours
  -   Working out the unsolved cases similar to the case from the similar RD related to the solved case *output_5HPO\stepC2.tsv*
  

- **FOLDER :d_script_minitsv** :
- Create sub_tsv for each phenopacket on each steps
  -   **build_minidf.py** : 
  -   input all tsv steps from the step folder stored in *output_5HPO\* 

- **FOLDER :e_script_cytoscape** :
- Description : exe depend on number order
  -   **1_cytoscape_all.py** : 
  -   input are output from *build_minidf.py*
  -   Convert a dataframe into a dataframe compatible with Cytoscape output store on the folder *output_5HPO\cytoscape*
 
  -   **2_wikipathway_all_py4cy.py** : 
  -   input are cytoscape dataframes 
  -   use a tools CyTargetLinker from Wikipathway database to create add cytoscape information on the network,output store on the folder *output_5HPO\cytoscape_wikipathways*
  -   It s done on python thanks to py4cytoscape librairy. 
    
  -   **3_get_ALL_nodetype.py** : 
  -   For each phenopacket information of their statue SOLVE or UNSOLVED
  -   For each phenopacket information of their their ERN origin 
  -   For each ORPHAcode information of their nomemclature (disorder or subtype)
  -   input are product6 (*input\RunSolveRD\en_product6_sep2022.xml*), phenopacket description to get the statue and the ERN file (*input\2022_10_14_cohort_ERN.csv*)

- **FOLDER :f_script_json** :
Build a json file for each dataframe output from **2_wikipathway_all_py4cy.py** to use them on the application done with cytoscape JS.
-   **buildjson_all.py** : 
  -   Input are output of **2_wikipathway_all_py4cy.py** all dataframe.
  -   For each elements(cases,genes,orphacodes,pathways) a node is build.
  -   For each interactions orphacode/case or case/case or orphacode/gene or case/gene or orphacode/orphacode an edge is build
  -   Output are stored on the folder *output_5HPO\json*
 
 
- **FOLDER :g_script_homepageJS** :
Build the homepage for the application on javascript (done with the librairy cytoscape JS)
  -   **homepage.py** : 
  -   Check the similarity between a case related to another case for all phenopackets.
  -   Build the dataframe then the cytoscape dataframe and the json file for cytoscape js.


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







