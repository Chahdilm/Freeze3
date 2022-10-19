### Auteur
Maroua CHAHDIL  


# Description

The service is coded in order to viasualise phenotypic similatiry network throught biological network.Its offer to users like clinician who belong to ERNs an interactive and user-friendly experience and can drive to diagnostic hypothesis. It's related to Solve rd project.
This is a portal for phenotipic similarity analysis,it will be deploy for clinicians  and guide to diagnostics hypothesis.

The data are generate thanks to python, the visualisation is done thanks to JS throught cytoscape JS.
Data file formats are : tsv,json,xlsx,txt and so on.


## SolveRD projet
Solve rd project aim to helped through diagnostics hypothesis on unsolved patients with rare disease done by clinicians who belong to ERNs. WP1 contains unsolved and solved-rd cases data in a standardised way, unsolved one are related to the RD ontology RDCO. Orphanet develop a three steps workflow featuring with CNAG. Data visualisation of cases based on the three steps is settle thanks to biological network through cytoscape js library. This aims to code an interactive js platform. Explorations, analysis and filtration of genes, cases and Orphacodes drive to diagnostic hypothesis. Several options for filtration, analysis, exploration or exportation are available on the platform like edges score filtration, ERN section, build a sub-graph, json export and so on Here, we specifically present jamborees cases network based on step C1. Cytoscape js lead to easily explore data and it can be link to other databases like GO, WikiPathways and so on.  




# Version : 

- Project made on Window 10 21H2
- Compatible with Mac OS Catalina Version 10.15.7 (19H15)  <- a voir encore 

  - R 4.2.1
      - R Language for IntelliJ on pycharm (223.6646—223.6646.*)

  - Python 3.9

# Package needed  :
## Python
  Install librairie on python console : use conda or pip or pip3
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
## R
   Install librairie on R console :
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
# /!\ Code to change  /!\
  -   R script *curation_excluding_redundant_terms_with_ontologyIndex.R*
      -   Change the variable path : **folder**
```
folder = "C:\\Users\\mchahdil\\Documents\\Freeze3\\\\output_5HPO\\\\curation_input\\"
setwd(folder)

input_json = "Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete\\"
output_json = "Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\\"
```
  -   Python script *path_variable.py*
      -   Change the path of the variable **PATH** Exmemple :
```
PATH= r"C:\Users\mchahdil\Documents\Freeze3\\"
```
            -   PATH describe the localisation of the main file named Freeze3. All  variable which store path depend on  the PATH.
      -   Change the path of the variable **PATH_R** Exemple : PATH_R = r'C:\Program Files\R\R-4.2.1\bin\Rscript' Exemple :
            -   PATH_R indicate where R is stored on the computer. Need to change it depending on R localisation and OS.
```
PATH= r"C:\Users\mchahdil\Documents\Freeze3\\"
```
 
# Prise en main :
  -   Extract the project on docker [ mettre la commande ]
  -   Download the zip folder 
  -   The folder is named Freeze3 it containt input (all input files needed for the process) ,protocole curation (documentation) and script (contain the pyhton code) folder
  -   Go to the file named *a_script_curation/curation_excluding_redundant_terms_with_ontologyIndex.R* and change the path on the variable folder 
  -   Go to the file named *path_variable.py* and change the variable PATH put the actual path where the folder Freeze3 is
  -   Go to script folder and lauche the main.py file 
  -   main.py execute all script on the right order 
      - a) Filter and extract phenopackets of interest, execute on R the curation thanks to ontologyIndex and lauch the RunsolveRD  (phenotipic similarity score algorithm)
      - b) Set all genes 
      - c) lauch all steps
      - d) build for each phenopacket of each steps its sub-dataframe
      - e) Generate dataframes of each phenopacket of each steps compatible with cytoscape. Then add wikipathway information
      - f) From dataframe cytoscape build a json file compatible with cytoscape JS for the service 
      - g) Folder which set the home page of the service from scratch.


# Remarque:

  -   Phenopacket are json files which describe hpo for one patient
  -   Cases is a phenopacket synonyme. Case is use to describe patients on the algorithme runsolveRD,went case-case similarity algorithme is lauch. Case describes also           unsolved or solved patient.



# Description
## **FOLDER :a_script_curation** :
From brut phenopacket to filtered one without parent and with curation phenopacket 
Output : *output_5HPO\curation_input* folder
  -   **1_filter_file_phenopacket.py** : 
      -   Filter phenopacket remove parent. All parent are stored in output file *output_5HPO\curation_input\parent_phenopacket.tsv*
      -   Filter phenopacket keep only phenopacket 5HPO without parents, stored in output folder *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete\*

  -   **2_curation.py** : 
      -    Lauch R script *curation_excluding_redundant_terms_with_ontologyIndex.R* to execute the **ontologyIndex** R librairy to make curations on phenopackets outputs stored in folder :  *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\*
 
  -   **3_transform_phenopackets_cleaned_for_RunSolveRD.py** : 
      -   Set json structure compatible with the RunSolveRD algo
      -   Rename phenopacket file output stored in folder  *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\*

  -   **4_RUN_runSolvedRD.py** : 
      -    Bash command run algo **RunSolveRD** output stored in folder *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results*

  -   **curation_excluding_redundant_terms_with_ontologyIndex.R** : 
      -   R file run in script 2_curation.py 
      -   /!\ need to change path variable.

- Phenopackets used for la suite is on folder *Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete* 

## **FOLDER :b_script_get_gene** :
pas de chiffre devant car pas d'ordre need ici
  -   **get_gene_from_case.py** : 
      -   Get genes from phenopackets output save in folder *output_5HPO\gene_info\gene_from_json_case.tsv*
      -   Get genes from phenopackets solved only GPAP output save in folder *output_5HPO\gene_info\gene_from_solved_case.tsv*
  
  -   **get_gene_from_orpha.py** : 
      -   Get genes from orphanet website based on folder *input\RunSolveRD* *en_product1_sep2022.xml* file and *en_product6_sep2022.xml* file and *input\product_ORPHApackets_childs.xml* file output stored in *output_5HPO\gene_info\pd_1_6_child.tsv*
  
## **FOLDER :c_script_step** :
Exe depend on number order
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
  

## **FOLDER :d_script_minitsv** :
Create sub_tsv for each phenopacket on each steps
  -   **build_minidf.py** : 
      -   input all tsv steps from the step folder stored in *output_5HPO\* 

## **FOLDER :e_script_cytoscape** :
Exe depend on number order
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

## **FOLDER :f_script_json** :
Build a json file for each dataframe output from **2_wikipathway_all_py4cy.py** to use them on the application done with cytoscape JS.
  -   **buildjson_all.py** : 
      -   Input are output of **2_wikipathway_all_py4cy.py** all dataframe.
      -   For each elements(cases,genes,orphacodes,pathways) a node is build.
      -   For each interactions orphacode/case or case/case or orphacode/gene or case/gene or orphacode/orphacode an edge is build
      -   Output are stored on the folder *output_5HPO\json*
 
 
## **FOLDER :g_script_homepageJS** :
Build the homepage for the application on javascript (done with the librairy cytoscape JS)
  -   **homepage.py** : 
      -   Check the similarity between a case related to another case for all phenopackets.
      -   Build the dataframe then the cytoscape dataframe and the json file for cytoscape js.


## **runsolvedRD** :
Phenotypic similarity is done using Resnik symmetric alorithm on RunsolveRD.Similarity depend on phenotypic annotations thanks to Human phenotype onthology (HPO). Similarity betweend cases and/or betwend cases and ORPHAcodes is done.The top 50 first ranked results are selected and store on the json resuls.


## **Docker** :
[...] **WIP**







