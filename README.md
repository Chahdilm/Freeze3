### Auteur
Maroua CHAHDIL  


# Description

The service is coded in order to viasualise phenotypic similatiry network throught biological network. Its offer to users like clinician who belong to ERNs an interactive and user-friendly experience and can drive to diagnostic hypothesis. It's related to Solve rd project.
This is a portal for phenotipic similarity analysis, it will be deploy for clinicians and guide to diagnostics hypothesis.

Data are generated with python, visualization is done with JS with cytoscape JS. 
The data file formats are: tsv, json, xlsx, txt and so on.


## SolveRD projet
Solve rd project aim to helped through diagnostics hypothesis on unsolved patients with rare disease done by clinicians who belong to ERNs. WP1 contains unsolved and solved-rd cases data in a standardised way, unsolved one are related to the RD ontology RDCO. Orphanet develop a three steps workflow featuring with CNAG. Data visualisation of cases based on the three steps is settle thanks to biological network through cytoscape js library. This aims to code an interactive js platform. Explorations, analysis and filtration of genes, cases and Orphacodes drive to diagnostic hypothesis. Several options for filtration, analysis, exploration or exportation are available on the platform like edges score filtration, ERN section, build a sub-graph, json export and so on Here, we specifically present jamborees cases network based on step C1. Cytoscape js lead to easily explore data and it can be link to other databases like GO, WikiPathways and so on.  


# Version : 

- Project done on Window 10 21H2

  - R 4.2.1
      - R Language for IntelliJ on pycharm (223.6646—223.6646.*)

  - Python 3.9


# Package needed  :
## Python
  Install libraries on the python console : use conda or pip or pip3
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
   Install libraries on R console :
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
  -   R script *curation_excluding_redundant_terms_with_ontologyIndex.R*.
      -   Change the path of the variable : **folder**.
```
folder = "C:\\Users\\username\\...\\Freeze3\\output_5HPO\\curation_input\\"
setwd(folder)

input_json = "Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete\\"
output_json = "Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\\"
```

  -   Python script *path_variable.py*.
      -   Change the path of the variable **PATH**. 
            -   PATH describes the location of the main file named Freeze3. All variables that store paths depend on the PATH variable.
            -   Exemple :
                  ```
                  PATH= r"C:\Users\username\Documents\Freeze3\\"
                  ```
      -   Change the path of the variable **PATH_R**. 
            -   PATH_R indicates where R is stored on the computer. It is necessary to modify it depending on the location of R and the OS use (window,linux or mac).
            -   Exemple : 
                  ```
                  PATH_R= r"C:\Program Files\R\R-4.2.1\bin\Rscript"
                  ```
                  
  - Python script *main.py*
      -   If needed uncomment the list.remove('__init__.py') part for all list. Sometime __init__.py is included on the list thus it should be removed, if it's not into the list to need to change the code.
           -   Exemple : 
                  ```
                  # files_all_step.remove('__init__.py') # uncomment if needed
                  ```  

 
 
# Getting started :
  -   Extract the project on docker [command ]  **WIP**
  -   Download the zip folder.
  -   The folder is called Freeze3 and contains the folders input (all the input files needed for the process), protocol curation (documentation) and script (containing the pyhton code).
  -   Go to the file named *a_script_curation/curation_excluding_redundant_terms_with_ontologyIndex.R* and change the path.
  -   Go to the file named *path_variable.py* and change the variable **PATH** put the actual path where the folder Freeze3 is.
  -   Go to the script folder and run the file main.py. 
  -   *main.py* execute all scripts on the right order. 
      - a) Filter and extract the phenopackets  of interest, run the curation on R with **OntologyIndex** and launch the RunsolveRD (phenopackets similarities score algorithm).
      - b) Set all genes. 
      - c) Lauch all steps.
      - d) Build for each phenopacket of each steps its sub-dataframe.
      - e) Generate the dataframes of each phenopacket of each step compatible with Cytoscape. Then add informations from **Wikipathway**.
      - f) From the dataframes, Cytoscape builds json files compatible with Cytoscape JS for the service. 
      - g) Folder that defines the homepage of the service from scratch.


# Note:
  -   Phenopackets are json files that describe HPO's of one patient.
  -   Cases is a synonym for phenopacket. Case is used to describe patients on the runsolveRD algorithm, the case-case similarity algorithm is run and or the case-ORPHAcode similarity algorithm is run. Cases also describe unsolved or solved patients.


# Description
## **FOLDER - a_script_curation** 
File are numbering according to the execution order:
From raw phenopacket to filtered phenopacket without parent and with curation and min 5HPOs. 
Output : *output_5HPO\curation_input* folder.
  -   **1_filter_file_phenopacket.py** : 
      -   Filter phenopackets by removing the parents. All parents are stored in the output file *output_5HPO\curation_input\parent_phenopacket.tsv*.
      -   Filter phenopacket keep only phenopacket 5HPO without parents, stored in output folder *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete\*.

  -   **2_curation.py** : 
      -    Run R script *curation_excluding_redundant_terms_with_ontologyIndex.R* to execute the **ontologyIndex** R library to make curations on phenopackets outputs stored in folder :  *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\*
 
  -   **3_transform_phenopackets_cleaned_for_RunSolveRD.py** : 
      -   Define a json structure compatible with the RunSolveRD algorithm
      -   Rename phenopacket file output stored in folder  *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\*

  -   **4_RUN_runSolvedRD.py** : 
      -    Bash command run algo **RunSolveRD** output stored in the folder *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results*

  -   **curation_excluding_redundant_terms_with_ontologyIndex.R** : 
      -   Run the R file trhought python script : *2_curation.py*. 
      -   /!\ need to change the folder path variable.

#### Phenopackets used for next steps are in the folder : *Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete* 

## **FOLDER - b_script_get_gene** :
No execution order-   **get_gene_from_case.py** : 
      -   Get genes from phenopackets, file is save on file *output_5HPO\gene_info\gene_from_json_case.tsv*
      -   Get genes from phenopackets solved only GPAP output save in folder *output_5HPO\gene_info\gene_from_solved_case.tsv*
 
  -   **get_gene_from_orpha.py** : 
      -   Get genes from orphanet website based on folder *input\RunSolveRD* *en_product1_sep2022.xml* file and *en_product6_sep2022.xml* file and *input\product_ORPHApackets_childs.xml* file output stored in *output_5HPO\gene_info\pd_1_6_child.tsv*
  
## **FOLDER - c_script_step** :
File are numbering according to the execution order:
  -   **1_stepA1_A2.py** : 
      -   Input : phenopackets results from RunSolveRD algo folder : *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results\results_noduplicates\resultsORDO*.
      -   Check the similarity between a case related to rare diseases with the same gene output stored in folder *output_5HPO\stepA1.tsv*. 
      -   Check the similarity between a case related to rare diseases with different genes output stored in folder *output_5HPO\stepA2.tsv*.
  
  -   **2_stepB1.py** : 
      -   Input : phenopackets result from RunSolveRD algo folder : *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results\results_noduplicates\resultsPhenopackets* .
      -   Find the corresponding between cases output stored in folder *output_5HPO\stepB1.tsv*.
  
  -   **3_stepB2.py** : 
      -   Input : phenopackets result from RunSolveRD algo folder : *output_5HPO\curation_input\Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\results\results_noduplicates\resultsORDO* and output B1 cases.
      -   Finding the unsolved cases similar to the solved case in folder *output_5HPO\stepB2.tsv*.
  
  -   **4_stepC1.py** : 
      -   Input : phenopackets result from RunSolveRD algo folder : 
      -   Working out unsolved cases similar to the RD related to the solved case *output_5HPO\stepC1.tsv*.
  

  -   **5_stepC2.R** : 
      -   input phenopacket result from RunSolveRD algo folder : 
      -   Working out  cases similar to the case from the similar RD related to the solved case *output_5HPO\stepC2.tsv*.
  

## **FOLDER - d_script_minitsv** :
No execution order
Create sub_tsv for each phenopacket on each steps.
  -   **build_minidf.py** : 
      -   input all tsv steps from the step folder stored in *output_5HPO\*.

## **FOLDER - e_script_cytoscape** :
File are numbering according to the execution order:
  -   **1_cytoscape_all.py** : 
      -   Input  are output from *build_minidf.py*.
      -   Convert a dataframe into a dataframe compatible with Cytoscape output store on the folder *output_5HPO\cytoscape*
 
  -   **2_wikipathway_all_py4cy.py** : 
      -   Input are Cytoscape dataframes 
      -   Use a tools CyTargetLinker from Wikipathway database information on the network on cytoscape. Output stored on the folder *output_5HPO\cytoscape_wikipathways*.
      -   It s done on python thanks to py4cytoscape librairy. 
    
  -   **3_get_ALL_nodetype.py** : 
      -   For each phenopacket informations of their statue SOLVE or UNSOLVED
      -   For each phenopacket informations of their their ERN origin 
      -   For each ORPHAcode information of their nomemclature (disorder or subtype)
      -   Input are product6 (*input\RunSolveRD\en_product6_sep2022.xml*), phenopacket description to get the statue and the ERN file (*input\2022_10_14_cohort_ERN.csv*)

## **FOLDER - f_script_json** :
No execution order
Build a json file for each dataframe output from **2_wikipathway_all_py4cy.py** to use them on the application done with cytoscape JS.
  -   **buildjson_all.py** : 
      -   Input are output of **1_cytoscape_all.py** all dataframe.
      -   For each elements(cases,genes,orphacodes,pathways) a node is build.
      -   For each interactions orphacode/case or case/case or orphacode/gene or case/gene or orphacode/orphacode an edge is build.
      -   Output are stored on the folder *output_5HPO\json*.
 
 
## **FOLDER - g_script_homepageJS** :
Build the homepage of the application on javascript (done with the librairy cytoscape JS).
  -   **homepage.py** : 
      -   Check the similarity between a case related to another case for all phenopackets..
      -   Build the dataframe then the cytoscape dataframe and the json file for cytoscape js.


## **runsolvedRD** :
Phenotypic similarity is done using Resnik symmetric alorithm on RunsolveRD.Similarity depend on phenotypic annotations thanks to Human phenotype onthology (HPO). Similarity betweend cases and/or betwend cases and ORPHAcodes is done.The top 50 first ranked results are selected and store on the json resuls.


## **Docker** :
[...] **WIP**







