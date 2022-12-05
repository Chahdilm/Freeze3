### Auteur
Maroua CHAHDIL  


# Description

The service is coded in order to visualise phenotypic similarity network thought 
biological networks. Its offer to users like clinician who belong to
ERNs an interactive and user-friendly experience and can drive to 
diagnostic hypothesis. It's related to Solve rd project.
This is a portal for phenotipic similarity analysis, it will be 
deployed for clinicians and guide to diagnostics hypothesis.

Data are generated with python, visualization is done with JS with Cytoscape JS.

## SolveRD projet
Solve rd project aim to helped through diagnostics hypothesis on unsolved patients with rare disease done by clinicians who belong to ERNs. WP1 contains unsolved and solved-rd cases data in a standardised way, unsolved one are related to the RD ontology RDCO. Orphanet develops a three steps workflow featuring with CNAG. Data visualisation of cases based on the three steps is settle thanks to biological network through Cytoscape js library. This aims to code an interactive js platform. Explorations, analysis and filtration of genes, cases and Orphacodes drive to diagnostic hypothesis. Several options for filtration, analysis, exploration or exportation are available on the platform like edges score filtration, ERN section, build a sub-graph, json exports and so on. Here, we specifically present jamborees cases network based on step C1. Cytoscape js lead to easily explore data and it can be link to other databases like GO, WikiPathways and so on.


# Version : 

- Project done on Window 10 21H2

  - R 4.2.1
  - Python 3.8
  - Java (openjdk-8-jre-headless)


# Package   :
## Python
```
pandas==1.3.4
numpy==1.21.3
openpyxl==3.0.9

# Library system
os 
json
subprocess
datetime
time
multiprocessing
threading
sys
```
## R
```
install.packages("data.table")
install.packages("dplyr")
install.packages("jsonlite")
install.packages("ontologyIndex")
```


# Getting started :
  -  docker pull or unzip orphanet_gpap_workflow.zip 
  -  go to orphanet_gpap_workflow folder
      - two folder on orphanet_gpap_workflow : docker and SolveRD
  - go to the file docker/workflow_docker.sh
      - Change the section in bold of path on lines 3;49;53;54
      - Exemple on line 49 : **/path/to/change/** SolveRD/output_files/log/docker_logs.log
      - close the file 
  - go to the file docker/docker-compose.yml
      - Change the path on lines 8;19;20;33
      - Exemple fo the volumes section on the python_process service : **/path/to/change/** SolveRD/output_files/log/docker_logs.log 
      - close the file
  - go to bash terminal
      - go to the docker folder location/path
      - run the command : **sh workflow_docker.sh**


## **Docker** :

**Docker** is a platform that organize each element (languages like R, python…) of the project. Thanks to docker the project can be shared and implemented easily.


- Docker is used on the workflow thought **Docker compose** 
  - Docker Compose is a tool which defines multi-images for each element. Dockerfile is used in the docker compose file.
  - Docker compose contains threes images from Docker hub (an open source images libraries).
    -   rocker/rstudio:4.2.1 
        - To run the R script which contains OntholofyIndex 
    - openjdk-8-jre-headless
        - To run the runsolveRD  algorithm (*S. Köhler RunSolveRD.jar package*)
    - python:3.8-slim
       - for the others bin files
       -  A **Dockerfile** is used to load pip and python libraries
          - Dockerfile contains all instructions command to describe elements of the project. 



# Folder descriptions :
### SolveRD                               
  - docker-compose.yml                    
    + docker-compose file with all images needed
  - Dockerfile                            
    + Dockerfile for the python part dockerfile is connected to docker-compose file
  - requirements.r  *cf Package section*                         
    + R libraries
  - requirements.txt  *cf Package section*                        
    + Python libraries
  - workflow_docker.sh                    
    + shell script with all command to load scripts 
    
### SolveRD folder :
  - binall
    - contains scripts in python, R and java this will be explain later
    - 
- input_files
  - gpap_variants_files : variants files from GPAP **need upload often**
  - hpo_files : hpo file needed for the runsolveRD algo **need upload often**
  - orphanet_files : products from orphanet **need upload often**
  - ped_files : files with parents/children info for each phenopackets **need upload often**
  - phenopackets_input_files : raw phenopackets from GPAP **need upload often**
  - wikipathways_files : wikipathways file **need upload often**
  - 2022_10_14_cohort_ERN.csv ERN temps **need upload often** (tmp file wait for the ERN info on phenopacket files)
  -
- --> an output_files will be created with all output folder and file once the workflow is over.
- readme file  
- howto file
 

    
# Note:
  -   Phenopackets are json files describing HPO's of one patient.
  -   Cases is a synonym for phenopacket. Case is used to describe patients on the runsolveRD algorithm, the case-case similarity algorithm is run and or the case-ORPHAcode similarity algorithm is run. Cases also describe unsolved or solved patients.


# Description of bin folders and files 
## FOLDER : a_script_curation
### Curation  
- The curation process filters the raw phenopaquets  by excluding parents, removing redundant HPO terms and keeping only phenopaquets containing more than 5 HPOs.
- Python, R and the runsolveRD algorithm in java are used.  The algorithm performs similarity score between 2 entities (orpha/pheno or pheno/pheno). In R, the ontologyIndex package is used to remove HPOs redundant terms. Thought the HPO classification in the aim to keep the lowest HPO annotation.
- Input : Raw phenopackets
- Output : Final study population phenopackets (*output_files/study_population*) and Results of runsolveRD computations (*output_files/study_population/results*)
- 
-   **filter_file_phenopacket.py** : 
    -   Filter phenopackets by removing the parents. 
    -   Filter phenopacket keep only phenopacket 5 HPOs without parents.
    -   Phenopackets files output stored on folder *output_files/json_filterHPO* (the folder will be removed, it's a temp folder)

 
  - **curation_excluding_redundant_terms_with_ontologyIndex.R** : 
      -   Run the R file with OnthologyIndex package to remove redunbant HPO terms.
      -   Phenopackets files output stored on folder *output_files/json_curation_tmp* (the folder will be removed, it's a temp folder)


  -   **transform_phenopackets_cleaned_for_RunSolveRD.py** : 
      -   Define a json structure compatible with the RunSolveRD algorithm
      -   Rename phenopacket files output
      -   Final study population files ready
      -   Phenopackets files output stored on folder *output_files/json_curation_tmp* ( the folder will be removed its a temp folder)
      
  -   **runSolveRdAnalysis.jar** :
      -  Bash command run algo **RunSolveRD**.
      -  Results of runsolveRD computations files ready
      -  Phenopackets files output stored on folder *output_files/study_population*
 
## FOLDER - b_script_get_gene
### Genes selections
- Get all genes from phenopackets and Orphacodes.
- Python is used.  
- Input :  Final study population phenopackets and Orphanet product 6 and 1.
- Output : Excels files with genes Orphacode and gene phenopackets descriptions
-
-   **get_gene_from_case.py** : 
      -   Get genes from phenopackets, file output stored in *output_files/gene_info/gene_from_solved_case.tsv* and *output_files/gene_info/gene_from_solved_json.tsv* (for B1 step only)
      - 
  -   **get_gene_from_orpha.py** : 
-   Get genes from Orphanet website based on  *en_product1_sep2022.xml* file and *en_product6_sep2022.xml* file and *input\product_ORPHApackets_childs.xml*. File output stored in *output/gene_info/pd_1_6_child.tsv*
  -   
  - **3_get_ALL_nodetype.py** : 
    -   For each phenopacket get the SOLVE or UNSOLVED information, file output stored in *output_files/node_type/nodetype_PhenoType.tsv*
    -   For each phenopacket get the ERN information, file output stored in *output_files/node_type/nodetype_ERN.tsv*
    -   For each ORPHAcode get their nomemclature information  (disorder or subtype) excluding group of disorder, file output stored in *output_files/node_type/nodetype_DisorderType.tsv* 
    -   For each genes  get the HGNC id and the Orphanet gene id, file output stored in *output_files/node_type/nodetype_gene_id.tsv* 

## **FOLDER - c_script_step** :
### 3-step programmatic methodology.
- Python is used.  
- Input :  Final study population phenopackets (*output_files/study_population*)
- Output : Excels files for each steps per phenopackets 
-
-   **1_stepA1_A2.py** : 
    -   Check the similarity between a case related to rare diseases with the same gene output stored in folder *output_files\stepA1.tsv*. 
    -   Check the similarity between a case related to rare diseases with different genes output stored in folder *output_files\stepA2.tsv*.
  
-   **2_stepB1.py** : 
    -   Find the corresponding between cases output stored in folder *output_files\stepB1.tsv*.
  
-   **3_stepB2.py** : 
    -   Finding the unsolved cases similar to the solved case in folder *output_files\stepB2.tsv*.
  
-   **4_stepC1.py** : 
    -   Working out unsolved cases similar to the RD related to the solved case *output_files\stepC1.tsv*.

-   **5_stepC2.py** : 
    -   Working out  cases similar to the case from the similar RD related to the solved case *output_files\stepC2.tsv*.
    
  -   **build_minidf.py** :
    - Create a dataframe for each phenopacket on each step.
    - input all tsv steps from the step folder stored in *output_files\*. 
    - output files stored in *output_files\tsv*.



## **FOLDER - d_script_cytoscape_json** :
  -   **cytoscape_all.py** : 
      -   Input  are output from *build_minidf.py*.
      -   Convert a dataframe into a dataframe compatible with Cytoscape.
      -   Then add Wikipathways database genes/pathways interactions for each genes in the network.
      -   output stored on the folder *output_files\cytoscape_wikipathways*

  -   **buildjson_all.py** : 
      -   Input are output of **1_cytoscape_all.py** all dataframe.
      -   Build a json file for each dataframe output from **cytoscape_all.py** to use them on the application done with cytoscape JS.
      -   For each element (cases,genes,Orphacodes,pathways) a node is build.
      -   For each interaction Orphacode/case or case/case or orphacode/gene or case/gene or orphacode/orphacode an edge is build.
      -   Output are stored on the folder *output_files\json*.
      
  -   **homepage.py** : 
      - Build the homepage of the application on javascript (done with the library Cytoscape JS).
      - Check the similarity between a case related to another case for all phenopackets.
      - Build the dataframe then the Cytoscape dataframe and the json file for Cytoscape js.









