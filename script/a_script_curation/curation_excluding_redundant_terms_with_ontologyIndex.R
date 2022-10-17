rm(list=ls())

# Libraries


library(tidyverse)
library(rlist)
library(data.table)
library(dplyr)
library(XML)
library(methods)
library(xlsx)
library(jsonlite)
library(writexl)
library(ontologyIndex)
data(hpo)
hpo_table = as.data.frame(hpo)


folder = "C:\\Users\\mchahdil\\Documents\\Freeze3\\\\output_5HPO\\\\curation_input\\"
setwd(folder)

input_json = "Freezes1_2_3_noduplicates_noparents_with_5phenotypes_complete\\"
output_json = "Freezes1_2_3_noduplicates_noparents_with_5phenotypes_aftercuration_complete\\"



setwd(paste(folder, input_json ,sep="/"))

# # Open phenotype datas json
fichiers = list.files(path = paste(folder,input_json,sep="/"), 
                      pattern = ".json", all.files = FALSE, 
                      full.names = FALSE, recursive = TRUE,
                      include.dirs = FALSE, no.. = FALSE)
file=list()
indice_list=1
n=length(fichiers)
names=matrix()
for (i in 1:n)
{
  #Nom du fichier/dossier
  names[[indice_list]]=fichiers[i]
  #Fichier
  file[[indice_list]]=fromJSON(readLines(con=fichiers[i]), flatten = T)
  indice_list=indice_list+1
}


# We remove the HPO negated
file_hpo_list = list()
indice_list = 1
for (i in 1:length(file))
  if (length(file[[i]]$phenopacket$phenotypicFeatures > 0))
  {
    file_hpo_list[[indice_list]] = cbind.data.frame(Pid = file[[i]]$phenopacket$id,
                                                    HPO = file[[i]]$phenopacket$phenotypicFeatures)
    indice_list = indice_list + 1
  }
file_hpo = rbindlist(file_hpo_list, fill = T)


file_hpo = file_hpo[is.na(file_hpo$HPO.negated), ]
file_hpo = select(file_hpo, -c(HPO.negated))

# We remove the wrong annotations
file_hpo = file_hpo[!file_hpo$HPO.type.id == file_hpo$HPO.type.label, ]
file_hpo = filter(file_hpo, !is.na(file_hpo$HPO.type.label))
file_hpo = filter(file_hpo, HPO.type.label != 'Invalid id')


### Remove redundant/implied terms from a set of terms

# We have a list with all Pid with its HPO
file_hpo_by_pid = split(file_hpo, by = 'Pid')

# Function Remove redundant/implied terms from a set of terms
indice_list = 1
a = list()
for (i in 1:length(file_hpo_by_pid))
  for (j in 1:length(file_hpo_by_pid[[i]]$HPO.type.id))
  {
    a[[i]] = cbind.data.frame(Pid = file_hpo_by_pid[[i]]$Pid[[1]],
                              HPO.type = minimal_set(hpo, file_hpo_by_pid[[i]]$HPO.type.id))
    
    indice_list = indice_list + 1
  }

a = rbindlist(a)


phenopacket_hpo = select(file_hpo, c(HPO.type.id, HPO.type.label))
phenopacket_hpo = filter(phenopacket_hpo, !duplicated(phenopacket_hpo))


# We merge the table without redundant terms with the labels of ids
b = merge(a, phenopacket_hpo, by.x = 'HPO.type', by.y = 'HPO.type.id', all.x = T, allow.cartesian = T)
b = filter(b, !duplicated(b[, 1:2]))
b = select(b, c(Pid = Pid, HPO.type, HPO.label = HPO.type.label))

# We remove all the labels that are not in the HPO database
b = b %>%
  arrange(Pid)


# We put this table in a list for each Pid
c = split(b, by = 'Pid')


# We merge the gene for the output
`%notin%` <- Negate(`%in%`)
file_output = list()
for (i in 1:length(file))
  if ("genes" %in% names(file[[i]]$phenopacket)) {
    if ("resolutionStatus" %in% names(file[[i]]))
      {
       file_output[[i]] = list(c[[i]], 
                               gene = file[[i]]$phenopacket$genes,
                               status = file[[i]]$resolutionStatus,
                               disease = file[[i]]$phenopacket$diseases)
       
      } else if ("resolutionStatus" %in% names(file[[i]]$interpretation)) {
    file_output[[i]] = list(c[[i]], 
                            gene = file[[i]]$phenopacket$genes,
                            status = file[[i]]$interpretation$resolutionStatus,
                            disease = file[[i]]$phenopacket$diseases)
      }

  } else if ("genes" %notin% names(file[[i]])) {
    if ("resolutionStatus" %in% names(file[[i]]))
    {
      file_output[[i]] = list(cbind.data.frame(c[[i]]),
                             genes = list(),
                             status = file[[i]]$resolutionStatus,
                             disease = file[[i]]$phenopacket$diseases)
      
      } else if ("resolutionStatus" %in% names(file[[i]]$interpretation)) {
        file_output[[i]] = list(c[[i]], 
                                gene = file[[i]]$phenopacket$genes,
                                status = file[[i]]$interpretation$resolutionStatus,
                                disease = file[[i]]$phenopacket$diseases)
      }
    
  }


# We make a good format of phenopacket
phenopacket_output = list()
for (i in 1:length(file_output))
{
  phenopacket_output[[i]] = list(id = file_output[[i]][[1]]$Pid[[1]],
                                 phenotypes = cbind.data.frame(type = list(cbind.data.frame(id = file_output[[i]][[1]]$HPO.type,
                                                                                            label = file_output[[i]][[1]]$HPO.label))),
                                 genes = file_output[[i]]$gene,
                                 status = file_output[[i]]$status,
                                 disease = file_output[[i]]$disease)
}


# Multiple export json
for (i in 1:length(phenopacket_output))
{
  write_json(phenopacket_output[[i]],
             path = paste0(folder, "/",output_json, "/", fichiers[i]),
             auto_unbox = T, pretty = T, flatten = T)
}




### Remove cases with less 5 HPO after ontologyIndex curation

setwd(paste0(folder, "/", output_json))


### Open phenotypes
fichiers = list.files(path = paste0(folder, "/", output_json),
                      all.files = FALSE, full.names = FALSE, recursive = TRUE,
                      include.dirs = FALSE, no.. = FALSE, pattern = ".json")

# On recupere tous les fichiers
file=list()
indice_list=1
for (i in 1:length(fichiers)) 
{
  file[[indice_list]]=fromJSON(readLines(con=fichiers[i]), flatten = T)
  indice_list=indice_list+1
}



### File with HPO
file_HPO = list()
indice_list = 1
for (i in 1:length(file))
  if (length(file[[i]]$phenotypes) > 0)
  {
    file_HPO[[indice_list]] = cbind.data.frame(Pid = fichiers[[i]], 
                                               file[[i]]$phenotypes)
    indice_list = indice_list + 1
  }

file_HPO = rbindlist(file_HPO, fill = T)


# We count the number of HPO
HPO_count = data.frame(table(file_HPO$Pid))
names(HPO_count) = c('Pid', 'HPO')


# We filter with 5 HPO or more
file_5HPO = filter(HPO_count, HPO > 4)


# All files in the folder
all_file = list()
for (i in 1:length(fichiers))
{
  all_file[[i]] = cbind.data.frame(Pid = fichiers[[i]])
}

all_file = rbindlist(all_file)


### All Pid that we should remove
to_remove = anti_join(all_file, file_5HPO, by = 'Pid')


# And delete them
for (i in 1:nrow(to_remove))
{
  unlink(to_remove[i])
}

