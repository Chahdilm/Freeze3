

import subprocess
from script.path_variable import *
# C:\Program Files\R\R-4.2.1
# C:\Program Files\R\R-4.2.1>\Rscript Rcode.r
print("2_curation.py\tRun R script ontologyIndex ")

# subprocess.call([r'C:\Program Files\R\R-4.2.1\bin\Rscript', "/Users/mchahdil/Documents/Freeze3/script/a_script_curation/curation_excluding_redundant_terms_with_ontologyIndex.R"])

command = (PATH_SCRIPT + str(r"a_script_curation/curation_excluding_redundant_terms_with_ontologyIndex.R"))

subprocess.call([PATH_R, str(command)])
