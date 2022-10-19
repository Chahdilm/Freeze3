import subprocess
from script.path_variable import *

print("4_RUN_runSolvedRD.py\tRun algo ")

command = ('java -Xmx8G -jar ' + PATH_INPUT+str("RunSolveRD\\runSolveRdAnalysis.jar") + " --remove-duplicates --create-report -p " +
          PATH_OUPUT_5HPO_NOPARENT_AFTER_CURATION + " -o " +
          PATH_INPUT+str("RunSolveRD\\hp.obo")    + " -x " +
          PATH_INPUT+str("RunSolveRD\\product4_RunSolveRD.xml") +
          " --numberResultsDiseases 50 --numberResultsPhenopackets 50"
            )


subprocess.run(command, shell=True)


