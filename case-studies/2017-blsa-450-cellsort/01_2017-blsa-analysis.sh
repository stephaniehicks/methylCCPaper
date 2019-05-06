#$ -cwd
#$ -l mem_free=50G,h_vmem=50G

Rscript /users/shicks1/projects/methylCC/dataBLSA_knownProps450k/dataBLSA.R 

Rscript /users/shicks1/projects/methylCCPaper/case-studies/2017-blsa-450-cellsort/01_2013-liu-create-data-object.R
Rscript /users/shicks1/projects/methylCCPaper/case-studies/2017-blsa-450-cellsort/02_2013-liu-analysis.R
