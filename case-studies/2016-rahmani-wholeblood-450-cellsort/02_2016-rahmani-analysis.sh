#$ -cwd
#$ -l mem_free=20G,h_vmem=20G

Rscript /users/shicks1/projects/methylCCPaper/case-studies/2016-rahmani-wholeblood-450-cellsort/01_2016-rahmani-create-data-object.R
Rscript /users/shicks1/projects/methylCCPaper/case-studies/2016-rahmani-wholeblood-450-cellsort/02_2016-rahmani-analysis.R
