#$ -cwd
#$ -l mem_free=20G,h_vmem=20G

Rscript /users/shicks1/projects/methylCCPaper/case-studies/2016-rahmani-wholeblood-450-cellsort/rahmani2016-01-create-data-object.R
Rscript /users/shicks1/projects/methylCCPaper/case-studies/2016-rahmani-wholeblood-450-cellsort/rahmani2016-02-analysis.R
