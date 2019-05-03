#$ -cwd
#$ -l mem_free=20G,h_vmem=20G

Rscript /users/shicks1/projects/methylCCPaper/case-studies/2017-pinoyanes-wholeblood-450-cellsort/01_2017-pinoyanes-create-data-object.R
Rscript /users/shicks1/projects/methylCCPaper/case-studies/2017-pinoyanes-wholeblood-450-cellsort/02_2017-pinoyanes-analysis.R
