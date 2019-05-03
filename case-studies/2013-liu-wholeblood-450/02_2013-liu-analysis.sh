#$ -cwd
#$ -l mem_free=20G,h_vmem=20G

Rscript /users/shicks1/projects/methylCCPaper/case-studies/2013-liu-wholeblood-450/01_2013-liu-create-data-object.R
Rscript /users/shicks1/projects/methylCCPaper/case-studies/2013-liu-wholeblood-450/02_2013-liu-analysis.R
