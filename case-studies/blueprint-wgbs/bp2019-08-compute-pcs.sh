#$ -pe local 5
#$ -cwd
# -l mem_free=4G,h_vmem=4G

Rscript /users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs/bp2019-08-compute-pcs.R
