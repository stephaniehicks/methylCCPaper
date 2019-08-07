#$ -pe local 5
#$ -cwd
# -l mem_free=7G,h_vmem=7G

Rscript /users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs/bp2019-05-apply-qc.R
