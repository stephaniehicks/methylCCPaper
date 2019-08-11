#$ -pe local 10
#$ -cwd
# -l mem_free=5G,h_vmem=6G

Rscript /users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs/bp2019-07-findDMRs-cd4t.R
