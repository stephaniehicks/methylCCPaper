library(minfi)
library(quadprog)
library(methylCC)

################################################
### Target Data: 806 450K DNA methylation samples with known cell type proportions
################################################

dataPath_blsa <- "/users/shicks1/data/DNAm/2017-blsa-450-cellsort"
workingDir_blsa <- "/users/shicks1/projects/methylCCPaper/case-studies/2017-blsa-450-cellsort"


# Flow Sorted measured cell type proportions (806 samples)
trueProps <- read.csv(file.path(dataPath_blsa, 
                                "2017-blsa-celltype-counts.csv"), 
                      header = TRUE)
colnames(trueProps) <- c("sampleID", "baso", "eos", "neu", "lymph", "mono")

# Create a "Gran" column (Baso + Eos + Neu)
trueProps$gran <- trueProps$baso + trueProps$eos + trueProps$neu
trueProps[,-1] <- trueProps[,-1] / 100

# Load 800 samples run on 450K DNA methylation array 
#       with corresponding cell type proportions
load(file.path(dataPath_blsa, "2017-blsa-rgset.rda")) 
RGset <- updateObject(RGset) 

# need to fix "non-forming" names
pData(RGset)$Sex <- ifelse(pData(RGset)$Sex == 1, "M", "F")
pData(RGset)$Slide <- as.numeric(pData(RGset)$Slide)

# Need to filter: 806 samples in trueProps, but only 800 samples in RGset
trueProps <- trueProps[match(pData(RGset)$Sample.ID, trueProps$sampleID), ]
colnames(trueProps) <- c("samples", "Baso", "Eos", "Neu", "Lymph", "Mono", "Gran")

# run minfi::estimateCellCounts()
# counts450K <- estimateCellCounts(RGset)
# save(counts450K, file = file.path(workingDir_blsa, "ests/dataBLSA-minfi-estimateCellCounts.RData"))

# compare to compote::estimatecc method
output <- find_dmrs(verbose = TRUE, gr_target = NULL, 
                    num_regions = 50, num_cpgs = 50,
                    include_cpgs = FALSE, include_dmrs = TRUE, 
                    bumphunter_beta_cutoff = 0.2, # defined by bumphunter
                    dmr_up_cutoff = 0.50, # smaller is better
                    dmr_down_cutoff = 0.40, # smaller is better
                    dmr_pval_cutoff = 1e-11, # default of 1e-11
                    cpg_pval_cutoff = 1e-08, # default of 1e-08
                    cpg_up_dm_cutoff = 0, # ranges from -infinity to 0
                    cpg_down_dm_cutoff = 0, # ranges from 0 to infinity
                    pairwise_comparison = FALSE)

set.seed(12345)
est <- estimatecc(object = RGset,
                  find_dmrs_object = output, 
                  verbose = TRUE, epsilon = 0.01, 
                  max_iter = 100, init_param_method = "random")
save(est, file = file.path(workingDir_blsa, "ests/dataBLSA-methylCC-estimatecc.RData"))
rm(est) 


