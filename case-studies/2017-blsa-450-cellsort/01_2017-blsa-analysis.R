library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(minfi)
library(quadprog)
library(methylCC)
library(here)

################################################
### Target Data: 806 450K DNA methylation samples with known cell type proportions
################################################

workingDir_blsa2017 <- here("case-studies", "2017-blsa-450-cellsort")
dataPath_blsa <- "/users/shicks1/data/DNAm/2017-blsa-450-cellsort"

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
counts450K <- estimateCellCounts(RGset)
save(counts450K, file = file.path(dataPath_blsa, "ests/cellcounts/dataBLSA-minfi-estimateCellCounts.RData"))

# compare to compote::estimateCC method
est <- estimateCC(object = RGset, verbose = TRUE, epsilon = 0.01,
                  maxIter = 100, initParamMethod = "random")
save(est, file = file.path(dataPath_blsa, "ests/cellcounts/dataBLSA-hicks-estimateCC-410regions.RData"))
rm(est) 


