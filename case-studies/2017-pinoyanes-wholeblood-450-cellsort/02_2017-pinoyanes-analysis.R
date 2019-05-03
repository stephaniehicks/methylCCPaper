# library(tidyr)
# library(dplyr)
library(ggplot2)
library(minfi)
library(quadprog)
library(bsseq)
library(GenomicRanges)
library(here)

################################################
### Train Data: 
###       * load 35 450K DNA methylation samples with known cell type proportions
###       * identify top differentially methylated 600 CpGs using `estimateCellCounts.modified()` on `RGset_target450k` 
#   - Returns 600 Houseman CpGs. 

################################################
workingDir_flowSort <- "/users/shicks1/data/DNAm/FlowSortedBlood450k"

library(FlowSorted.Blood.450k)
# FlowSorted.Blood.450k <- updateObject(FlowSorted.Blood.450k)
# Mset_train450k <- preprocessQuantile(FlowSorted.Blood.450k)
# Mset_train450k <- mapToGenome(Mset_train450k, mergeManifest = FALSE)
# save(Mset_train450k, file=file.path(workingDir_flowSort, "FlowSorted-Mset-Quantile.rda"))
load(file.path(workingDir_flowSort, "FlowSorted-Mset-Quantile.rda"))

Mset_train450k


################################################
### Target Data: 
###       * 95 450K DNA methylation samples with known cell type proportions
################################################
workingDir_PinoYanes2015 <- here("case-studies", "2017-pinoyanes-wholeblood-450-cellsort")
dataPath <- "/users/shicks1/data/GEO/GSE77716"

# Load 95 samples run on 450K DNA methylation array
#   - see dataFromLPinoYanes2015-450k-loadData.R for how Mset_PinoYanes2015.rds was created
galaMset <- readRDS(file.path(dataPath, "Mset_PinoYanes2015.rds"))

# Select only the 78 files used in Rahamani et al. 2016
#  - see README.md (https://github.com/cozygene/refactor) says Figure 2 in Rahamani et al. 2016 was generated with these 78 samples for which they had cell counts at the time
keepNames <- readr::read_table("https://raw.githubusercontent.com/cozygene/refactor/master/assets/Sample_ids_for_Fig2.txt", 
                                col_names = FALSE)
galaMset <- galaMset[, match(keepNames$X1, colnames(galaMset) )]

# sanity check: see if cell composition information is available for these 78 samples
table(pData(galaMset)$cc_available) # yes

# read in beta values
pd <- pData(galaMset)
grObject <- granges(galaMset)
betaValues = getBeta(galaMset)

# extract true cell type proportions
trueProps <- as.data.frame(pd[, 46:50])
trueProps$cc_gran <- rowSums(as.matrix(pd[, 46:48]))
trueProps <- trueProps*.01
colnames(trueProps) <- c("Baso", "Eos", "Neu", "Lymph", "Mono", "Gran")

##### run Houseman method 
devtools::load_all()

# extract out 600 Houseman CpGs (do not search for regions)
output <- findDMRs(Mset_train450k = Mset_train450k, verbose = TRUE, 
                   gr_target = granges(galaMset),
                   numRegions = 50, 
                   includeCpGs = TRUE, includeDMRs = FALSE, 
                   betaCutoff = 0.2, 
                   pairwiseComparison = FALSE)
dim(output$housemanCpGs)

counts450KHouseman <- minfi:::projectCellType(
          betaValues[rownames(output$housemanCpGs),], 
          output$housemanCpGs)

counts450KHouseman <- as.data.frame(counts450KHouseman)

# boxplot(counts450KHouseman, y = c(0,1))
counts450KHouseman$Lymph <- 
  rowSums(counts450KHouseman[, colnames(counts450KHouseman) %in% 
                               c("CD4T", "CD8T", "NK", "Bcell")])
counts450KHouseman <- counts450KHouseman[, colnames(counts450KHouseman) %in% 
                                           c("Gran", "Mono", "Lymph")]






##### run methylCC

# identify DMRs 
regions <- findDMRs(Mset_train450k = Mset_train450k, verbose=TRUE, 
                      gr_target=granges(galaMset), 
                      includeCpGs = TRUE, includeDMRs = TRUE, 
                      numRegions=50, numProbes=50,
                      betaCutoff = 0.2, 
                      pairwiseComparison = FALSE)

est <- estimateCC(object = galaMset, findRegions = FALSE, 
                  verbose = TRUE, epsilon = 0.01, 
                  maxIter = 100, initParamMethod = "random")
counts450KHicks <- cellcounts(est)

counts450KHicks$Lymph <- 
  rowSums(counts450KHicks[, colnames(counts450KHicks) %in% c("CD4T", "CD8T", "NK", "Bcell")])
counts450KHicks <- counts450KHicks[, which(colnames(counts450KHicks) %in% 
                                                 c("Gran", "Mono", "Lymph"))]
rownames(counts450KHicks) <- rownames(counts450KHouseman)



##### plot cell composition estimates
df.trueProps <- tidyr::gather(cbind("samples" = rownames(counts450KHouseman), 
                                    trueProps[,c(4:6)]), celltype, est, -samples)

df.Houseman <- tidyr::gather(cbind("samples" = rownames(counts450KHouseman), 
                                   counts450KHouseman), celltype, est, -samples)
df.Hicks <- tidyr::gather(cbind("samples" = rownames(counts450KHicks), 
                                counts450KHicks), celltype, est, -samples)


M1.rmse = mean(sqrt(colMeans((trueProps[,c(4,5,6)] - 
                                counts450KHouseman[,c(3,2,1)])^2)))
M1.rmse

M2.rmse = mean(sqrt(colMeans((trueProps[,c(4,5,6)] - 
                                counts450KHicks[,c(3,2,1)])^2)))
M2.rmse

dfcombined <- dplyr::full_join(df.trueProps, df.Houseman,  by = c("samples", "celltype"))
dfcombined <- dplyr::full_join(dfcombined, df.Hicks,  by = c("samples", "celltype"))
colnames(dfcombined) <- c("samples", "celltype", "Truth", "Houseman", "Hicks")
dfcombined$celltype <- factor(dfcombined$celltype, levels = c("Lymph", "Mono", "Gran"))
dfcombined <- tidyr::gather(dfcombined, model, est, -c(samples, celltype, Truth))
dfcombined$model <- factor(dfcombined$model, levels = c("Houseman", "Hicks"))

dfcombined %>% 
  ggplot(aes(x=Truth, y = est, color=celltype)) + 
  geom_abline(intercept = 0, slope = 1) + xlim(0,1) + ylim(0,1) +
  facet_grid(~model) + geom_point() +
  xlab("True Cell Composition (measured with flow cytometry)") + 
  ylab("Model-based Cell Composition Estimates")

