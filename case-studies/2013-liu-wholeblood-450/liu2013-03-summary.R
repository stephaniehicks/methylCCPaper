
library(ggplot2)
library(cowplot)

#########################################################
### Target Data: 686 450K DNA methylation samples 
###              (without known cell type proportions)
##########################################################

workingDir_liu2013 <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/2013-liu-wholeblood-450"
dataPath <- "/users/shicks1/data/GEO/GSE42861"

# see liu2013-01-create-data-object.R for how 
#   * RGset_Liu2013.rds
#   * Mset_Liu2013.rds 
# were created.

counts450KHouseman <- read.csv(file.path(workingDir_liu2013, 
                          "ests/450k_minfi-estimateCellCounts.csv"), 
                          row.names = 1)

counts450KHicks <- write.csv(file.path(workingDir_liu2013, 
                          "ests/450k_methylCC-estimatecc.csv"),
                          row.names = 1)

df.Houseman <- tidyr::gather(cbind("samples" = rownames(counts450KHouseman), 
                                   data.frame(counts450KHouseman)), celltype, est, -samples)
df.Hicks <- tidyr::gather(cbind("samples" = rownames(counts450KHicks), 
                                data.frame(counts450KHicks)), celltype, est, -samples)


dfcombined <- dplyr::full_join( df.Houseman, df.Hicks, by = c("samples", "celltype"))
colnames(dfcombined) <- c("samples", "celltype", "Houseman", "methylCC")
dfcombined$model <- factor(dfcombined$model, levels = c("Houseman", "methylCC"))

gmat <- dfcombined %>% 
  ggplot(aes(x=Houseman, y = methylCC)) +  # 
  geom_abline(intercept = 0, slope = 1) + xlim(0,1) + ylim(0,1) +
  geom_point(aes(color=celltype)) +
  theme(legend.position="top") + 
  xlab("Houseman method") + ylab("methylCC") + 
  labs(title = "Model-based cell composition estimates from\nwhole blood samples (N=689, Li et al. 2013)", 
       color = "Cell type") + 
  # scale_color_discrete(name = "Cell type", 
  #                      labels = c("Lymphocytes", "Monocytes", "Granulocytes"))

dat_text <- data.frame(
  label = c(paste0("RMSE: ", round(M1.rmse, 4)),
            paste0("RMSE: ", round(M2.rmse, 4))),
  model   = c("Houseman", "methylCC"), 
  x = c(0, 0), 
  y = c(0.65, 0.65)
)
gmat <- gmat + geom_text(
  data    = dat_text,
  mapping = aes(x = x, y = y, label = label),
  hjust   = -0.1,
  vjust   = -1,
  size = 5, 
  fontface=2
)

pdf(file.path(workingDir_liu2013, 
              "figs/fig-liu2013-ests.pdf"), width = 6, height = 6)
print(gmat)
dev.off()


# to delete after this point
# Mset <- readRDS(file=file.path(dataPath, "Mset_Liu2013.rds"))

# pd = colData(Mset)
# grObject <- granges(Mset)
# betaValues = getBeta(Mset, type = "Illumina")

# message("Pick Target Positions")
# Sys.time()
# keepDMRs = methylCC:::pickTargetPositions(targetgRanges=grObject, targetObject = betaValues,
#                                           dmpRegions=celltypeSpecificDMRegions)$dmpGRanges

# zmatSub <- subsetByOverlaps(celltypeSpecificDMRegions, keepDMRs, type = "equal")
# rm(betaValues)


# ## run minfi::estimateCellCounts()
# counts.450k.Houseman <- minfi::estimateCellCounts(RGset)
# write.csv(counts.450k.Houseman,
#          file = here("case-studies", "2013-liu-wholeblood-450", "ests", 
#                      "cellCountsEsts", "450k_minfi-estimateCellCounts.csv"),
#          row.names = TRUE)

# counts.450k.Houseman$Tcell <- rowSums(counts.450k.Houseman[,1:3])
# counts.450k.Houseman <- counts.450k.Houseman[,c(7,4:6)]
# boxplot(counts.450k.Houseman)
# df.Houseman <- gather(cbind("samples" = rownames(counts.450k.Houseman), counts.450k.Houseman), celltype, est, -samples)



# ## Create a processed methylcc object
# ymat <- mcols(keepDMRs); colnames(ymat) <- sampleNames(Mset)
# rm(Mset)
# n <- ncol(mcols(keepDMRs))
# ids <- colnames(mcols(zmatSub))
# zmat <- as.matrix(as.data.frame(mcols(zmatSub)))

# dat <- list(ymat = ymat, n = n, ids = ids, zmat = zmat,
#             grObject = grObject, keepDMRs = keepDMRs)

# message("Saving methylCC Processed Data.")
# Sys.time()
# saveRDS(dat, file=file.path(workingDir_liu2013, "methylCCProcessed_Liu2013.rds"))

# dat <- readRDS(file = file.path(workingDir_liu2013, "methylCCProcessed_Liu2013.rds"))

## run methylCC::estimateCC()
# Sys.time()
# est <- estimateCC(object = as.matrix(dat$ymat), regionMat = as.matrix(dat$zmat),
#                   initParamMethod = "random", epsilon = 0.1, 
#                   verbose = TRUE, maxIter = 50)
# Sys.time()
# counts.450k.Hicks <- cellcounts(est)

# write.csv(counts.450k.Hicks,
#           file = here("case-studies", "2013-liu-wholeblood-450", "ests", 
#                      "cellCountsEsts", "450k_methylCC-estimateCC.csv"),
#          row.names = TRUE)

# counts.450k.Hicks <- read.csv(here("case-studies", "2013-liu-wholeblood-450", "ests", 
#                                    "cellCountsEsts", "450k_methylCC-estimateCC.csv"),  
#                               row.names = 1)
# boxplot(counts.450k.Hicks)
# df.Hicks <- gather(cbind("samples" = rownames(counts.450k.Hicks), counts.450k.Hicks), celltype, est, -samples)



##### combine into main figure 

# dfcombined <- full_join(df.Houseman, df.Hicks,  by = c("samples", "celltype"))
# dfcombined$celltype <- factor(dfcombined$celltype, levels = c("Bcell", "Tcell", "Mono", "Gran"))

# gmat <- dfcombined %>% 
#   ggplot(aes(x=est.x, y = est.y)) + 
#     geom_abline(intercept = 0, slope = 1) + xlim(0,1) + ylim(0,1) +
#     facet_grid( ~celltype) + geom_point() +
#     xlab("Reference-based model (Houseman)") + ylab("Proposed model (Hicks)") + 
#     labs(title = "Model-based cell composition estimates from whole blood samples (n=689, Li et al. 2013)")
# pdf(here("case-studies", "2013-liu-wholeblood-450", "figs", 
#          "array_HousemanHicks_scatterplot_Liu2013_noTruth.pdf"), 
#     width = 12, height = 4)
# print(gmat)
# dev.off()

# gmat <- dfcombined %>% 
#  ggplot(aes(x=est.x, y = est.y, color = celltype)) + 
#  geom_abline(intercept = 0, slope = 1) + xlim(0,1) + ylim(0,1) +
#  scale_color_discrete(name = "Cell type") +
#  geom_point() + theme(legend.position = "top", legend.justification= "center") + 
#  xlab("Houseman method") + ylab("methylCC") + 
#  labs(title = "Model-based cell composition estimates from\nwhole blood samples (N=689, Li et al. 2013)")


# pdf(here("case-studies", "2013-liu-wholeblood-450", "figs", 
#         "array_HousemanHicks_scatterplot_Liu2013_noTruth_notfaceted.pdf"), 
#     width = 6, height = 6)
# print(gmat)
# dev.off()


