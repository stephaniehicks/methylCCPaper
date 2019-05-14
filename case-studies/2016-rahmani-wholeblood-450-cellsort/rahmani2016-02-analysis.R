library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(minfi)
library(quadprog)
library(methylCC)


################################################
### Target Data: 
###       * 95 450K DNA methylation samples with known cell type proportions
################################################
workingDir_rahmani2016 <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/2016-rahmani-wholeblood-450-cellsort"
dataPath <- "/users/shicks1/data/GEO/GSE77716"

# Load 95 samples run on 450K DNA methylation array
#   - see 01_2016-rahmani-create-data-object.R for how Mset_rahmani2016.rds was created
galaMset <- readRDS(file.path(dataPath, "Mset_rahmani2016.rds"))


# Select only the 78 files used in Rahamani et al. 2016
#  - see README.md (https://github.com/cozygene/refactor) says Figure 1 in Rahamani et al. 2016 was generated with these 78 samples for which they had cell counts at the time
#  - data available from https://raw.githubusercontent.com/cozygene/refactor/master/assets/Sample_ids_for_Fig1.txt
#  - downloaded the 'Sample_ids_for_Fig1.txt' file and stored in the data/ folder
keepNames <- readr::read_table(file.path(workingDir_rahmani2016, 
                                         "data/Sample_ids_for_Fig1.txt"), col_names = FALSE)
galaMset <- galaMset[, match(keepNames$X1, colnames(galaMset) )]

# sanity check: see if cell composition information is available for these 78 samples
table(pData(galaMset)$cc_available) # yes

# read in beta values
pd <- pData(galaMset)
grObject <- granges(galaMset)
betaValues <- getBeta(galaMset, type = "Illumina")

# extract true cell type proportions
trueProps <- as.data.frame(pd[, names(pd) %in% 
            c("cc_baso", "cc_eso", "cc_neu", "cc_lymph", "cc_mono")]) # 46:50
trueProps$cc_gran <- rowSums(as.matrix(pd[, names(pd) %in% 
                                c("cc_baso", "cc_eso", "cc_neu")])) # 46:48
trueProps <- trueProps*.01
colnames(trueProps) <- c("Baso", "Eos", "Neu", "Lymph", "Mono", "Gran")


##### find regions
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

table(output$regions_all$cellType, output$regions_all$L, 
      output$regions_all$dmr_status, output$regions_all$status)

##### compare to houseman approach using 600 Houseman CpGs 
counts450KHouseman <- 
  minfi:::projectCellType(
          betaValues[rownames(output$cpgs_houseman),], output$cpgs_houseman)

counts450KHouseman <- as.data.frame(counts450KHouseman)

# boxplot(counts450KHouseman, y = c(0,1))
counts450KHouseman$Lymph <- 
  rowSums(counts450KHouseman[, colnames(counts450KHouseman) %in% 
                               c("CD4T", "CD8T", "NK", "Bcell")])
counts450KHouseman <- counts450KHouseman[, colnames(counts450KHouseman) %in% 
                                           c("Gran", "Mono", "Lymph")]

M1.rmse = mean(sqrt(colMeans((trueProps[,c(4,5,6)] - 
                                counts450KHouseman[,c(3,2,1)])^2)))
M1.rmse


##### run methylCC
set.seed(12345)
est <- estimatecc(object = galaMset, 
                  find_dmrs_object = output, 
                  verbose = TRUE, epsilon = 0.01, 
                  max_iter = 100, init_param_method = "random")
counts450KHicks <- cell_counts(est)

counts450KHicks$Lymph <- 
  rowSums(counts450KHicks[, colnames(counts450KHicks) %in% c("CD4T", "CD8T", "NK", "Bcell")])
counts450KHicks <- counts450KHicks[, which(colnames(counts450KHicks) %in% 
                                                 c("Gran", "Mono", "Lymph"))]
rownames(counts450KHicks) <- rownames(counts450KHouseman)

M2.rmse = mean(sqrt(colMeans((trueProps[,c(4,5,6)] - 
                                counts450KHicks[,c(3,2,1)])^2)))
M2.rmse


##### plot cell composition estimates
df.trueProps <- tidyr::gather(cbind("samples" = rownames(counts450KHouseman), 
                                    trueProps[,c(4:6)]), celltype, est, -samples)

df.Houseman <- tidyr::gather(cbind("samples" = rownames(counts450KHouseman), 
                                   counts450KHouseman), celltype, est, -samples)
df.Hicks <- tidyr::gather(cbind("samples" = rownames(counts450KHicks), 
                                counts450KHicks), celltype, est, -samples)


dfcombined <- dplyr::full_join(df.trueProps, df.Houseman,  by = c("samples", "celltype"))
dfcombined <- dplyr::full_join(dfcombined, df.Hicks,  by = c("samples", "celltype"))
colnames(dfcombined) <- c("samples", "celltype", "Truth", "Houseman", "methylCC")
dfcombined$celltype <- factor(dfcombined$celltype, levels = c("Lymph", "Mono", "Gran"))
dfcombined <- tidyr::gather(dfcombined, model, est, -c(samples, celltype, Truth))
dfcombined$model <- factor(dfcombined$model, levels = c("Houseman", "methylCC"))

gmat <- dfcombined %>% 
  ggplot(aes(x=Truth, y = est)) +  # 
  geom_abline(intercept = 0, slope = 1) + xlim(0,1) + ylim(0,1) +
  facet_grid(model~.) + geom_point(aes(color=celltype)) +
  theme(legend.position="top") + 
  xlab("True cell composition (measured with flow cytometry)") + 
  ylab("Model-based cell composition estimates") + 
  labs(title = "Whole blood samples measured on\nIllumina 450K platform (N=78, Rahmani et al. 2016)", 
       color = "Cell type") + 
  scale_color_discrete(name = "Cell type", 
                       labels = c("Lymphocytes", "Monocytes", "Granulocytes"))

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

pdf(file.path(workingDir_rahmani2016, 
              "figs/fig-rahmani2016-ests.pdf"), width = 6, height = 7)
print(gmat)
dev.off()
