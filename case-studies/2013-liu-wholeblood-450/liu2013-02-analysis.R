library(minfi)
library(quadprog)
library(methylCC)

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

# # Load 689 samples run on 450K DNA methylation array
Sys.time()
RGset <- readRDS(file.path(dataPath, "RGset_liu2013.rds")) # use this for minfi::estimateCellComp
RGset <- updateObject(RGset)

pd <- colData(RGset)
colData(RGset) <- DataFrame(pd[,1:2])
colData(RGset)$title <- as.character(colData(RGset)$title)
Sys.time()

## run minfi::estimateCellCounts()
# counts450KHouseman <- minfi::estimateCellCounts(RGset)
# write.csv(counts450KHouseman,
#          file.path(workingDir_liu2013, 
#                    "ests/450k_minfi-estimateCellCounts.csv"),
#          row.names = TRUE)



## run methylCC::estimateCC()

## find regions
Sys.time()
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
Sys.time()

table(output$regions_all$cellType, output$regions_all$L, 
      output$regions_all$dmr_status, output$regions_all$status)

## run methylCC::estimatecc
Sys.time()
set.seed(12345)
est <- estimatecc(object = RGset, 
                  find_dmrs_object = output, 
                  verbose = TRUE, epsilon = 0.01, 
                  max_iter = 100, init_param_method = "random")
Sys.time()

counts450KHicks <- cell_counts(est)
rownames(counts450KHicks) <- rownames(counts450KHicks)

write.csv(counts450KHicks,
          file.path(workingDir_liu2013, 
                    "ests/450k_methylCC-estimatecc.csv"),
          row.names = TRUE)
