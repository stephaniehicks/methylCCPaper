library(bsseq)
library(methylCC)
library(dplyr)

workingDir_blueprint <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs"
dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"

custom_table <- readRDS(file = file.path(dataPath,"blueprint_blood_custom_table.RDS"))
custom_table$file_name <- unlist(lapply(stringr::str_split(custom_table$file_path, "CNAG/"), 
                                        function(x){ x[[2]]}))
custom_table$biosource_short <- unlist(lapply(stringr::str_split(custom_table$biosource, "-"), 
                                              function(x){ x[[1]]}))
fnames <- list.files(file.path(dataPath, "files_bsseq/"), pattern = "blood_bs")
fnames <- unlist(lapply(stringr::str_split(fnames, "_"), 
                        function(x){ x[[4]]}))


pull_sH<- function(x,y){ 
  as.data.frame(findOverlaps(x, y)) %>% 
    dplyr::group_by(queryHits) %>% 
    dplyr::slice(1) %>% 
    dplyr::pull(subjectHits)
}




####################################################
# run methylCC::find_dmrs to find DMRs
output_nogr <- find_dmrs(gr_target = NULL,
                         include_cpgs = TRUE, include_dmrs = TRUE)

output_nogr$regions_all
output_nogr$y_regions
output_nogr$profiles

####################################################
# 
output_wgbs_profiles = output_450k_profiles <- vector("list", length(fnames))
pdat <- NULL
for(smps in seq_along(fnames))
{
  mes <- "[load bsseq data] Sample %s"
  message(sprintf(mes,smps))
  BS <- readRDS(file = file.path(dataPath, "files_bsseq",
                                 paste0("blueprint_blood_bs_", fnames[smps])))
  pdat <- rbind(pdat, pData(BS))
  cvg_targetbs <- getBSseq(BS, type = "Cov")
  M_targetbs <- getBSseq(BS, type = "M")
  output_wgbs_profiles[[smps]] <- .pick_target_positions(target_granges = granges(BS), 
                                     target_object = M_targetbs,
                                     target_cvg = cvg_targetbs, 
                                     dmp_regions = output_nogr$regions_all)$dmp_granges
  pull_subjectHits <- pull_sH(output_wgbs_profiles[[smps]], output_nogr$regions_all)
  output_450k_profiles[[smps]] <- output_nogr$profiles[pull_subjectHits, ]
} 

pdat$biosource <- 
  factor(pdat$biosource,
         levels = c("CD4-positive, alpha-beta T cell",
                    "CD8-positive, alpha-beta T cell", 
                    "mature neutrophil", 
                    "mature eosinophil", 
                    "CD38-negative naive B cell", 
                    "CD14-positive, CD16-negative classical monocyte", 
                    "cytotoxic CD56-dim natural killer cell"))

levels(pdat$biosource)

pdat$celltype <- factor(pdat$biosource_short)
levels(pdat$celltype) <- c("Mono", "Bcell", "CD4T", "CD8T", "NK", "Gran", "Gran")

x <- methylCC:::.splitit(pdat$celltype)

profiles_wgbs <- vector("list", length(levels(pdat$celltype)))
names(profiles_wgbs) <- levels(pdat$celltype)
for(ind in levels(pdat$celltype)){ 

    ct_granges <- Reduce(intersect, output_wgbs_profiles[x[[ind]]])
    data_profiles <- data.frame(matrix(nrow = length(ct_granges), ncol =  2))
    pull_subjectHits <- pull_sH(ct_granges, output_nogr$regions_all)
    data_profiles[, 1] <- output_nogr$profiles[pull_subjectHits, ind]
    
    data_prof <- NULL
    for(ii in seq_along(x[[ind]])){
      pull_subjectHits <- pull_sH(ct_granges, output_wgbs_profiles[x[[ind]]][[ii]])
      data_prof <- cbind(data_prof, output_wgbs_profiles[x[[ind]]][[ii]][pull_subjectHits, ]$X)
    } 
    data_profiles[,2] <- rowMeans(as.matrix(data_prof))
    profiles_wgbs[[ind]] <- data_profiles
}

data_tidy <- plyr::ldply(profiles_wgbs)
saveRDS(data_tidy, 
        file = file.path(dataPath, "blueprint_celltype-450k-wgbs-withCpGs.RDS"))










####################################################
# run methylCC::find_dmrs to find DMRs
output_onlydmr <- find_dmrs(gr_target = NULL,
                         include_cpgs = FALSE, include_dmrs = TRUE)

output_onlydmr$regions_all
output_onlydmr$y_regions
output_onlydmr$profiles

####################################################
# 
output_wgbs_profiles = output_450k_profiles <- vector("list", length(fnames))
pdat <- NULL
for(smps in seq_along(fnames))
{
  mes <- "[load bsseq data] Sample %s"
  message(sprintf(mes,smps))
  BS <- readRDS(file = file.path(dataPath, "files_bsseq",
                                 paste0("blueprint_blood_bs_", fnames[smps])))
  pdat <- rbind(pdat, pData(BS))
  cvg_targetbs <- getBSseq(BS, type = "Cov")
  M_targetbs <- getBSseq(BS, type = "M")
  output_wgbs_profiles[[smps]] <- .pick_target_positions(target_granges = granges(BS), 
                                                         target_object = M_targetbs,
                                                         target_cvg = cvg_targetbs, 
                                                         dmp_regions = output_onlydmr$regions_all)$dmp_granges
  pull_subjectHits <- pull_sH(output_wgbs_profiles[[smps]], output_onlydmr$regions_all)
  output_450k_profiles[[smps]] <- output_onlydmr$profiles[pull_subjectHits, ]
} 

pdat$biosource <- 
  factor(pdat$biosource,
         levels = c("CD4-positive, alpha-beta T cell",
                    "CD8-positive, alpha-beta T cell", 
                    "mature neutrophil", 
                    "mature eosinophil", 
                    "CD38-negative naive B cell", 
                    "CD14-positive, CD16-negative classical monocyte", 
                    "cytotoxic CD56-dim natural killer cell"))

levels(pdat$biosource)

pdat$celltype <- factor(pdat$biosource_short)
levels(pdat$celltype) <- c("Mono", "Bcell", "CD4T", "CD8T", "NK", "Gran", "Gran")

x <- methylCC:::.splitit(pdat$celltype)

profiles_wgbs <- vector("list", length(levels(pdat$celltype)))
names(profiles_wgbs) <- levels(pdat$celltype)
for(ind in levels(pdat$celltype)){ 
  
  ct_granges <- Reduce(intersect, output_wgbs_profiles[x[[ind]]])
  data_profiles <- data.frame(matrix(nrow = length(ct_granges), ncol =  2))
  pull_subjectHits <- pull_sH(ct_granges, output_onlydmr$regions_all)
  data_profiles[, 1] <- output_onlydmr$profiles[pull_subjectHits, ind]
  
  data_prof <- NULL
  for(ii in seq_along(x[[ind]])){
    pull_subjectHits <- pull_sH(ct_granges, output_wgbs_profiles[x[[ind]]][[ii]])
    data_prof <- cbind(data_prof, output_wgbs_profiles[x[[ind]]][[ii]][pull_subjectHits, ]$X)
  } 
  data_profiles[,2] <- rowMeans(as.matrix(data_prof))
  profiles_wgbs[[ind]] <- data_profiles
}

data_tidy <- plyr::ldply(profiles_wgbs)
# saveRDS(data_tidy, 
#         file = file.path(dataPath, "blueprint_celltype-450k-wgbs-noCpGs.RDS"))
# saveRDS(data_tidy, 
#         file = file.path(dataPath, "blueprint_celltype-450k-wgbs-withCpGs.RDS"))

rm(data_tidy)



