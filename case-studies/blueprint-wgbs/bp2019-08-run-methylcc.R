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

####################################################
# run methylcc
output_nogr <- find_dmrs(gr_target = NULL,
                         include_cpgs = FALSE, include_dmrs = TRUE)

counts_methylcc <- NULL
for(smps in fnames)
{
  mes <- "[load bsseq data] Sample %s"
  message(sprintf(mes,smps))
  BS <- readRDS(file = file.path(dataPath, "files_bsseq",
                                 paste0("blueprint_blood_bs_", smps)))

  print("[starting parameter estimation]")
  set.seed(12345)
  est_nogr <- estimatecc(object = BS,
                         find_dmrs_object = output_nogr,
                         verbose = TRUE, epsilon = 0.01,
                         max_iter = 100, init_param_method = "random")

  counts.WGBS.Hicks <- cell_counts(est_nogr)
  rownames(counts.WGBS.Hicks) <- pData(BS)$sample_name
  counts_methylcc <- rbind(counts_methylcc, counts.WGBS.Hicks)
}

saveRDS(counts_methylcc,
        file = file.path(dataPath, "blueprint_blood_cc_methylcc.RDS"))


####################################################
# run methylcc (USING ONLY CPGs)
output_nogr <- find_dmrs(gr_target = NULL, 
                         include_cpgs = TRUE, include_dmrs = FALSE)

counts_methylcc <- NULL
for(smps in fnames)
{
  mes <- "[load bsseq data] Sample %s"
  message(sprintf(mes,smps))
  BS <- readRDS(file = file.path(dataPath, "files_bsseq",
                                 paste0("blueprint_blood_bs_", smps)))
  
  print("[starting parameter estimation]")
  set.seed(12345)
  est_nogr <- estimatecc(object = BS,
                         find_dmrs_object = output_nogr,
                         verbose = TRUE, epsilon = 0.01,
                         max_iter = 100, init_param_method = "random")
  
  counts.WGBS.Hicks <- cell_counts(est_nogr)
  rownames(counts.WGBS.Hicks) <- pData(BS)$sample_name
  counts_methylcc <- rbind(counts_methylcc, counts.WGBS.Hicks)
}

saveRDS(counts_methylcc, 
        file = file.path(dataPath, "blueprint_blood_cc_methylcc-only-cpgs.RDS"))
