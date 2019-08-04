library(rtracklayer)
library(bsseq)
library(liftOver)
library(GenomicRanges)
library(rhdf5) 
library(HDF5Array)

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

workingDir_blueprint <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs"
dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"

custom_table <- readRDS(file = file.path(dataPath,"blueprint_blood_custom_table.RDS"))
custom_table$file_name <- unlist(lapply(stringr::str_split(custom_table$file_path, "CNAG/"), 
                                        function(x){ x[[2]]}))
custom_table$biosource_short <- unlist(lapply(stringr::str_split(custom_table$biosource, "-"), 
                                              function(x){ x[[1]]}))

# smps <- levels(custom_table$sample_name)[4]
for(smps in levels(custom_table$sample_name))
{
  custom_table_call <- custom_table[custom_table$sample_name %in% smps &
                                      grepl("_calls.bs_call", custom_table$file_name), ]
  custom_table_cov  <- custom_table[custom_table$sample_name %in% smps &
                                      grepl("_calls.bs_cov", custom_table$file_name), ]
  
  gr_call <- rtracklayer::import(file.path(dataPath, "files_bigwig", 
                                           custom_table_call$file_name), 
                                 format = "bigWig")
  colnames(mcols(gr_call)) <-  "score_call" 
  
  gr_cov <- rtracklayer::import(file.path(dataPath, "files_bigwig", 
                                          custom_table_cov$file_name), 
                                format = "bigWig")
  colnames(mcols(gr_cov)) <-  "score_cov" 
  
  gr_call$score_cov <- gr_cov$score_cov
  rm(gr_cov)
  
  genome(gr_call) <- "GRCh38"
  gr_call_19 <- liftOver(gr_call, ch)
  rm(gr_call)
  
  gr_call_19 <- unlist(gr_call_19)
  genome(gr_call_19) <- "hg19"
  gr_call_19 <- unique(gr_call_19)
  
  mcols(gr_call_19)$score_meth <- 
    round(mcols(gr_call_19)$score_call * mcols(gr_call_19)$score_cov)
   
  if(smps == levels(custom_table$sample_name)[1]){
    gr_complete <- gr_call_19[1:2e5,]
    
    colnames(mcols(gr_complete))[2:3] <- paste(smps, c("cov", "meth"), sep="_")
    gr_complete <- gr_complete[,-1]
  }
  
  if(smps != levels(custom_table$sample_name)[1]){
    gr_call_19 <- gr_call_19[1:2e5,]
    gr_complete = sort(unique(c(gr_complete, granges(gr_call_19))))

    mcols(gr_complete)[match(gr_call_19, gr_complete), 
                       paste(smps, c("cov"), sep="_")] = gr_call_19$score_cov
    mcols(gr_complete)[match(gr_call_19, gr_complete), 
                       paste(smps, c("meth"), sep="_")] = gr_call_19$score_meth
  }
 print(smps)
}  
  
tmp <- as.matrix(mcols(gr_complete))
tmp[is.na(tmp)] <- 0

cov_ids <- grepl("cov", colnames(tmp))
meth_ids <- grepl("meth", colnames(tmp))
  
hdf5_cov_path <- file.path(dataPath, "files_bsseq", 
                           "blueprint_blood_cov.h5")
hdf5_cov <- writeHDF5Array(x = tmp[,cov_ids], 
                           filepath = hdf5_cov_path, 
                           name = "cov")

hdf5_meth_path <- file.path(dataPath, "files_bsseq", 
                            "blueprint_blood_meth.h5")
hdf5_meth <- writeHDF5Array(x = tmp[, meth_ids],
                            filepath = hdf5_meth_path, 
                            name = "meth")

saveRDS(granges(gr_complete), 
        file = file.path(dataPath, "files_bsseq",
                         "blueprint_blood_granges.RDS"))

