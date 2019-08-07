library(rtracklayer)
library(bsseq)
library(liftOver)
library(GenomicRanges)
library(rhdf5) 
library(HDF5Array)

rhdf5::h5disableFileLocking()
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

workingDir_blueprint <-
  "/users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs"
dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"
scratchPath <- "/fastscratch/myscratch/shicks1"

custom_table <- readRDS(file = file.path(dataPath,"blueprint_blood_custom_table.RDS"))
custom_table$file_name <- unlist(lapply(stringr::str_split(custom_table$file_path, "CNAG/"),
                                        function(x){ x[[2]]}))
custom_table$biosource_short <- unlist(lapply(stringr::str_split(custom_table$biosource, "-"),
                                              function(x){ x[[1]]}))

#######################################
### Creating HDF5 files for data package
#######################################

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
    gr_complete <- gr_call_19

    colnames(mcols(gr_complete))[2:3] <- paste(smps, c("cov", "meth"), sep="_")
    gr_complete <- gr_complete[,-1]
  }

  if(smps != levels(custom_table$sample_name)[1]){
    gr_complete = unique(c(gr_complete, granges(gr_call_19)))

    mcols(gr_complete)[match(gr_call_19, gr_complete),
                       paste(smps, c("cov"), sep="_")] = gr_call_19$score_cov
    mcols(gr_complete)[match(gr_call_19, gr_complete),
                       paste(smps, c("meth"), sep="_")] = gr_call_19$score_meth
   }
 print(smps)
}

gr_complete <- sort(gr_complete)

# save everything as an RDS (for safety)
message("save gr_complete")
saveRDS(gr_complete,
        file = file.path(dataPath, "files_bsseq",
                         "blueprint_blood_gr_complete.RDS"))
# deleted this file to save space (see below)

# Create Cov and M datasets and save to a HDF5 file
tmp_colnames <- colnames(mcols(gr_complete))
cov_ids <- grepl("cov", tmp_colnames)
meth_ids <- grepl("meth", tmp_colnames)

hdf5_bs_path <- file.path(dataPath, "files_bsseq",
                          "blueprint_blood_00.h5")
rhdf5::h5disableFileLocking()
hdf5_cov <- writeHDF5Array(x = mcols(gr_complete)[, cov_ids],
                           filepath = hdf5_bs_path, name = "cov")
rhdf5::h5disableFileLocking()
hdf5_meth <- writeHDF5Array(x = mcols(gr_complete)[, meth_ids],
                           filepath = hdf5_bs_path, name = "meth")

# replace NAs with 0s (converts from an HDF5Array to a DelayedArray)
hdf5_cov[is.na(hdf5_cov)] <- 0
hdf5_meth[is.na(hdf5_meth)] <- 0

# realize new files on disk
hdf5_bs_path_clean <- file.path(dataPath, "files_bsseq",
                          "blueprint_blood_01.h5")
hdf5_cov <- writeHDF5Array(x = hdf5_cov,
                           filepath = hdf5_bs_path_clean, name = "cov")
hdf5_meth <- writeHDF5Array(x = hdf5_meth, 
                            filepath = hdf5_bs_path_clean, name = "meth")

message("save granges only from gr_complete")
gr_complete <- granges(gr_complete)
saveRDS(gr_complete,
        file = file.path(dataPath, "files_bsseq",
                         "blueprint_blood_gr_01.RDS"))

# cleaning up
file.remove(file.path(dataPath, "files_bsseq",
                      "blueprint_blood_gr_complete.RDS"))