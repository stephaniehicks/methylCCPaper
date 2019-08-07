library(bsseq)
library(GenomicRanges)

library(rhdf5) 
library(HDF5Array)
rhdf5::h5disableFileLocking()

library(BiocParallel)
workers <- 5
DelayedArray::setAutoBPPARAM(MulticoreParam(workers))

DelayedArray:::set_verbose_block_processing(TRUE)
DelayedArray:::set_verbose_block_processing(TRUE)

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
### Applying quality control
#######################################

# load GRanges object
message("load GRanges object")
gr_complete <- readRDS(file = file.path(dataPath, "files_bsseq",
                                        "blueprint_blood_gr_01.RDS"))

# load Cov and M datasets
message("load Cov and M datasets")
hdf5_bs_path <- file.path(dataPath, "files_bsseq",
                          "blueprint_blood_01.h5")
hdf5_cov  <- HDF5Array(filepath = hdf5_bs_path, name = "cov")
hdf5_meth <- HDF5Array(filepath = hdf5_bs_path, name = "meth")

getAutoBlockSize()
block_size <- 5e6
setAutoBlockSize(block_size)

# dmrseq requires at least 1 read per group
# min group size is 4 samples, so we require less than or equal to three 0s per group
message("Subsetting CpGs")
# keep_ids <- which(DelayedMatrixStats::rowSums2(hdf5_cov > 3) == 44)
keep_ids <- which(DelayedArray::rowSums(hdf5_cov > 3) == 44)
keep_ids <- sort(keep_ids)
length(keep_ids)

hdf5_cov_sub <- hdf5_cov[keep_ids,]
hdf5_meth_sub <- hdf5_meth[keep_ids,]
gr_complete_sub <- gr_complete[keep_ids,]

message("save granges only from gr_complete_sub")
saveRDS(gr_complete_sub,
        file = file.path(dataPath, "files_bsseq",
                         "blueprint_blood_gr_02.RDS"))

# getHDF5DumpChunkDim(dim(hdf5_cov))
# getHDF5DumpChunkDim(dim(hdf5_cov_sub))

# realize new files on disk
hdf5_bs_path_02 <- file.path(dataPath, "files_bsseq",
                             "blueprint_blood_02.h5")

message("Saving Cov file.")
hdf5_cov_02  <- writeHDF5Array(x = hdf5_cov_sub, filepath = hdf5_bs_path_02,
                               chunkdim = c(2e3, ncol(hdf5_cov_sub)),
                               name = "cov", verbose=TRUE)
message("Saving Meth file.")
hdf5_meth_02 <- writeHDF5Array(x = hdf5_meth_sub, filepath = hdf5_bs_path_02,
                               chunkdim = c(2e3, ncol(hdf5_meth_sub)),
                               name = "meth", verbose=TRUE)

