library(GenomicRanges)
library(rhdf5) 
library(HDF5Array)
library(bsseq)
library(dmrseq)
library("BiocParallel")
# register(MulticoreParam(4))

workingDir_blueprint <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs"
dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"

custom_table <- readRDS(file = file.path(dataPath,"blueprint_blood_custom_table.RDS"))
custom_table$file_name <- unlist(lapply(stringr::str_split(custom_table$file_path, "CNAG/"), 
                                        function(x){ x[[2]]}))
custom_table$biosource_short <- unlist(lapply(stringr::str_split(custom_table$biosource, "-"), 
                                              function(x){ x[[1]]}))

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

gr_complete <- readRDS(file = file.path(dataPath, "files_bsseq",
                                        "blueprint_blood_granges.RDS"))

bs <- BSseq(gr = gr_complete,
            M = hdf5_meth, 
            Cov = hdf5_cov,
            sampleNames = levels(custom_table$sample_name)[1:4])
pData(bs)$Condition <- factor(rep(c("control", "case"), each=2))
pryr::object_size(bs)

loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")==0) == 0)

bs.filtered <- bs[loci.idx,]

regions <- dmrseq(bs=bs.filtered, 
                  testCovariate="Condition", 
                  minNumRegion = 1,
                  cutoff = 0.1, verbose = TRUE)

show(regions)

table(regions$L)

