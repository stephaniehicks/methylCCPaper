library(GenomicRanges)
library(rhdf5) 
library(HDF5Array)
library(bsseq)
library(dmrseq)
library(BiocSingular) 
library(BiocParallel)
workers <- 10
register(MulticoreParam(workers))

workingDir_blueprint <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs"
dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"
scratchPath <- "/fastscratch/myscratch/shicks1"

print("Reading in BSseq object")
Sys.time()
hdf5_bs_se_path <- file.path(dataPath, "files_bsseq_hdf5_se")
bs <- loadHDF5SummarizedExperiment(hdf5_bs_se_path)
Sys.time()

pryr::object_size(bs)

print("Creating design matrix")
xmat = model.matrix(~bs$cell_type - 1)
colnames(xmat) = levels(bs$cell_type)
pData(bs) <- cbind(pData(bs), xmat)

print("Identifying DMRs with dmrseq")
Sys.time()
regions <- dmrseq(bs=bs, 
                  testCovariate="CD8T", minNumRegion = 3,
                  cutoff = 0.2, verbose = TRUE)
Sys.time()

print("save DMRs found with dmrseq")
saveRDS(regions,
        file = file.path(dataPath, "files_dmrs",
                         "blueprint_blood_regions_dmrseq_cd8t.RDS"))


