# load libraries
library(minfi)
library(bsseq)
library(GEOquery)
library(GenomicRanges)

dataPath <- "/users/shicks1/data/GEO/GSE42861"
setwd(dataPath)

################################################################
#### Target Data: 686 blood samples run on 450K (without known cell type proportions)
####            ** 450K Samples ** 
################################################################
message("[dataFromLiu2013-450k] Downloading 450k data.")

# # download all supplemental files (including CEL files)
# dataPath <- "/users/shicks1/data/GEO"
# setwd(dataPath)
# getGEOSuppFiles("GSE42861", makeDirectory = TRUE, baseDir = dataPath)
# untar("GSE42861/GSE42861_RAW.tar", exdir="GSE42861/GSE42861_RAW")
cels <- list.files(file.path(dataPath, "GSE42861_RAW/"), pattern = ".idat")
tmp = paste0("GSE42861_RAW/", unique(gsub(".....idat.gz", "", cels)))

message("Downloading pData.")
gse <- getGEO("GSE42861",GSEMatrix=TRUE, destdir = dataPath)
pd <- pData(gse[[1]])
save(pd, file = file.path(dataPath, "pd_Liu2013.rds"))
rm(gse)

message("Creating RGset.")
RGset <- read.metharray(file.path(dataPath, tmp), 
                         extended = FALSE, verbose = FALSE, force = FALSE)
pData(RGset) <- pd

message("Saving RGset.")
# save RGset object
saveRDS(RGset, file=file.path(dataPath, "RGset_Liu2013.rds"))

# Load 689 samples run on 450K DNA methylation array
# RGset <- readRDS(file.path(dataPath, "RGset_Liu2013.rds")) # use this for minfi::estimateCellComp
RGset <- updateObject(RGset)

Sys.time()
Mset <- preprocessIllumina(RGset)
Sys.time()
Mset <- mapToGenome(Mset, mergeManifest = FALSE)

message("Save Mset object.")
Sys.time()
saveRDS(Mset, file=file.path(dataPath, "Mset_Liu2013.rds"))

