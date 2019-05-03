# load libraries
library(minfi)
library(bsseq)
library(GEOquery)
library(GenomicRanges)
library(stringr)

dataPath <- "/users/shicks1/data/GEO/GSE77716"
setwd(dataPath)

################################################################
#### Target Data: 686 blood samples run on 450K
####
####            ** 450K Samples ** 
################################################################
# message("Downloading pData.")
# gse <- getGEO("GSE77716",GSEMatrix=TRUE, destdir = dataPath)
# pd <- pData(gse[[1]])
# rownames(pd) <- pd$title
# saveRDS(pd, file = file.path(dataPath, "pd_PinoYanes2015.rds"))
# rm(gse)
pd <- readRDS(file.path(dataPath, "pd_PinoYanes2015.rds"))
pd$cc_baso <- as.numeric(dplyr::na_if(pd$`gran_basophils:ch1`, "null"))
pd$cc_eso <- as.numeric(dplyr::na_if(pd$`gran_eosinophils:ch1`, "null"))
pd$cc_neu <- as.numeric(dplyr::na_if(pd$`gran_neutrophils:ch1`, "null"))
pd$cc_lymph <- as.numeric(dplyr::na_if(pd$`lymphocytes:ch1`, "null"))
pd$cc_mono <- as.numeric(dplyr::na_if(pd$`monocytes:ch1`, "null"))

# is cell composition information available? 
pd$cc_available <- !is.na(pd$cc_baso) & !is.na(pd$cc_eso) & !is.na(pd$cc_neu) & 
  !is.na(pd$cc_lymph) & !is.na(pd$cc_mono)

message("[dataFromPinoYanes2015-450k] Downloading 450k data.")
# # download all supplemental files (including CEL files)
# dataPath <- "/users/shicks1/data/GEO"
# setwd(dataPath)
# getGEOSuppFiles("GSE77716", makeDirectory = TRUE, baseDir = dataPath)
# gunzip("GSE77716/GSE77716_Matrix_signal_intensities_updated_2018-Jan08.tsv.gz")
dat <- readr::read_tsv(file.path(dataPath, "GSE77716_Matrix_signal_intensities_updated_2018-Jan08.tsv"))

meth <- as.matrix(dat[, seq(3, 1720, by = 3)])
colnames(meth) <- gsub(" Methylated Signal", "", colnames(meth)) 
meth <- meth[, match(pd$title, colnames(meth))]
rownames(meth) <- dat$ID_REF

unmeth <- as.matrix(dat[, seq(2, 1720, by = 3)])
colnames(unmeth) <- gsub(" Unmethylated Signal", "", colnames(unmeth)) 
unmeth <- unmeth[, match(pd$title, colnames(unmeth))]
rownames(unmeth) <- dat$ID_REF

# Create MethylSet
galaMset <- MethylSet(Meth = meth, Unmeth = unmeth, colData = pd)

# emailed the authors to ask about annotation and manifest
#   same as Bioc package `FlowSorted.Blood.450k`
annotation(galaMset) <- annotation(FlowSorted.Blood.450k::FlowSorted.Blood.450k)

# add genome information
galaMset <- mapToGenome(galaMset, mergeManifest = TRUE)

message("Saving Mset")
# save Mset object
saveRDS(galaMset, file=file.path(dataPath, "Mset_PinoYanes2015.rds"))

