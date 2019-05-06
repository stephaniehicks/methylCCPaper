# load libraries
library(tidyr)
library(dplyr)
library(minfi)
library(bsseq)
library(methylCC) # on Github
library(GenomicRanges)

dataPath <- "/net/irizarryfs01/srv/export/irizarryfs01_backed_up/share_root/shicks/projects/methylCC/dataFromCarmona-RRBS-450k"

################################################################
#### Target Data: 12 blood samples run on 450K and RRBS
####
####            ** 450K Samples **
################################################################
message("[dataFromCarmona-RRBS-450k] Loading 450k data.")

load(file.path(dataPath, "mset450krawall.RData"))
mset450krawall <- updateObject(mset450krawall)
RGset_target450k <- mset450krawall
rm(mset450krawall)

# This is needed to run the estimateCellCounts() function
colData(RGset_target450k)$Sentrix_Position <- as.character(colData(RGset_target450k)$Sentrix_Position)

# add ID to pData()
library(stringr)
newnames <- sapply(pData(RGset_target450k)$sampleID, function(x){
  str_split(x, pattern = "_")[[1]][1] })
newnames[1:3] <- "D07"      # these are three tech reps from same individual
pData(RGset_target450k)$ID <- newnames

# create target450k datasets
Mset_target450k <- preprocessIllumina(RGset_target450k)
Mset_target450k <- mapToGenome(Mset_target450k, mergeManifest = FALSE)

pd_target450k = as.data.frame(pData(Mset_target450k))
p_target450k = getBeta(Mset_target450k , type = "Illumina")
gr_target450k <- granges(Mset_target450k)

mcols(gr_target450k) <- p_target450k


################################################################
#### Target Data: 12 blood samples run on 450K and WGBS
####
####              ** WGBS Samples **
################################################################
message("[dataFromCarmona-RRBS-450k] Loading RRBS data.")

# load(file.path(dataPath, "combined.data.Best12.RData"))
# dat <- combined.data.Best12 # list of three data.frames with 6million rows each (DO NOT print by itself!)
# rm(combined.data.Best12)
#
# dat$methyl <- round(((0.01) * dat$percentmethyl[,-c(1,2)]) * dat$numberreads[,-c(1,2)] )
# dat$methyl[is.na(dat$methyl)] <- 0 # remove NAs
# dat$numberreads[is.na(dat$numberreads)] <- 0 # remove NAs
#
# newnames <- sapply(colnames(dat$methyl), function(x){
#   paste(str_c(str_split(x, pattern = "\\.")[[1]][2:4], collapse = "_"), "R1.fastq", sep="_") })
# newnames[11:12] <- c("D07repA_RRBS", "D07repB_RRBS")
# colnames(dat$numberreads)[-c(1,2)] = colnames(dat$methyl) = newnames
#
# message("[dataFromCarmona-RRBS-450k] Creating BSseq object.")
# # create BSeq Object
# BS <- BSseq(pos = dat$numberreads$start,
#             chr = dat$numberreads$chr,
#             M = as.matrix(dat$methyl),
#             Cov = as.matrix(dat$numberreads[,-c(1,2)]))
#
# # add pData() object; create ID column
# newIDs <- sapply(colnames(dat$methyl), function(x){
#   str_split(x, pattern = "_")[[1]][1] })
# newIDs[11:12] <- c("D07", "D07")        # these are three tech reps from same individual
# pData(BS) <- DataFrame("sampleID" = colnames(dat$methyl), ID = newIDs)
# rm(dat)
# save(BS, file=file.path(dataPath, "bsseqkrawall.RData"))

load(file.path(dataPath, "bsseqkrawall.RData"))

# cvg_targetbs <- as.data.frame(getBSseq(BS, type = "Cov"))
# M_targetbs <- as.data.frame(getBSseq(BS, type = "M"))
# p_targetbs <- as.data.frame(getMeth(BS, type="raw")) # meth / cvg
# p_targetbs <- as.data.frame(p_targetbs)
# pd_targetbs <- pData(BS)
# gr_targetbs <- granges(BS)
#
# mcols(gr_targetbs) <- p_targetbs


# # length(gr_target450k) = 485512
# # length(gr_targetbs) = 6823620
# # How many CpGs overlap?
# intersectProbes = subsetByOverlaps(gr_target450k, gr_targetbs)
# intersectProbes # length(intersectProbes) = 142002

```



