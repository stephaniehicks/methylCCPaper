library(GenomicRanges)
library(rhdf5) 
library(HDF5Array)
library(bsseq)
library(dmrseq)

library(BiocParallel)
workers <- 1
register(MulticoreParam(workers))

DelayedArray:::set_verbose_block_processing(TRUE)
DelayedArray:::set_verbose_block_processing(TRUE)

workingDir_blueprint <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs"
dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"
scratchPath <- "/fastscratch/myscratch/shicks1"


getAutoBlockSize()
block_size <- 5e6 
setAutoBlockSize(block_size)

## Search for regions using `dmrseq`
regions <- dmrseq(bs=bs, 
                  testCovariate="cell_type", 
                  minNumRegion = 3,
                  cutoff = 0.1, verbose = TRUE)

show(regions)

table(regions$L)

# get annotations for hg18
annoTrack <- getAnnot("hg18")

plotDMRs(bs, regions=regions[1,], testCovariate="CellType",
         annoTrack=annoTrack)

