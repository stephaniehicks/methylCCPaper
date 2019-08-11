library(GenomicRanges)
library(rhdf5) 
library(HDF5Array)
library(bsseq)
library(BiocSingular) 
library(BiocParallel)
workers <- 1
register(MulticoreParam(progressbar = TRUE, workers = workers))

# DelayedArray:::set_verbose_block_processing(TRUE)
# DelayedArray:::set_verbose_block_processing(TRUE)
# getAutoBlockSize()
# block_size <- 5e6
# setAutoBlockSize(block_size)

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
bs <- bs[!(seqnames(bs) %in% c("chrY", "chrX", "chrM")), ]

Sys.time() 
print("transpose")
# meth <- log(t(getCoverage(bs, type = "M")[seqnames(bs) %in% c("chr1"), ] / 
#           getCoverage(bs, type = "Cov")[seqnames(bs) %in% c("chr1"), ]) + 1)
meth  <- getCoverage(bs, type = "M")[seqnames(bs) %in% c("chr1"), ]
Sys.time() 


col_means <- DelayedMatrixStats::rowMeans2(meth)
col_sds <- DelayedMatrixStats::rowSds(meth)

Sys.time() 
print("runPCA")
set.seed(1234)
s <- runPCA(t(meth), rank = 2,
            center = col_means, 
            scale = col_sds,
            BSPARAM=IrlbaParam())
Sys.time() 

# print("save PCs")
# saveRDS(s$x,
#         file = file.path(dataPath, "files_pca",
#                          "blueprint_blood_pca_irlba.RDS"))

library(dplyr)
library(ggplot2)

dat <- data.frame(as.data.frame(pData(bs)), s$x[,1:2]) %>%
  mutate(begin = lubridate::dmy_hms(first_submission_date)) %>%
  mutate(year = factor(lubridate::year(begin)))

p1 <- dat %>% 
        ggplot(aes(x = PC1, y = PC2, color = cell_type)) + 
          geom_point(size=2) + xlab("Principal Component 1") +
          ylab("Principal Component 2") +
          labs(title = "BLUEPRINT reference methylomes (Chr 1)\n(N=44 WGBS samples colored by cell type)", 
               color = "Cell type") + 
          theme(legend.position = "top", legend.justification= "center") + 
          scale_color_discrete(name = "Cell type", 
                  labels = c("B-cells", "CD4 T-cells", "CD8 T-cells", 
                             "Granulocytes","Monocytes", "Natural killer cells"))


p2 <- dat %>%
        ggplot(aes(x = PC1, y = PC2, color = year)) + 
          geom_point(size=2) + xlab("Principal Component 1") +
          ylab("Principal Component 2") +
          labs(title = "BLUEPRINT reference methylomes (Chr 1)\n(N=44 WGBS samples colored by year processed)", 
              color = "Year") + 
          theme(legend.position = "top", legend.justification= "center")# + 
          # scale_color_discrete(name = "Cell type", 
          #         labels = c("B-cells", "CD4 T-cells", "CD8 T-cells", 
          #                   "Granulocytes","Monocytes", "Natural killer cells"))
  
gmat <- cowplot::plot_grid(p1, p2)

pdf(file.path(workingDir_blueprint, "figs/fig-blueprint-wgbs-pca.pdf"), width=11, height = 5)
print(gmat)
dev.off()
