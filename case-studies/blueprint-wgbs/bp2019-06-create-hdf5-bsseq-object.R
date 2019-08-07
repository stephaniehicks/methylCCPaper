library(bsseq)
library(GenomicRanges)

library(rhdf5) 
library(HDF5Array)
rhdf5::h5disableFileLocking()

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
### Creating final BSseq object
#######################################

# Create phenotype table
message("creating pheno table")
pheno_table <- NULL
for(smps in levels(custom_table$sample_name))
{
  custom_table_call <- custom_table[custom_table$sample_name %in% smps &
                                      grepl("_calls.bs_call", custom_table$file_name), ]
  custom_table_cov  <- custom_table[custom_table$sample_name %in% smps &
                                      grepl("_calls.bs_cov", custom_table$file_name), ]
  custom_table_call$file_path_cov <- as.character(custom_table_cov$file_path)
  custom_table_call$file_name_cov <- as.character(custom_table_cov$file_name)
  pheno_table <- rbind(pheno_table, custom_table_call)
}
table(pheno_table$biosource_short)

# add `cell_type`` to match the Flow.Sorted.Blood.450K labels
pheno_table$cell_type <-
  factor(plyr::revalue(pheno_table$biosource_short,
                       c("CD14"="Mono", "CD38"="Bcell", "CD4" = "CD4T",
                         "CD8"="CD8T", "cytotoxic CD56" = "NK",
                         "mature eosinophil" = "Gran",
                         "mature neutrophil" = "Gran")))
table(pheno_table$cell_type)

# load GRanges object
message("load GRanges object")
gr_complete <- readRDS(file = file.path(dataPath, "files_bsseq",
                                        "blueprint_blood_gr_02.RDS"))

# load Cov and M datasets 
message("load Cov and M datasets")
hdf5_bs_path <- file.path(dataPath, "files_bsseq",
                          "blueprint_blood_02.h5")
hdf5_cov  <- HDF5Array(filepath = hdf5_bs_path, name = "cov")
hdf5_meth <- HDF5Array(filepath = hdf5_bs_path, name = "meth")

pryr::object_size(gr_complete)
pryr::object_size(hdf5_cov)
pryr::object_size(hdf5_meth)

# Create BSseq object
message("Creating BSseq object")
Sys.time()
bs <- BSseq(gr = gr_complete, 
            M = hdf5_meth, 
            Cov = hdf5_cov, 
            sampleNames = pheno_table$sample_name)
Sys.time()

pData(bs) <- DataFrame(pheno_table)

message("Updating BSseq object")
Sys.time()
bs <- updateObject(bs) # maybe not necessary? 
Sys.time()

pryr::object_size(bs)

message("Saving BSseq object")
Sys.time()
hdf5_bs_se_path <- file.path(dataPath, "files_bsseq_hdf5_se")
saveHDF5SummarizedExperiment(bs, hdf5_bs_se_path, replace = TRUE,
                             chunkdim = c(2e3, ncol(bs)), verbose = TRUE)
Sys.time()
