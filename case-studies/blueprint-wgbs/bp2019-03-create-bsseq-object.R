library(rtracklayer)
library(bsseq)
library(liftOver)

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

# qrsh -l mem_free=90G,h_vmem=100G
workingDir_blueprint <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs"
dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"

custom_table <- readRDS(file = file.path(dataPath,"blueprint_blood_custom_table.RDS"))
custom_table$file_name <- unlist(lapply(stringr::str_split(custom_table$file_path, "CNAG/"), 
                                        function(x){ x[[2]]}))
custom_table$biosource_short <- unlist(lapply(stringr::str_split(custom_table$biosource, "-"), 
                                              function(x){ x[[1]]}))

for(smps in levels(custom_table$sample_name))
{
  custom_table_call <- custom_table[custom_table$sample_name %in% smps &
                                      grepl("_calls.bs_call", custom_table$file_name), ]
  custom_table_cov  <- custom_table[custom_table$sample_name %in% smps &
                                      grepl("_calls.bs_cov", custom_table$file_name), ]
  
  gr_call <- import(file.path(dataPath, "files_bigwig", 
                              custom_table_call$file_name), 
                    format = "bigWig")
  colnames(mcols(gr_call)) <-  "score_call" 
  
  gr_cov <- import(file.path(dataPath, "files_bigwig", 
                             custom_table_cov$file_name), 
                   format = "bigWig")
  colnames(mcols(gr_cov)) <-  "score_cov" 
  
  sh <- findOverlaps(gr_call, gr_cov)
  gr_call <- gr_call[queryHits(sh),]
  gr_cov <- gr_cov[subjectHits(sh),]
  gr_call$score_cov <- gr_cov$score_cov
  genome(gr_call) <- "GRCh38"
  gr_call_19 <- liftOver(gr_call, ch)
  gr_call_19_rbind <- unlist(gr_call_19)
  genome(gr_call_19_rbind) <- "hg19"
  gr_call_19_rbind <- unique(gr_call_19_rbind)
  rm(gr_call, gr_cov, gr_call_19)
  
  # bsseq objects
  methyl <- round(as.matrix(mcols(gr_call_19_rbind)$score_call) * 
                    as.matrix(mcols(gr_call_19_rbind)$score_cov))
  colnames(methyl) <- as.character(custom_table_cov$sample_name)
  
  BS <- BSseq(gr = gr_call_19_rbind,
              M = methyl,
              Cov = as.matrix(mcols(gr_call_19_rbind)$score_cov))
  custom_table_call$file_path_cov <- custom_table_cov$file_path
  custom_table_call$file_name_cov <- custom_table_cov$file_name
  pData(BS) <- DataFrame(custom_table_call)
  BS <- updateObject(BS)
  
  saveRDS(BS, file = file.path(dataPath, "files_bsseq",
        paste0("blueprint_blood_bs_", 
               custom_table_cov$sample_name, ".RDS")))
}

