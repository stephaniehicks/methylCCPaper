library(DeepBlueR)
library(foreach)
library(data.table)
library(readr)

workingDir_blueprint <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/blueprint-wgbs"
dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"

# Test installation and connectivity by saying hello to the DeepBlue server:
deepblue_info("me")

# ok this works
deepblue_list_genomes()
deepblue_list_projects() # "BLUEPRINT Epigenome", "DEEP (IHEC)"
deepblue_list_techniques() # "WGBS", "RRBS", "BisulfiteSeq"
deepblue_list_epigenetic_marks() # "DNA Methylation"
deepblue_list_biosources() # e.g. blood, muscle, etc
deepblue_list_experiments() 



# next we search experiments
keep_biosource <- c("CD14-positive, CD16-negative classical monocyte", 
                    "CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell", 
                    "CD38-negative naive B cell", 
                    "cytotoxic CD56-dim natural killer cell", 
                    "mature neutrophil", "mature eosinophil")

blueprint_DNA_meth <- deepblue_list_experiments(genome = "GRCh38",
                                                epigenetic_mark = "DNA Methylation",
                                                technique = "BisulfiteSeq",
                                                biosource = keep_biosource,
                                                project = "BLUEPRINT Epigenome")

# Then we remove `.bed` files and only lookg at `.wig` files
blueprint_DNA_meth <- 
  blueprint_DNA_meth %>% 
  filter(!grepl(".bed", name)) %>%
  data.table()

# To get more information about one experiment, use `deepblue_info()`
deepblue_info("e93346")


## Extract meta-data

# Using the experiment IDs, extract meta data about each sample, 
# including the biosource, etc. 
custom_table = do.call("rbind", apply(blueprint_DNA_meth, 1, function(experiment){
  experiment_id = experiment[1]

  # Obtain the information about the experiment_id
  info = deepblue_info(experiment_id)

  # Print the experiment name, project, biosource, and epigenetic mark.
  with(info, { data.frame(id = `_id`, name = name, project = project,
    technique = technique, epigenetic_mark = epigenetic_mark,
    biosource = sample_info$biosource_name,
    tissue_type = sample_info$TISSUE_TYPE,
    disease_status = extra_metadata$DISEASE,
    donor_id = sample_info$DONOR_ID,
    donor_age = extra_metadata$DONOR_AGE,
    donor_sex = extra_metadata$DONOR_SEX,
    experiment_id = extra_metadata$EXPERIMENT_ID,
    sample_id = sample_id,
    sample_name = sample_info$SAMPLE_NAME,
    ample_barcode = extra_metadata$SAMPLE_BARCODE,
    sample_description = extra_metadata$SAMPLE_DESCRIPTION,
    sample_source = sample_info$source,
    file_path = extra_metadata$FILE,
    first_submission_date = extra_metadata$FIRST_SUBMISSION_DATE,
    instrument_model = extra_metadata$INSTRUMENT_MODEL)
      })
}))
saveRDS(custom_table, file = file.path(dataPath,"blueprint_blood_custom_table.RDS"))

dim(custom_table)
head(custom_table)

# we also write a file with the paths to the bigwigs to download directly
write_csv(data.frame(paste0("ftp://ftp.ebi.ac.uk/pub/databases/", custom_table$file_path)), 
          file.path(dataPath,"blueprint_blood_ftp_paths.csv"), 
          col_names = FALSE)

# **note**After much effort, I failed to download the WGBS data 
# from the deepblueR bioconductor package or from the API directly. 
# Instead, I decided to use the blueprint_blood_ftp_paths.csv 
# file to download the data with wget. However, the code below this point 
# I wrote to try and download the data. But my requests kept failing.
# Maybe it will be useful for someone else. 

# Create two tables: one for `.call` files, and 
# one for `.cov` files. The  `.call` file contains
# the methylation signal (or percent of reads that
#                        are methylated). The `.cov` file contains the 
# coverage of the methylation signal (or how many
#                                    reads cover the CpG). 

table_call <- 
  custom_table %>% 
  filter(grepl("calls.bs_call", name)) %>% 
  data.table()

table_cov <- 
  custom_table %>% 
  filter(grepl("calls.bs_cov", name)) %>% 
  data.table()

###### Parallelizing it
# We can also split this up for all the chromosomes so we do not
# hit the download limit of DeepBlue. We also break up the 
# processing of the `.cov` and the `.call` files. 

# list all available chromosomes in GRCh38
chromosomes_GRCh38 <- deepblue_extract_ids(
  deepblue_chromosomes(genome = "GRCh38") )

# keep only the essential ones
chromosomes_GRCh38 <- 
  grep(pattern = "chr([0-9]{1,2}|X)$", chromosomes_GRCh38, 
       value = TRUE)

# We create `query_id`s, one for each chromosome to avoid 
# hitting the limits of deepblue. First we process the 
# `.call` files.

blueprint_regions_call <- 
  foreach(chr = chromosomes_GRCh38, .combine = c) %do% 
  {
    query_id = deepblue_select_experiments(
      experiment_name = 
        deepblue_extract_names(table_call), 
      chromosome = chr) 
  }
blueprint_regions_call # these are query_id's

# Then we process the `.cov` files. 
blueprint_regions_cov <- 
  foreach(chr = chromosomes_GRCh38, .combine = c) %do% 
  {
    query_id = deepblue_select_experiments(
      experiment_name = 
        deepblue_extract_names(table_cov), 
      chromosome = chr)
  }
blueprint_regions_cov # these are query_id's

# Next, we prepare to create the score matrix for the 
# `.call` and `.cov` files. 
exp_columns_call <- deepblue_select_column(table_call, "VALUE")
exp_columns_cov <- deepblue_select_column(table_cov, "VALUE")

# Then we submit requests for both the `.call` files 
request_ids_call <- foreach(query_id = blueprint_regions_call, 
                            .combine = c) %do%  
  {
    deepblue_score_matrix(
      experiments_columns = exp_columns_call,
      aggregation_function = "max",
      aggregation_regions_id = query_id)
  }
request_ids_call

# check to see if the requests are done 
foreach(request = request_ids_call, .combine = c) %do% {
  deepblue_info(request)$state
}

# And `.cov` files
request_ids_cov <- foreach(query_id = blueprint_regions_cov , 
                           .combine = c) %do%  
  {
    deepblue_score_matrix(
      experiments_columns = exp_columns_cov ,
      aggregation_function = "max",
      aggregation_regions_id = query_id)
  }
request_ids_cov 

# check to see if the requests are done 
foreach(request = request_ids_cov, .combine = c) %do% {
  deepblue_info(request)$state
}


# Once the requests are complete, we can create the score matrix.
list_score_matrices_call <-
  deepblue_batch_export_results(request_ids_call)

score_matrix_call <- data.table::rbindlist(
  list_score_matrices_call, use.names = TRUE)

score_matrix_call[, 1:5, with=FALSE]

saveRDS(score_matrix_call, 
        file = file.path(dataPath, "blueprint_blood_call.RDS"))
rm(score_matrix_call)



list_score_matrices_cov <-
  deepblue_batch_export_results(request_ids_cov)

score_matrix_cov <- data.table::rbindlist(
  list_score_matrices_cov, use.names = TRUE)

score_matrix_cov[, 1:5, with=FALSE]

saveRDS(score_matrix_cov, 
        file = file.path(dataPath, "blueprint_blood_cov.RDS"))
rm(score_matrix_cov)

