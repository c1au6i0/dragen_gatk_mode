library(targets)
library(future)
library(future.batchtools)
library(batchtools)
library(tarchetypes)
library(here)
library(ggplot2)

# Initiate functions
source(here("code", "targets_functions", "functions.R"))

ggplot2::theme_set(theme_bw())

options(tidyverse.quiet = TRUE)


# Initiate file for defining slurm resources to use in targets
source("_targets_resources.conf.R")

tar_option_set(
  resources = resources_all,
  packages = c(
    "here",
    "tidyverse"
  ),
  workspace_on_error = TRUE,
  storage = "worker",
  retrieval = "worker",
  # debug = "qc_res"
  garbage_collection = TRUE,
  memory = "transient"
)


list(

  # @@@@@@@@@@@@
  # Set up -----
  # @@@@@@@@@@@@

  # Download reference files from gs://gcp-public-data--broad-references/hg38/v0/

  tar_target(
    reference_folder,
    "/athena/marchionnilab/scratch/clz4002/ppcg_wgs/dragenmode/data/big_files/references"
  ),
  # This is the file that allows to map PCAWG to PPCG id and it is used later on to map file in the mark duplex stage
  tarchetypes::tar_file_read(
    ppcg_mapping,
    here::here("data", "id_mapping.csv"),
    read_csv(file = !!.x, col_types = cols())
  ),

  tar_target(
    ubam_folder,
    "/athena/marchionnilab/scratch/lab_data/sharedData/ppcg/ppcg_wgs"
    ),
  tar_target(
    reference_files,
    download_references(path_download = reference_folder),
    format = "file"
  ),
  tarchetypes::tar_files(
    ubams,
    list.files(path = ubam_folder, pattern = "unmapped", recursive = FALSE, include.dirs = FALSE, full.names = TRUE)
  ),

  # @@@@@@@@@@@@@@@@
  # Allignment -----
  # @@@@@@@@@@@@@@@@

  tar_target(
    bams,
    dragen_align(reference_folder,
                 ubams,
                 output_dir = here("data", "big_files", "bam_files"),
                 output_prefix = gsub("_unmapped\\.bam", "", basename(ubams)),
                 num_threads = 50,
                 interleaved = 1
                 ),
    resources = resources_alignment,
    map(ubams),
    format = "file"
  ),

 tar_target(
    sorted_bams,
    sort_fixmate_collate_bam(
               bam = bams,
               tmp_path = here("tmp"),
               output_dir = here("data", "big_files", "sorted_bam_files"),
               num_threads = 50
    ),
    resources = resources_alignment,
    map(bams),
    format = "file"
  ),

  # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Duplication and Contamination -----
  # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#   tar_target(
    # dup_bams,
    # remapid_mark_duplex(
    #            sorted_bam = sorted_bams,
    #            ppcg_mapping = ppcg_mapping,
    #            tmp_path = here::here("tmp"),
    #            output_dir = here("data", "big_files", "dup_bam_files"),
    #            num_threads = 50
    # ),
    # resources = resources_alignment,
    # map(sorted_bams),
    # format = "file"
  # ),

  tar_target(
    contamination_references,
    download_contamination_resources(here::here("data", "big_files", "1000_genomes")),
    format = "file"
  )
)
