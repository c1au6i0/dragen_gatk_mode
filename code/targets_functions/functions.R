ggplot2::theme_set(ggplot2::theme_bw())

#' download references
#'
#' Download references files from gs://gcp-public-data--broad-references/hg38/v0/ and track files.
#'
#' @param path_download Where to save data.
download_references <- function(path_download) {
  gsutil_installed <- system("gsutil --version") == 0


  if (!dir.exists(path_download)){
    stop("path_download dir does not exist!")
  }

  if (!gsutil_installed) {
    stop("gsutil installation not found! Please install macs2...")
  }


  system(
         paste0(
               "gsutil -m cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.str ",
              path_download
              )
  )
  system(
         paste0(
              "gsutil -m rsync  gs://gcp-public-data--broad-references/hg38/v0/dragen_reference/ ",
              path_download
              )
  )



  # path are saved in a v0 folder

  list.files(path = path_download, include.dirs = FALSE, full.names = TRUE, recursive = TRUE)

}

#' align
#'
#' Use `dragen-os` to align bam to reference.
#'
#' @param reference_path Path to reference folder.
#' @param ubam Path to ubam.
#' @param output_dir Path to output folder.
#' @param output_prefix String of prefix name.
#' @param nun_threads Max 56 for alignment and `samtools view`
#' @param interleaved 0 for SE, 1 for PE.
dragen_align <- function(reference_path, ubam, output_dir, output_prefix, num_threads, interleaved = 1){

     bam_out_path <- c(paste0(file.path(output_dir, output_prefix), ".bam"))

     align_bash <- paste0(
       "dragen-os -r ", reference_path,
       " -b ", ubam,
       " --num-threads ", num_threads,
       " --interleaved ", interleaved,
       " | samtools view -h -@ ", num_threads, " -O BAM - > ", bam_out_path
     )
     system(align_bash)
     bam_out_path
}



#' sort and index
#'
#' Use `samtools` to sort and index bams
#'
#' @param tmp_path for tmp_files.
#' @param bam Path to bam.
#' @param output_dir Path to output folder.
#' @param num_threads
sort_index <- function(bam, tmp_path, output_dir, num_threads){


     sorted_bam_path <- file.path(output_dir, paste0("sorted_", basename(bam)))
     sorted_bai_path <- paste0(sorted_bam_path, ".bai")

     message("Sorting ...")
     sort_bash <- paste0(
       "samtools sort -T ", tmp_path,
       " ", bam,
       " -@ ", num_threads,
       " -o ", sorted_bam_path
     )

     system(sort_bash)

     message("Indexing ...")
    index_bash <- paste0(
       "samtools index -@ ", num_threads,
       " ", sorted_bam_path
     )

     system(index_bash)

    c(sorted_bam_path, sorted_bai_path)
}

#' collate fixmate sort
#'
#' Use `samtools` to collating, fixmate and sorting bams
#'
#' @param bam Path to bam.
#' @param output_dir Path to output folder.
#' @param num_threads
#' @param tmp_path for tmp_files.
sort_fixmate_collate_bam :w<- function(bam, output_dir, num_threads, tmp_path){

    bam_paths <- sapply(c(collated = "collated", fixmated = "fixmated", sorted = "sorted"),
           function(x) file.path(output_dir, paste0(x, "_", basename(bam))),
           simplify = FALSE)

     message(paste0("Collating, Fixmating and Sorting ", basename(bam)))
     command_bash <- paste0(
       "set -e;",
       "echo 'Collating ...'; ",
       "samtools collate ", bam, " -o ", bam_paths$collated,
       " -@ ", num_threads, " ; ",
       "echo 'Fixmating...' ; ",
       "samtools fixmate -m -@ ",  num_threads, " ",  bam_paths$collated, " ", bam_paths$fixmated, " ; ",
       "echo 'Sorting...' ; ",
       "samtools sort ", bam_paths$fixmated, " -T ", tmp_path, " ", " -@ ", num_threads, " -o ", bam_paths$sorted
     )

     system(command_bash)
     bam_paths$sorted
}



#' remap id and find dup
#'
#' use `samtools` to `markdup`
#'
#' @param sorted_bam path to bam.
#' @param ppcg_mapping Dataframe with 2 columns, first with PPCG_id and second with PCAWG_ids.
#' @param output_dir Path to output folder.
#' @param num_threads Number of threads.
#' @param tmp_path For tmp_files.
remapid_mark_duplex <- function(sorted_bam, ppcg_mapping = ppcg_mapping, output_dir, num_threads, tmp_path){
     pcawg_name <- basename(sorted_bam)
     ppcg_name <- as.character(ppcg_mapping[ppcg_mapping[, 2] == pcawg_name, 1])

     ppcg_file_path <- file.path(output_dir, paste0(ppcg_name, ".bam"))

     message("Marking duplicate allignment and changing name of file...")
     command_bash <- paste0(
       "samtools markdup ", sorted_bam, " ", ppcg_file_path,
       " -@ ", num_threads,
       " -T ", file.path(tmp_path,  paste0(ppcg_name, ".bam"))
     )

     system(command_bash)
     ppcg_file_path
}





#' download_contamination_resources
#'
#' @param path_out
download_contamination_resources <- function(path_out = here::here("data", "big_files", "1000_genomes")){

  dir.create(path_out, showWarnings = FALSE)

  download.file(
                url = "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.UD",
                destfile = file.path(path_out, "1000g.phase3.100k.b38.vcf.gz.dat.UD"),
                method = "wget"
  )

  download.file(
                url = "https://storage.cloud.google.com/gcp-public-data--broad-references/hg38/v0/contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.bed",
                destfile = file.path(path_out, "1000g.phase3.100k.b38.vcf.gz.dat.bed"),
                method = "wget"
  )

  download.file(
                url = "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.mu",
                destfile = file.path(path_out, "1000g.phase3.100k.b38.vcf.gz.dat.mu"),
                method = "wget"
  )
  download.file(
                url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
                destfile = file.path(path_out, "GRCh38_full_analysis_set_plus_decoy_hla.fa"),
                method = "wget"
  )

  list.files(path = path_out, include.dirs = FALSE, full.names = TRUE, recursive = TRUE)

}


# library(here)
# library(targets)
# tar_load(bams)
#
# sort_fixmate_collate_bam(
#                bam = bams[1],
#                tmp_path = here("tmp"),
#                output_dir = here("data", "big_files", "sorted_bam_files"),
#                num_threads = 50
    # )


