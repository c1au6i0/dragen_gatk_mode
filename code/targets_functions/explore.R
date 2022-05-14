# This is used as file for experiment with code that will became   functions

library(ChIPQC)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(targets)
library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(ggupset)
library(here)


sorted_bam <- tar_read(sorted_bams_e5a093da)
atac_reads_open <- tar_read(atac_reads_open_4d336326)

tar_load(atac_reads_open_d22787ac)


export_bigWig(sorted_bam, atac_reads_open)

str(atac_reads_open_d22787ac)


export_open_bam2 <- function(sorted_bam){

  open_region_bam_path <- gsub("\\.bam", "_open_regions\\.bam", sorted_bam)
  sys_command <- paste0(
                        "samtools view -h ",
                        sorted_bam,
                        " | ",
                        "awk 'substr($0,1,1)==\"@\" || ($9 >= 0 && $9 <= 100) ||  ($9 <= 0 && $9 >= -100)' | samtools view -b > ",
                        open_region_bam_path
                        )

  system(sys_command)
  # return path for targets tracking
  open_region_bam_path

}

tar_load(sorted_bams)

 call_peaks("/athena/marchionnilab/scratch/clz4002/inghirami/atac_seq/210318_A00814_0381_AH5HYLDRXY/inghirami_atacseq/data/big_files/sorted_bam/sorted_BELLI_DMSO_1_open_regions.bam",  here("data", "big_files", "peak_files"))


   qc_res <- ChIPQC::ChIPQCsample(
                 open_region_bams,
                 peaks = grep("narrowPeak$", peaks_files, value = TRUE), # take path of file ending in narrowPeak
                 annotation = "hg38",
                 chromosomes = chr_check,
                 blacklist = black_list_hg38,
                 verboseT = FALSE)




export_bw(sorted_bama[1], atac_reads_open[1])


export_bw <- function(sorted_bam, atac_reads_open, folder_bw = "bigWigs"){
   
  browser()
  if (!dir.exists(here::here("data", "big_files", folder_bw))) stop("There is not a folder with that name in data/big_files!")

  open_region_bw_path_pre <- gsub("\\.bam", "_open_regions\\.bw", sorted_bam)
  open_region_bw_path <- gsub("sorted_bam", folder_bw, open_region_bw_path_pre)


  atac_fragments_opens <- GenomicRanges::granges(atac_reads_open)
  rtracklayer::export.bw(GenomicRanges::coverage(atac_fragments_opens), open_region_bw_path)

  # return path for targets tracking
  open_region_bw_path

}
sorted_bams[1]
