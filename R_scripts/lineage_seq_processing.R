# for taking the .stats output from SingleCellLineage pipeline, extracting UMIs (ftag, clonetag, etc), subsetting for valid 10x cell IDs

##### packages #####
library(ggplot2)
library(pals)
library(cowplot)
library(ggpubr)
library(devtools)
library(ShortRead) 
library(data.table)
library(dplyr)
library(stringr)
library(BiocManager)
library(viridis)
library(radiant.data)
library(RColorBrewer)
library(bayesbio)
library(grid)

##### processing output from SingleCellLineage pipeline #####
recorder_stats_file <- fread("/dartfs-hpc/rc/lab/M/McKennaLab/projects/saxe/data/illumina/240123_RS21_RS22_initial/data/pipeline_output_fixed_castag_ref/rs21_rs22_lin/rs21_rs22_lin.stats")

# load cell ID whitelist associated with 10X chemistry
rosetta <- read.table("/dartfs/rc/lab/M/McKennaLab/resources/cellranger_versions/cellranger-7.2.0/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz") #pre GEM-X 10x whitelist
rosetta <- read.table("/oldhome/rachel/updated_cellranger/cellranger-8.0.1/lib/python/cellranger/barcodes/translation/3M-3pgex-may-2023.txt.gz") #GEM-X 10x whitelist



extract_static_barcode = function(read,ref) {
  ref_n_locations = str_locate_all(ref, "N")
  n_starts = ref_n_locations[[1]][1,]["start"]
  n_stops = ref_n_locations[[1]][nrow(ref_n_locations[[1]]),]["end"]
  return(substr(read,n_starts,n_stops))
}

extract_merged_or_paired_barcode = function(x) {
  if (x["merged"] == "MERGED") {
    return(extract_static_barcode(x["mergedRead"],x["mergedReadRef"]))
  } else if (x["merged"] == "PAIR") {
    return(extract_static_barcode(x["fwdRead"],x["fwdReadRef"]))
  } else {
    return("UNKNOWN")
  }
}


get_UMI_cellID_ftag_10xdata = function(df) {
  df = df[df$keep == "PASS",]
  ftag <- apply(df,1,extract_merged_or_paired_barcode)
  umi_cellID <- sapply(strsplit(df$readName, '_'), `[`, 2)
  df <- df %>%
    mutate(ftag, umi_cellID)
  cell_ID <- substr(df$umi_cellID,1,16)
  umi_10x <- substr(df$umi_cellID,17,28)
  df <- df %>%
    mutate(umi_10x, cell_ID)
  #df <- as.data.frame(df[,-c(2:18, 20:22, 40)])
}

get_UMI_cellID_castag_10xdata = function(df) {
  df = df[df$keep == "PASS",]
  cas_tag <- apply(df,1,extract_merged_or_paired_barcode)
  umi_cellID <- sapply(strsplit(df$readName, '_'), `[`, 2)
  df <- df %>%
    mutate(cas_tag, umi_cellID, umi)
  cell_ID <- substr(df$umi_cellID,1,16)
  umi_10x <- substr(df$umi_cellID,17,28)
  df <- df %>%
    mutate(umi_10x, cell_ID)
  df <- df[-grep("-{2,14}",df$cas_tag),]
  #df <- df[,c(1, 19, 25, 27, 28)]
}

get_UMI_non10x = function(df) {
  df = df[df$keep == "PASS",]
  ftag <- apply(df,1,extract_merged_or_paired_barcode)
  df <- df %>%
    mutate(ftag)
  #df <- as.data.frame(df[,-c(2:18, 20:21)])
}


bulk_ftag_clonetag_reads <- get_UMI_non10x(bulk_stats_file)
recorder_ftag_reads_10x <- get_UMI_cellID_ftag_10xdata(recorder_stats_file)
clonetag_reads_10x <- get_UMI_cellID_castag_10xdata(clonetag_stats_file)

# join with whitelist to get pA capture/GEX cell IDs
recorder_ftag_reads_10x <- recorder_ftag_reads_10x %>% inner_join(rosetta, 
                                    by=c('cell_ID'='V2'))
clonetag_reads_10x <- clonetag_reads_10x %>% inner_join(rosetta, 
                                    by=c('cell_ID'='V2'))


