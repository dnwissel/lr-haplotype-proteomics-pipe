log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

filter_orfanage <- function(orfanage_gtf_file_path, output_path_to_be_predicted, output_path_filtered, minimum_orf_length) {
  transcriptome <- import(orfanage_gtf_file_path)
  transcriptome_non_cds <- transcriptome[transcriptome$type != "CDS"]
  transcriptome_cds <- transcriptome[transcriptome$type == "CDS"]
  above_threshold_cds <- data.frame(
    width = width(transcriptome_cds),
    transcript_id = transcriptome_cds$transcript_id
  ) %>%
    group_by(transcript_id) %>%
    summarise(total_cds_length = sum(width)) %>%
    filter(total_cds_length >= minimum_orf_length)
  passing_cds <- unique(above_threshold_cds$transcript_id)
  non_passing_cds <- setdiff(unique(transcriptome$transcript_id), passing_cds)
  export(transcriptome_non_cds[transcriptome_non_cds$transcript_id %in% non_passing_cds], output_path_to_be_predicted)
  export(transcriptome[transcriptome$transcript_id %in% passing_cds], output_path_filtered)
  return(0)
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
})

status <- filter_orfanage(
  orfanage_gtf_file_path = snakemake@input[["orfanage_gtf_file_path"]],
  output_path_to_be_predicted = snakemake@output[["output_path_to_be_predicted"]],
  output_path_filtered = snakemake@output[["output_path_filtered"]],
  minimum_orf_length = snakemake@params[["minimum_orf_length"]]
)

sessionInfo()

sink()
sink()
