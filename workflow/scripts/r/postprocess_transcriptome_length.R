log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

postprocess_transcriptome_length <- function(transcriptome, unfiltered_transcriptome, output_path, minimum_orf_length) {
  unfiltered_transcriptome <- import(unfiltered_transcriptome)
  unfiltered_transcriptome <- unfiltered_transcriptome[
    seqnames(unfiltered_transcriptome) %in% paste0("chr", 1:22) 
    & unfiltered_transcriptome$gene_type == "protein_coding"
  ]
  transcriptome <- import(transcriptome)
  transcriptome <- transcriptome[transcriptome$gene_id %in% unfiltered_transcriptome$gene_id]
  transcriptome_cds <- transcriptome[transcriptome$type == "CDS"]
  above_threshold_cds <- data.frame(
    width = width(transcriptome_cds),
    transcript_id = transcriptome_cds$transcript_id
  ) %>%
    group_by(transcript_id) %>%
    summarise(total_cds_length = sum(width)) %>%
    filter(total_cds_length >= minimum_orf_length)

  filtered_transcriptome <- transcriptome[transcriptome$transcript_id %in% unique(above_threshold_cds$transcript_id)]
  export(filtered_transcriptome, output_path)
  return(0)
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
})

status <- postprocess_transcriptome_length(
  transcriptome = snakemake@input[["transcriptome"]],
  unfiltered_transcriptome = snakemake@input[["unfiltered"]],
  output_path = snakemake@output[["output_path"]],
  minimum_orf_length = snakemake@params[["minimum_orf_length"]]
)

sessionInfo()

sink()
sink()
