log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

adjust_transcriptome_ids <- function(tsv, transcriptome, output_path) {
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(vroom)
  })
  transcriptome <- import(transcriptome)
  disambiguated_proteome <- vroom::vroom(tsv)

  homo_transcriptome <- transcriptome
  a_transcriptome <- transcriptome
  b_transcriptome <- transcriptome
  a_transcriptome$transcript_id <- paste0(a_transcriptome$transcript_id, "A")
  b_transcriptome$transcript_id <- paste0(b_transcriptome$transcript_id, "B")

  haplotype_resolved_transcriptome <- c(transcriptome, a_transcriptome, b_transcriptome)
  export(
    haplotype_resolved_transcriptome[haplotype_resolved_transcriptome$transcript_id %in% disambiguated_proteome$unique_protein_id],
    output_path
  )
  return(0)
}

status <- adjust_transcriptome_ids(
  tsv = snakemake@input[["tsv"]],
  transcriptome = snakemake@input[["transcriptome"]],
  output_path = snakemake@output[[1]]
)

sessionInfo()

sink()
sink()
