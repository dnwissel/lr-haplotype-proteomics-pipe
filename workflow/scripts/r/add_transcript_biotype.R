log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

add_transcript_biotype <- function(input_path, output_path) {
  transcriptome <- import(input_path)
  transcriptome$biotype <- "protein_coding"
  export(transcriptome, output_path, "gff3")

  return(0)
}

suppressPackageStartupMessages({
  library(rtracklayer)
})

status <- add_transcript_biotype(
  input_path = snakemake@input[[1]],
  output_path = snakemake@output[[1]]
)

sessionInfo()

sink()
sink()
