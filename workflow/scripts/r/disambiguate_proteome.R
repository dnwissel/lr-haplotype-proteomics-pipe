log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

disambiguate_proteome_gencode <- function(input_path_fasta, output_path_fasta, output_path_tsv) {
  proteome <- readAAStringSet(input_path_fasta)
  proteome <- sort(proteome[nchar(proteome) > 34])

  unique_proteome <- unique(proteome)

  unique_mapping <- sfsmisc::uniqueL(as.character(proteome))

  disambiguation_frame <- data.frame(
    unique_protein_id = names(unique_proteome),
    all_contained_protein_ids = unlist(lapply(sapply(1:length(unique_proteome), function(x) names(proteome)[which(unique_mapping[[1]] == x)]), function(y) paste0(y, collapse = "|")))
  )
  mask <- which(grepl("Bambu", disambiguation_frame$unique_protein_id) & grepl("ENST", disambiguation_frame$all_contained_protein_ids))


  mask_names <- unlist(lapply(strsplit(disambiguation_frame[mask, ]$all_contained_protein_ids, "\\|"), function(x) {
    first(grep("ENST", x, value = TRUE))
  }))

  disambiguation_frame$unique_protein_id[mask] <- mask_names
  names(unique_proteome)[mask] <- mask_names

  readr::write_tsv(disambiguation_frame, output_path_tsv)
  export(unique_proteome, output_path_fasta)

  return(0)
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(readr)
  library(sfsmisc)
  library(dplyr)
})

status <- disambiguate_proteome_gencode(
  input_path_fasta = snakemake@input[["input_path_fasta"]],
  output_path_fasta = snakemake@output[["output_path_fasta"]],
  output_path_tsv = snakemake@output[["output_path_tsv"]]
)

sessionInfo()

sink()
sink()
