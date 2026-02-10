log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

create_protein_db <- function(splice_proteome, haplosaurus_json, output_path) {
  suppressPackageStartupMessages({
    library(jsonlite)
    library(readr)
    library(Biostrings)
  })
  splice_proteome_sequences <- Biostrings::readAAStringSet(splice_proteome)
  ### Approach from SO: https://stackoverflow.com/questions/60298776/convert-json-file-with-multiple-lines-to-r-dataframe
  json_raw <- readr::read_file(haplosaurus_json)
  json_lines <- unlist(strsplit(json_raw, "\\n"))
  json_parsed <- lapply(json_lines, jsonlite::fromJSON)
  mask <- !(names(splice_proteome_sequences) %in% sapply(json_parsed, function(x) x$transcript_id))

  haplosaurus_sequences <- unlist(lapply(json_parsed, function(x) substr(x$protein_haplotypes$seq, 1, nchar(x$protein_haplotypes$seq) - 1)))

  haplosaurus_names <- unlist(lapply(json_parsed, function(x) {
    if (length(x$protein_haplotypes$seq) < 2) {
      return(x$transcript_id)
    } else {
      if (length(x$protein_haplotypes$contributing_variants[[1]]) >= length(x$protein_haplotypes$contributing_variants[[2]])) {
        return(
          c(
            paste0(x$transcript_id, "A"),
            paste0(x$transcript_id, "B")
          )
        )
      }
      else {
        return(
          c(
            paste0(x$transcript_id, "B"),
            paste0(x$transcript_id, "A")
          )
        )
      }

    }
  }))
  finalized_proteome <- c(
    splice_proteome_sequences[mask],
    setNames(AAStringSet(haplosaurus_sequences), haplosaurus_names)
  )
  writeXStringSet(finalized_proteome, output_path)
  return(0)
}

status <- create_protein_db(
  splice_proteome = snakemake@input[["splice_proteome"]],
  haplosaurus_json = snakemake@input[["haplosaurus_json"]],
  output_path = snakemake@output[[1]]
)

sessionInfo()

sink()
sink()
