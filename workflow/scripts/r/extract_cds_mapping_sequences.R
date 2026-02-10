log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

extract_cds_mapping_sequences <- function(sqanti_protein, haplosaurus_json, splice_fasta, genome, variants, transcriptome, output_path_homozygous,
                                          output_path_first_haplotype, output_path_second_haplotype) {
  suppressPackageStartupMessages({
    library(vroom)
    library(dplyr)
    library(rjson)
    library(readr)
    library(VariantAnnotation)
    library(transmogR)
  })

  json_raw <- readr::read_file(haplosaurus_json)
  json_lines <- unlist(strsplit(json_raw, "\\n"))
  json_parsed <- lapply(json_lines, jsonlite::fromJSON)

  splice_fasta <- readDNAStringSet(splice_fasta)
  protein_class <- vroom::vroom(sqanti_protein)
  relevant_ids <- unique(protein_class$pb)
  var <- readVcf(variants)
  genome <- readDNAStringSet(genome)
  transcriptome <- import(transcriptome)

  names(genome) <- sapply(strsplit(names(genome), "\\ "), function(x) x[[1]])
  splice_fasta <- splice_fasta[names(splice_fasta) %in% relevant_ids]

  isoform_id <- unlist(lapply(json_parsed, function(x) {
    if (length(x$protein_haplotypes$seq) < 2) {
      return(x$transcript_id)
    } else {
      return(
        c(
          paste0(x$transcript_id, "A"),
          paste0(x$transcript_id, "B")
        )
      )
    }
  }))

  cds_mapping <- rep("", length(isoform_id))

  ix <- 1
  for (transcript in 1:length(json_parsed)) {
    id <- json_parsed[[transcript]]$transcript_id
    if (length(json_parsed[[transcript]]$protein_haplotypes$seq) > 1) {
      cds_mapping[ix] <- as.character(transmogrify(
        x = genome,
        var = rowRanges(var[which(names(rowRanges(var)) %in% json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]]), ]),
        exons = transcriptome[which(transcriptome$transcript_id == id & transcriptome$type == "CDS"), ],
        verbose = FALSE
      ))
      ix <- ix + 1
      cds_mapping[ix] <- as.character(transmogrify(
        x = genome,
        var = rowRanges(var[which(names(rowRanges(var)) %in% json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[2]]), ]),
        exons = transcriptome[which(transcriptome$transcript_id == id & transcriptome$type == "CDS"), ],
        verbose = FALSE
      ))
      ix <- ix + 1
    } else {
      cds_mapping[ix] <- as.character(transmogrify(
        x = genome,
        var = rowRanges(var[which(names(rowRanges(var)) %in% json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]]), ]),
        exons = transcriptome[which(transcriptome$transcript_id == id & transcriptome$type == "CDS"), ],
        verbose = FALSE
      ))
      ix <- ix + 1
    }
  }

  cds_fasta <- DNAStringSet(cds_mapping)
  names(cds_fasta) <- isoform_id

  cds_fasta <- c(
    cds_fasta,
    splice_fasta[!names(splice_fasta) %in% isoform_id]
  )
  cds_fasta <- cds_fasta[names(cds_fasta) %in% relevant_ids]
  last_letter_cds_name <- substr(names(cds_fasta), nchar(names(cds_fasta)), nchar(names(cds_fasta)))
  writeXStringSet(cds_fasta[!(last_letter_cds_name %in% c("A", "B"))], output_path_homozygous)
  writeXStringSet(cds_fasta[(last_letter_cds_name %in% c("A"))], output_path_first_haplotype)
  writeXStringSet(cds_fasta[(last_letter_cds_name %in% c("B"))], output_path_second_haplotype)

  return(0)
}

status <- extract_cds_mapping_sequences(
  sqanti_protein = snakemake@input[["sqanti_protein"]],
  haplosaurus_json = snakemake@input[["haplosaurus_json"]],
  splice_fasta = snakemake@input[["splice_fasta"]],
  genome = snakemake@input[["genome"]],
  variants = snakemake@input[["variants"]],
  transcriptome = snakemake@input[["transcriptome"]],
  output_path_homozygous = snakemake@output[["output_path_homozygous"]],
  output_path_first_haplotype = snakemake@output[["output_path_first_haplotype"]],
  output_path_second_haplotype = snakemake@output[["output_path_second_haplotype"]]
)

sessionInfo()

sink()
sink()
