log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

get_directly_observed_and_linked_variants <- function(variants,
                                                      inferred_proteins,
                                                      variant_annotation,
                                                      linked_variants_path,
                                                      observed_variants_path) {
  variant_evidence <- readVcf(variants)
  observed_proteins <- vroom::vroom(inferred_proteins, delim = ",")
  variant_annotation <- vroom::vroom(variant_annotation)

  mask <- unlist(info(variant_evidence)$homo_ac > 0) | unlist(info(variant_evidence)$first_ac > 0) | unlist(info(variant_evidence)$second_ac > 0)

  observed_variants <- variant_evidence[mask]

  observed_variants



  all_associated_variants <- sapply(
    observed_proteins$protein_group,
    function(x) {
      Reduce(intersect, as.list(variant_annotation[variant_annotation$isoform_id %in% unlist(strsplit(x, "-")), ]$contributing_variants))
    }
  )

  all_associated_variants <- unique(unlist(all_associated_variants))
  all_associated_variants <- all_associated_variants[!is.na(all_associated_variants)]
  all_associated_variants <- unlist(strsplit(all_associated_variants, ";"))
  linked_variants <- all_associated_variants[!all_associated_variants %in% names(rowRanges(observed_variants))]

  mask <- names(rowRanges(variant_evidence)) %in% linked_variants
  linked_variants <- variant_evidence[mask]

  writeVcf(linked_variants, linked_variants_path)
  writeVcf(observed_variants, observed_variants_path)
  return(0)
}

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(readr)
  library(VariantAnnotation)
})

status <- get_directly_observed_and_linked_variants(
  variants = snakemake@input[["variants"]],
  inferred_proteins = snakemake@input[["inferred_proteins"]],
  variant_annotation = snakemake@input[["variant_annotation"]],
  linked_variants_path = snakemake@output[["linked_variants_path"]],
  observed_variants_path = snakemake@output[["observed_variants_path"]]
)

sessionInfo()

sink()
sink()
