log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

differentiate_per_isoform_statistics <- function(input_path, sqanti_protein, output_path) {
  suppressPackageStartupMessages({
    library(vroom)
    library(dplyr)
    library(rjson)
    library(readr)
  })
  variant_tracking <- vroom::vroom(input_path)
  sqanti_protein <- vroom::vroom(sqanti_protein)
  variant_tracking$id <- substr(variant_tracking$isoform_id, 1, nchar(variant_tracking$isoform_id) - 1)

  filtered_heterozygous <- variant_tracking[which(substr(variant_tracking$isoform_id, nchar(variant_tracking$isoform_id), nchar(variant_tracking$isoform_id)) %in% c("A", "B")), ]

  heterozygous_indels <- rep(0, nrow(filtered_heterozygous))
  homozygous_indels <- rep(0, nrow(filtered_heterozygous))

  heterozygous_snvs <- rep(0, nrow(filtered_heterozygous))
  homozygous_snvs <- rep(0, nrow(filtered_heterozygous))

  for (ix in 1:nrow(filtered_heterozygous)) {
    id <- filtered_heterozygous[ix, ]$id
    isoform_id <- filtered_heterozygous[ix, ]$isoform_id

    local_variants <- strsplit(filtered_heterozygous[ix, ]$contributing_variants, ";")[[1]]

    other_variants <- strsplit(filtered_heterozygous[filtered_heterozygous$id == id & filtered_heterozygous$isoform_id != isoform_id, ]$contributing_variants, ";")[[1]]
    if (all(is.na(other_variants)) & all(is.na(local_variants))) {
      next
    } else if (all(is.na(other_variants))) {
      local_indel_mask <- nchar(sapply(strsplit(local_variants, "_"), function(x) x[[2]])) > 3
      local_indels <- local_variants[local_indel_mask]
      local_snvs <- local_variants[!local_indel_mask]
      heterozygous_indels[ix] <- length(local_indels)
      homozygous_indels[ix] <- 0
      heterozygous_snvs[ix] <- length(local_snvs)
      homozygous_snvs[ix] <- 0
    } else if (all(is.na(local_variants))) {
      heterozygous_indels[ix] <- 0
      homozygous_indels[ix] <- 0
      heterozygous_snvs[ix] <- 0
      homozygous_snvs[ix] <- 0
    } else {
      local_indel_mask <- nchar(sapply(strsplit(local_variants, "_"), function(x) x[[2]])) > 3
      local_indels <- local_variants[local_indel_mask]
      local_snvs <- local_variants[!local_indel_mask]
      other_indel_mask <- nchar(sapply(strsplit(other_variants, "_"), function(x) x[[2]])) > 3
      other_indels <- other_variants[other_indel_mask]
      other_snvs <- other_variants[!other_indel_mask]
      heterozygous_indels[ix] <- length(local_indels[!(local_indels %in% other_indels)])
      homozygous_indels[ix] <- length(local_indels[(local_indels %in% other_indels)])
      heterozygous_snvs[ix] <- length(local_snvs[!(local_snvs %in% other_snvs)])
      homozygous_snvs[ix] <- length(local_snvs[(local_snvs %in% other_snvs)])
    }
  }

  heterozygous_variants <- heterozygous_snvs + heterozygous_indels
  homozygous_variants <- homozygous_snvs + homozygous_indels

  variant_tracking_non_heterozygous <- variant_tracking[-which(substr(variant_tracking$isoform_id, nchar(variant_tracking$isoform_id), nchar(variant_tracking$isoform_id)) %in% c("A", "B")), ]
  variant_tracking_non_heterozygous$n_heterozygous_variants <- 0
  variant_tracking_non_heterozygous$n_heterozygous_indels <- 0
  variant_tracking_non_heterozygous$n_heterozygous_snvs <- 0

  variant_tracking_non_heterozygous$n_homozygous_variants <- variant_tracking_non_heterozygous$n_variants
  variant_tracking_non_heterozygous$n_homozygous_indels <- variant_tracking_non_heterozygous$n_indels
  variant_tracking_non_heterozygous$n_homozygous_snvs <- variant_tracking_non_heterozygous$n_snvs

  filtered_heterozygous$n_heterozygous_variants <- heterozygous_variants
  filtered_heterozygous$n_heterozygous_indels <- heterozygous_indels
  filtered_heterozygous$n_heterozygous_snvs <- heterozygous_snvs

  filtered_heterozygous$n_homozygous_variants <- homozygous_variants
  filtered_heterozygous$n_homozygous_indels <- homozygous_indels
  filtered_heterozygous$n_homozygous_snvs <- homozygous_snvs

  rbind(variant_tracking_non_heterozygous, filtered_heterozygous) %>%
    select(
      isoform_id,
      protein_splice_category,
      n_variants,
      n_snvs,
      n_indels,
      n_heterozygous_variants,
      n_heterozygous_snvs,
      n_heterozygous_indels,
      n_homozygous_variants,
      n_homozygous_indels,
      n_homozygous_snvs,
      contributing_variants,
      diff,
      has_indel,
      has_changed_stop,
      has_frameshift,
      has_resolved_frameshift
    ) %>%
    filter(isoform_id %in% sqanti_protein$pb) %>%
    write_tsv(output_path)
  return(0)
}

status <- differentiate_per_isoform_statistics(
  input_path = snakemake@input[["input_path"]],
  sqanti_protein = snakemake@input[["sqanti_protein"]],
  output_path = snakemake@output[[1]]
)

sessionInfo()

sink()
sink()
