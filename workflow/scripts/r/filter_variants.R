log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

filter_variants <- function(variants, transcriptome, min_gq, output_path) {
  variants <- readVcf(variants)
  transcriptome <- import(transcriptome)

  cds_transcriptome <- transcriptome[transcriptome$type == "CDS"]

  cds_variant_ids <- (unique(subjectHits(findOverlaps(cds_transcriptome, variants))))
  cds_variants <- variants[cds_variant_ids]

  gq_mask <- as.vector(geno(cds_variants)$GQ >= min_gq)
  cds_variants <- cds_variants[gq_mask]

  writeVcf(cds_variants, output_path)
  return(0)
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(VariantAnnotation)
})

status <- filter_variants(
  variants = snakemake@input[["variants"]],
  transcriptome = snakemake@input[["transcriptome"]],
  min_gq = snakemake@params[["min_gq"]],
  output_path = snakemake@output[["output_path"]]
)

sessionInfo()

sink()
sink()
