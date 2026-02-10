log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

postprocess_phasing <- function(vcf, query_transcriptome, output_path) {
  vcf <- VariantAnnotation::readVcf(vcf)
  bambu <- import(query_transcriptome)

  bambu_transcripts <- (unique(bambu$transcript_id))
  variant_state <- rep(1, nrow(vcf))
  for (id in 1:length(bambu_transcripts)) {
    transcript_id <- bambu_transcripts[id]
    matches <- GenomicRanges::findOverlaps(vcf, bambu[bambu$transcript_id == transcript_id]) #
    if (length(matches) > 0) {
      unique_query_matches <- unique(queryHits(matches))
      genotypes <- geno(vcf[unique_query_matches])$GT
      phase_sets <- geno(vcf[unique_query_matches])$PS

      hetero_mask <- as.vector(!(genotypes == "1/1"))

      unique_hetero_query_matches <- c(unique_query_matches)[hetero_mask]
      hetero_genotypes <- genotypes[hetero_mask]
      hetero_phase_sets <- phase_sets[hetero_mask]

      if (sum(hetero_mask) == 0) {
        next
      } else if (sum(hetero_mask) == 1) {
        if (!grepl("\\|", hetero_genotypes)) {
          variant_state[unique_hetero_query_matches] <- min(variant_state[unique_hetero_query_matches], 0)
        } else {
          next
        }
      } else if (sum(hetero_mask) > 1) {
        if (all(grepl("\\|", hetero_genotypes)) & length(unique(hetero_phase_sets)) == 1 & sum(is.na(hetero_phase_sets)) == 0) {
          next
        } else {
          variant_state[unique_hetero_query_matches] <- -1
        }
      }
    }
  }

  geno(vcf)$GT[variant_state == 0, 1] <- "0|1"
  mask <- variant_state %in% c(0, 1)
  VariantAnnotation::writeVcf(vcf[mask], output_path)
  return(0)
}

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(rtracklayer)
})

status <- postprocess_phasing(
  vcf = snakemake@input[["variants"]],
  query_transcriptome = snakemake@input[["query_transcriptome"]],
  output_path = snakemake@output[["resolved_vcf"]]
)

sessionInfo()

sink()
sink()
