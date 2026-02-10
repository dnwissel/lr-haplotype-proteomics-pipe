log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

unify_references <- function(vcf, peptide_bam_homo, peptide_bam_first, peptide_bam_second) {
  suppressPackageStartupMessages({
    library(Biostrings)
    library(rtracklayer)
    library(Rsamtools)
    library(VariantAnnotation)
  })
  genome <- readDNAStringSet(genome)
  transcriptome <- import(transcriptome)
  new_levels <- levels(seqnames(transcriptome))
  for (contig in unique(seqnames(transcriptome))[!(unique(seqnames(transcriptome)) %in% c(paste0("chr", c(1:22, "X", "Y", "M"))))]) {
    new_levels <- c(new_levels, grep(contig, sapply(strsplit(names(genome), " "), function(x) x[[1]]), value = TRUE))
  }
  seqlevels(transcriptome) <- new_levels
  for (contig in unique(seqnames(transcriptome))[!(unique(seqnames(transcriptome)) %in% c(paste0("chr", c(1:22, "X", "Y", "M"))))]) {
    seqnames(transcriptome[seqnames(transcriptome) == contig]) <- factor(grep(contig, sapply(strsplit(names(genome), " "), function(x) x[[1]]), value = TRUE), levels = seqlevels(transcriptome))
  }

  seqlevels(transcriptome) <- levels(droplevels(seqnames(transcriptome)))
  export(transcriptome, output_path_transcriptome)
  writeXStringSet(genome, output_path_genome)
  return(0)
}

status <- unify_references(
  vcf = snakemake@input[["vcf"]],
  peptide_bam_homo = snakemake@input[["peptide_bam_homo"]],
  peptide_bam_first = snakemake@input[["peptide_bam_first"]],
  peptide_bam_second = snakemake@input[["peptide_bam_second"]]
)

sessionInfo()

sink()
sink()
