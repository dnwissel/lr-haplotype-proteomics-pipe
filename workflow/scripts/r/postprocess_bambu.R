log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

postprocess_bambu <- function(novel_transcriptome, filtered_transcriptome, output_path) {
  novel_transcriptome <- import(novel_transcriptome)
  filtered_transcriptome <- import(filtered_transcriptome)

  novel_transcriptome <- novel_transcriptome[
    seqnames(novel_transcriptome) %in% paste0("chr", 1:22) &
    novel_transcriptome$gene_id %in% filtered_transcriptome$gene_id
  ]

  novel_isoforms <- unique(novel_transcriptome$transcript_id[which(substr(novel_transcriptome$transcript_id, 1, 7) == "BambuTx")])
  known_isoforms_kept <- unique(novel_transcriptome$transcript_id[which(substr(novel_transcriptome$transcript_id, 1, 7) != "BambuTx")])
  known_isoforms_kept_cds <- known_isoforms_kept[known_isoforms_kept %in% filtered_transcriptome$transcript_id]

  filtered_transcriptome <- filtered_transcriptome[filtered_transcriptome$transcript_id %in% known_isoforms_kept_cds, c("source", "type", "score", "phase", "gene_id", "transcript_id")]
  filtered_transcriptome <- filtered_transcriptome[filtered_transcriptome$type %in% c("transcript", "CDS", "exon")]
  novel_transcriptome <- novel_transcriptome[novel_transcriptome$transcript_id %in% novel_isoforms, c("source", "type", "score", "phase", "gene_id", "transcript_id")]
  joint_transcriptome <- c(filtered_transcriptome, novel_transcriptome)
  export(joint_transcriptome, output_path)
  return(0)
}

suppressPackageStartupMessages(library(rtracklayer))

status <- postprocess_bambu(
  novel_transcriptome = snakemake@input[["novel_transcriptome"]],
  filtered_transcriptome = snakemake@input[["filtered_transcriptome"]],
  output_path = snakemake@output[["output_path"]]
)

sessionInfo()

sink()
sink()
