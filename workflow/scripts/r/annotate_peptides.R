log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

annotate_peptides <- function(input_peptides, input_cds_homo,
                              input_cds_first, input_cds_second,
                              reference_transcriptome,
                              reference_extracted_cds,
                              output_path) {
  suppressPackageStartupMessages({
    library(jsonlite)
    library(readr)
    library(Biostrings)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(stringr)
  })
  peptide_sequences <- readDNAStringSet(input_peptides)
  #print(peptide_sequences)
  #print(length(peptide_sequences))
  #stop("")
  if (length(peptide_sequences) == 0) {
  	data.frame() %>% write_tsv(output_path)
 	return(0)
  }
  cds_sequences <- c(
    readDNAStringSet(input_cds_homo),
    readDNAStringSet(input_cds_first),
    readDNAStringSet(input_cds_second)
  )
  reference_cds_sequences <- readDNAStringSet(reference_extracted_cds)
  transcriptome <- import(reference_transcriptome)

  non_bambu_transcriptome <- transcriptome[!grepl("Bambu", transcriptome$transcript_id) & transcriptome$type == "CDS"]
  bambu_transcriptome <- transcriptome[grepl("Bambu", transcriptome$transcript_id) & transcriptome$type == "CDS"]
  unique_bambu_cds_regions <- setdiff(bambu_transcriptome, non_bambu_transcriptome)

  output_peptides <- lapply(
    1:length(peptide_sequences),
    function(peptide_ix) {
      peptide <- peptide_sequences[peptide_ix]
      transcript_id <- sapply(strsplit(names(peptide), "\\."), function(x) {
        if (length(x) == 2) {
          return(paste0(x[1:2], collapse = "."))
        } else {
          return(paste0(x, collapse = "."))
        }
      })
      non_het_transcript_id <- sub("[AB]$", "", transcript_id)
      non_het_gene_id <- transcriptome$gene_id[which(transcriptome$transcript_id == non_het_transcript_id)][1]
      transcript_cds <- transcriptome[transcriptome$transcript_id == non_het_transcript_id & transcriptome$type == "CDS"]

      start_stop <- str_locate(
        as.character(cds_sequences[which(names(cds_sequences) == transcript_id)]),
        as.character(peptide)
      )

      reference_equivalent_peptide <- substr(as.character(reference_cds_sequences[which(names(reference_cds_sequences) == non_het_transcript_id)]), start_stop[1], start_stop[2])

      transcript_transcriptome <- transcriptome[transcriptome$type == "transcript"]
      gene_transcripts_without_peptide_transcript <- transcript_transcriptome[transcript_transcriptome$gene_id == non_het_gene_id & transcript_transcriptome$transcript_id != non_het_transcript_id]$transcript_id
      gene_transcripts_without_peptide_transcript <- gene_transcripts_without_peptide_transcript[gene_transcripts_without_peptide_transcript %in% names(reference_cds_sequences) & !grepl("Bambu", gene_transcripts_without_peptide_transcript)]
      joint_gene_cds <- paste0(reference_cds_sequences[gene_transcripts_without_peptide_transcript], collapse = "")
      c(
        !str_detect(
          joint_gene_cds,
          as.character(peptide)
        ) & grepl("Bambu", non_het_transcript_id),
        !str_detect(
          as.character(reference_cds_sequences[which(names(reference_cds_sequences) == non_het_transcript_id)]),
          as.character(peptide)
        )
      )
    }
  )
  names(output_peptides) <- names(peptide_sequences)
  data.frame(output_peptides) %>% write_tsv(output_path)
  return(0)
}

status <- annotate_peptides(
  input_peptides = snakemake@input[["input_peptides"]],
  input_cds_homo = snakemake@input[["input_cds_homo"]],
  input_cds_first = snakemake@input[["input_cds_first"]],
  input_cds_second = snakemake@input[["input_cds_second"]],
  reference_transcriptome = snakemake@input[["reference_transcriptome"]],
  reference_extracted_cds = snakemake@input[["reference_extracted_cds"]],
  output_path = snakemake@output[[1]]
)

sessionInfo()

sink()
sink()
