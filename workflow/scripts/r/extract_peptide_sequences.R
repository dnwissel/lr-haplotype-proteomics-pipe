log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

CODON_NT <- 3
CODON_OFFSET <- 2
MATCH_START <- 1
MATCH_END <- 2

extract_peptide_sequences <- function(search_output,
                                      proteome_fasta,
                                      homozygous_isoforms_fasta,
                                      first_haplotype_isoforms_fasta,
                                      second_haplotype_isoforms_fasta,
                                      contaminants_path,
                                      output_path_homozygous,
                                      output_path_first_haplotype,
                                      output_path_second_haplotype,
                                      fdr = 0.01) {
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(Biostrings)
    library(vroom)
    library(dplyr)
    library(stringr)
  })

  contam <- readAAStringSet(contaminants_path)
  contam_mask <- paste0("(", paste0(c(
  sapply(strsplit(names(contam)[-length(names(contam))], "\\|"), function(x) x[[2]]),
  names(contam)[length(names(contam))]

), collapse = "|"), ")")
  results <- vroom::vroom(search_output) %>% dplyr::filter(peptide_q < fdr & label == 1) %>% filter(!grepl(contam_mask, proteins))
  results$peptide <- str_replace_all(
    results$peptide,
    "\\[(.*?)\\]",
    ""
  )
  proteome <- readAAStringSet(proteome_fasta)

  homozygous_isoforms <- readDNAStringSet(homozygous_isoforms_fasta)
  first_haplotype_isoforms <- readDNAStringSet(first_haplotype_isoforms_fasta)
  second_haplotype_isoforms <- readDNAStringSet(second_haplotype_isoforms_fasta)
  unique_proteins <- lapply(lapply(strsplit(results$proteins, ";"), function(x) sapply(x, function(y) substr(y, nchar(y), nchar(y)))), unique)
  protein_assignment <- unlist(lapply(unique_proteins, function(x) {
    if (length(x) == 1) {
      if (x == "A") {
        return("A")
      } else if (x == "B") {
        return("B")
      }
    }
    if (length(x) == 2) {
      if (all(sort(x) == c("A", "B"))) {
        return("Homozygous")
      }
    }
    return("Homozygous")
  }))

  matched_proteins <- sapply(strsplit(results$proteins, ";"), function(x) x[[1]])
  homozygous_ix <- which(protein_assignment == "Homozygous")
  first_haplotype_ix <- which(protein_assignment == "A")
  second_haplotype_ix <- which(protein_assignment == "B")
  genomic_peptide_sequences_homozygous <- sapply(homozygous_ix, function(ix) {
    ix_matched_protein <- matched_proteins[ix]
    last_protein_letter <- substr(ix_matched_protein, nchar(ix_matched_protein), nchar(ix_matched_protein))
    protein_sequence <- as.character(proteome[which(names(proteome) == ix_matched_protein)])
    matched_locations <- stringr::str_locate(protein_sequence, results$peptide[ix])
    if (last_protein_letter == "A") {
      substr(as.character(first_haplotype_isoforms[names(first_haplotype_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    } else if (last_protein_letter == "B") {
      substr(as.character(second_haplotype_isoforms[names(second_haplotype_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    } else {
      substr(as.character(homozygous_isoforms[names(homozygous_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    }
  })

  genomic_peptide_sequences_first_haplotype <- sapply(first_haplotype_ix, function(ix) {
    ix_matched_protein <- matched_proteins[ix]
    last_protein_letter <- substr(ix_matched_protein, nchar(ix_matched_protein), nchar(ix_matched_protein))
    protein_sequence <- as.character(proteome[which(names(proteome) == ix_matched_protein)])
    matched_locations <- stringr::str_locate(protein_sequence, results$peptide[ix])
    if (last_protein_letter == "A") {
      substr(as.character(first_haplotype_isoforms[names(first_haplotype_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    } else if (last_protein_letter == "B") {
      substr(as.character(second_haplotype_isoforms[names(second_haplotype_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    } else {
      substr(as.character(homozygous_isoforms[names(homozygous_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    }
  })

  genomic_peptide_sequences_second_haplotype <- sapply(second_haplotype_ix, function(ix) {
    ix_matched_protein <- matched_proteins[ix]
    last_protein_letter <- substr(ix_matched_protein, nchar(ix_matched_protein), nchar(ix_matched_protein))
    protein_sequence <- as.character(proteome[which(names(proteome) == ix_matched_protein)])
    matched_locations <- stringr::str_locate(protein_sequence, results$peptide[ix])
    if (last_protein_letter == "A") {
      substr(as.character(first_haplotype_isoforms[names(first_haplotype_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    } else if (last_protein_letter == "B") {
      substr(as.character(second_haplotype_isoforms[names(second_haplotype_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    } else {
      substr(as.character(homozygous_isoforms[names(homozygous_isoforms) == ix_matched_protein]), matched_locations[MATCH_START] * CODON_NT - CODON_OFFSET, matched_locations[MATCH_END] * CODON_NT)
    }
  })

  genomic_peptide_sequences_first_haplotype <- (DNAStringSet((genomic_peptide_sequences_first_haplotype)))
  genomic_peptide_sequences_second_haplotype <- (DNAStringSet((genomic_peptide_sequences_second_haplotype)))
  names(genomic_peptide_sequences_first_haplotype) <- matched_proteins[first_haplotype_ix]
  names(genomic_peptide_sequences_second_haplotype) <- matched_proteins[second_haplotype_ix]

  writeXStringSet(unique(DNAStringSet(unlist(genomic_peptide_sequences_homozygous))), output_path_homozygous)
  writeXStringSet(unique(genomic_peptide_sequences_first_haplotype), output_path_first_haplotype)
  writeXStringSet(unique(genomic_peptide_sequences_second_haplotype), output_path_second_haplotype)
  return(0)
}

status <- extract_peptide_sequences(
  search_output = snakemake@input[["search_output"]],
  proteome_fasta = snakemake@input[["proteome_fasta"]],
  homozygous_isoforms_fasta = snakemake@input[["homozygous_isoforms_fasta"]],
  first_haplotype_isoforms_fasta = snakemake@input[["first_haplotype_isoforms_fasta"]],
  second_haplotype_isoforms_fasta = snakemake@input[["second_haplotype_isoforms_fasta"]],
  contaminants_path = snakemake@input[["contaminants_path"]],
  output_path_homozygous = snakemake@output[["output_path_homozygous"]],
  output_path_first_haplotype = snakemake@output[["output_path_first_haplotype"]],
  output_path_second_haplotype = snakemake@output[["output_path_second_haplotype"]],
  fdr = snakemake@params[["fdr"]]
)

sessionInfo()

sink()
sink()
