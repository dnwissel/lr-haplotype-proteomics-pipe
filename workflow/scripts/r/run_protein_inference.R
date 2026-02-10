log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

run_protein_inference <- function(sage_tsv, contaminants_path, output_path, fdr_cutoff, seed) {
  set.seed(seed)
  results <- vroom::vroom(sage_tsv) %>%
    filter(peptide_q < fdr_cutoff) %>%
    filter(label == 1)

  contaminants <- readAAStringSet(contaminants_path)
  grep_contaminant_names <- sapply(strsplit(names(contaminants), "\\|"), function(x) x[length(x)])
  contaminant_indices <- grep(paste0(grep_contaminant_names, collapse = "|"), results$proteins)
  results <- results[-contaminant_indices, ]

  # Collapse duplicates
  results <- results %>%
    select(peptide, proteins) %>%
    distinct()

  # Collapse peptides that group to the same group of proteins
  results_filtered <- results %>%
    group_by(proteins) %>%
    slice_head(n = 1) %>%
    ungroup()
    #summarise(peptide = first(peptide))

  print(results_filtered)
  possible_proteins <- unlist(strsplit(results_filtered$proteins, ";"))
  possible_peptides <- results_filtered$peptide

  protein_matches <- strsplit(results_filtered$proteins, ";")

  protein_mapping <- data.frame(lapply(possible_proteins, function(protein) sapply(protein_matches, function(x) as.integer(protein %in% x))))
  protein_mapping <- apply(protein_mapping, 2, function(x) which(x != 0))
  mapping_table <- table(sapply(protein_mapping, function(x) paste0(sort(x), collapse = ",")))
  group_connections <- lapply(strsplit(names(mapping_table), ","), as.numeric)
  set <- rep(1:length(group_connections), times = sapply(group_connections, length))
  element <- unlist(group_connections)
  X <- data.table::data.table(
    set = set,
    element = element,
    key = c("set", "element")
  )

  res <- RcppGreedySetCover::greedySetCover(X)

  fixed_names <- sapply(protein_mapping, function(x) paste0(sort(x), collapse = ","))

  inferred_proteins <- unlist(sapply(names(mapping_table[unique(res$set)]), function(x) {
    paste0(unique(possible_proteins[which(x == fixed_names)]), collapse = "-")
  }))
  data.frame(protein_group = inferred_proteins) %>% write_tsv(output_path)
  return(0)
}

suppressPackageStartupMessages({
  library(RcppGreedySetCover)
  library(vroom)
  library(dplyr)
  library(readr)
  library(Biostrings)
  library(data.table)
})

status <- run_protein_inference(
  sage_tsv = snakemake@input[["sage_tsv"]],
  contaminants_path = snakemake@input[["contaminants_path"]],
  output_path = snakemake@output[["output_path"]],
  fdr_cutoff = snakemake@params[["fdr_cutoff"]],
  seed = snakemake@params[["seed"]]
)

sessionInfo()

sink()
sink()
