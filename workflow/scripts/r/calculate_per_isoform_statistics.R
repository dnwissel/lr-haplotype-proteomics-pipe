log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

calculate_per_isoform_statistics <- function(sqanti_protein, haplosaurus_json, output_path) {
  suppressPackageStartupMessages({
    library(vroom)
    library(dplyr)
    library(rjson)
    library(readr)
  })
  ### Approach from SO: https://stackoverflow.com/questions/60298776/convert-json-file-with-multiple-lines-to-r-dataframe
  json_raw <- readr::read_file(haplosaurus_json)
  json_lines <- unlist(strsplit(json_raw, "\\n"))
  json_parsed <- lapply(json_lines, jsonlite::fromJSON)

  protein_class <- vroom::vroom(sqanti_protein)
  protein_class[which(!grepl("BambuTx", protein_class$pb)), ]$protein_classification_base <- "pFSM"

  protein_class[protein_class$protein_classification_base == "orphan_monoexon", ]$protein_classification_base <- "OE"


  n_variants <- unlist(lapply(json_parsed, function(x) {
    if (length(x$protein_haplotypes$contributing_variants) == 2) {
      return(c(
        length(x$protein_haplotypes$contributing_variants[[1]]),
        length(x$protein_haplotypes$contributing_variants[[2]])
      ))
    } else {
      return(length(x$protein_haplotypes$contributing_variants[[1]]))
    }
  }))


  n_indels <- unlist(lapply(
    1:length(json_parsed),
    function(transcript) {
      if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants) == 2) {
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]]) == 0) {
          first_indel_count <- 0
        } else {
          first_indel_count <- sum(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]], "_"), function(x) sum(max(nchar(strsplit(x[[length(x)]], "/")[[1]])) > 1)))
        }
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[2]]) == 0) {
          second_indel_count <- 0
        } else {
          second_indel_count <- sum(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[2]], "_"), function(x) sum(max(nchar(strsplit(x[[length(x)]], "/")[[1]])) > 1)))
        }
        return(
          c(
            first_indel_count,
            second_indel_count
          )
        )
      } else {
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]]) == 0) {
          return(0)
        } else {
          return(
            sum(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]], "_"), function(x) sum(max(nchar(strsplit(x[[length(x)]], "/")[[1]])) > 1)))
          )
        }
      }
    }
  ))

  n_snvs <- n_variants - n_indels



  diffs <- unlist(lapply(json_parsed, function(x) {
    if (length(x$protein_haplotypes$seq) > 1) {
      if (length(x$protein_haplotypes$diffs[[1]]) == 0) {
        first_diff <- ""
      } else {
        first_diff <- paste0(unlist(x$protein_haplotypes$diffs[[1]]), collapse = ";")
      }
      if (length(x$protein_haplotypes$diffs[[2]]) == 0) {
        second_diff <- ""
      } else {
        second_diff <- paste0(unlist(x$protein_haplotypes$diffs[[2]]), collapse = ";")
      }
      return(c(
        first_diff, second_diff
      ))
    } else {
      if (length(x$protein_haplotypes$diffs[[1]]) == 0) {
        first_diff <- ""
      } else {
        first_diff <- paste0(unlist(x$protein_haplotypes$diffs[[1]]), collapse = ";")
      }
      return(first_diff)
    }
  }))

  contributing_variants <- unlist(lapply(json_parsed, function(x) {
    if (length(x$protein_haplotypes$seq) > 1) {
      if (length(x$protein_haplotypes$diffs[[1]]) == 0) {
        first_diff <- ""
      } else {
        first_diff <- paste0(unlist(x$protein_haplotypes$contributing_variants[[1]]), collapse = ";")
      }
      if (length(x$protein_haplotypes$diffs[[2]]) == 0) {
        second_diff <- ""
      } else {
        second_diff <- paste0(unlist(x$protein_haplotypes$contributing_variants[[2]]), collapse = ";")
      }
      return(c(
        first_diff, second_diff
      ))
    } else {
      if (length(x$protein_haplotypes$diffs[[1]]) == 0) {
        first_diff <- ""
      } else {
        first_diff <- paste0(unlist(x$protein_haplotypes$contributing_variants[[1]]), collapse = ";")
      }
      return(first_diff)
    }
  }))

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

  has_changed_stop <- as.integer(unlist(lapply(json_parsed, function(x) {
    if (length(x$protein_haplotypes$seq) > 1) {
      if (length(x$protein_haplotypes$flags[[1]]) == 0) {
        first_flag <- FALSE
      } else {
        first_flag <- any(grepl("stop_change", x$protein_haplotypes$flags[[1]]))
      }
      if (length(x$protein_haplotypes$flags[[2]]) == 0) {
        second_flag <- FALSE
      } else {
        second_flag <- any(grepl("stop_change", x$protein_haplotypes$flags[[2]]))
      }
      return(c(
        first_flag, second_flag
      ))
    } else {
      if (length(x$protein_haplotypes$flags[[1]]) == 0) {
        first_flag <- FALSE
      } else {
        first_flag <- any(grepl("stop_change", x$protein_haplotypes$flags[[1]]))
      }
      return(first_flag)
    }
  })))

  has_frameshift <- unlist(lapply(
    1:length(json_parsed),
    function(transcript) {
      if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants) == 2) {
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]]) == 0) {
          first_indel_count <- FALSE
        } else {
          first_indel_count <- any(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]], "_"), function(x) any((abs(diff(nchar(strsplit(x[[length(x)]], "/")[[1]]))) %% 3) != 0)))
        }
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[2]]) == 0) {
          second_indel_count <- FALSE
        } else {
          second_indel_count <- any(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[2]], "_"), function(x) any((abs(diff(nchar(strsplit(x[[length(x)]], "/")[[1]]))) %% 3) != 0)))
        }
        return(
          c(
            first_indel_count,
            second_indel_count
          )
        )
      } else {
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]]) == 0) {
          return(FALSE)
        } else {
          return(
            any(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]], "_"), function(x) any((abs(diff(nchar(strsplit(x[[length(x)]], "/")[[1]]))) %% 3) != 0)))
          )
        }
      }
    }
  ))

  has_resolved_frameshift <- unlist(lapply(
    1:length(json_parsed),
    function(transcript) {
      if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants) == 2) {
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]]) == 0) {
          first_indel_count <- FALSE
        } else {
          first_indel_count <- sum(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]], "_"), function(x) ((diff(nchar(strsplit(x[[length(x)]], "/")[[1]])))))) == 0
        }
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[2]]) == 0) {
          second_indel_count <- FALSE
        } else {
          second_indel_count <- sum(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[2]], "_"), function(x) ((diff(nchar(strsplit(x[[length(x)]], "/")[[1]])))))) == 0
        }
        return(
          c(
            first_indel_count,
            second_indel_count
          )
        )
      } else {
        if (length(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]]) == 0) {
          return(FALSE)
        } else {
          return(
            sum(sapply(strsplit(json_parsed[[transcript]]$protein_haplotypes$contributing_variants[[1]], "_"), function(x) ((diff(nchar(strsplit(x[[length(x)]], "/")[[1]])))))) == 0
          )
        }
      }
    }
  ))



  has_resolved_frameshift <- has_resolved_frameshift & has_frameshift

  has_frameshift <- has_frameshift & (!has_resolved_frameshift)




  annotation_frame <- data.frame(
    isoform_id = isoform_id,
    n_variants = n_variants,
    n_snvs = n_snvs,
    n_indels = n_indels,
    contributing_variants = contributing_variants,
    diff = diffs,
    has_indel = as.integer(n_indels > 0),
    has_changed_stop = has_changed_stop,
    has_frameshift = as.integer(has_frameshift),
    has_resolved_frameshift = as.integer(has_resolved_frameshift),
    match_id = gsub("(A$|B$)", "", isoform_id)
  )

  has_indel <- as.integer(n_indels > 0)

  append_frame <- data.frame(
    isoform_id = protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)],
    protein_splice_category = protein_class$protein_classification_base[which(!protein_class$pb %in% annotation_frame$isoform_id)],
    n_variants = rep(0, length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)])),
    n_snvs = rep(0, length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)])),
    n_indels = rep(0, length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)])),
    contributing_variants = rep("", length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)])),
    diff = rep("", length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)])),
    has_indel = rep(0, length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)])),
    has_changed_stop = rep(0, length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)])),
    has_frameshift = rep(0, length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)])),
    has_resolved_frameshift = rep(0, length(protein_class$pb[which(!protein_class$pb %in% annotation_frame$isoform_id)]))
  )

  join_frame <- data.frame(
    isoform_id = protein_class$pb[which(protein_class$pb %in% annotation_frame$isoform_id)],
    protein_splice_category = protein_class$protein_classification_base[which(protein_class$pb %in% annotation_frame$isoform_id)]
  )

  annotation_frame %>%
    left_join(
      join_frame,
      by = c("isoform_id" = "isoform_id")
    ) %>%
    select(
      `isoform_id`,
      `protein_splice_category`,
      `n_variants`,
      `n_snvs`,
      `n_indels`,
      `contributing_variants`,
      `diff`,
      `has_indel`,
      `has_changed_stop`,
      `has_frameshift`,
      `has_resolved_frameshift`
    ) %>%
    rbind(append_frame) %>%
    write_tsv(output_path)

  return(0)
}

status <- calculate_per_isoform_statistics(
  haplosaurus_json = snakemake@input[["haplosaurus_json"]],
  sqanti_protein = snakemake@input[["sqanti_protein"]],
  output_path = snakemake@output[[1]]
)

sessionInfo()

sink()
sink()
