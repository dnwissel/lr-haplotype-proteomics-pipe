log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

install_rcppgreedysetcover <- function(version, outfile) {
  renv::install(paste0("RcppGreedySetCover@", version), type = "binary")
  # https://stackoverflow.com/questions/23922497/create-a-touch-file-on-unix
  write.table(data.frame(), file = outfile, col.names = FALSE)
  return(0)
}

suppressPackageStartupMessages({
  library(renv)
})


status <- install_rcppgreedysetcover(
  version = snakemake@params[["version"]],
  outfile = snakemake@output[[1]]
)

sessionInfo()

sink()
sink()
