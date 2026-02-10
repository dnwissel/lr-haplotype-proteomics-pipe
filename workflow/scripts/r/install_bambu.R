log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

install_bambu <- function(version, outfile) {
  # Install XML separately to avoid compilation issues
  #renv::install("XML", type = "binary")
  #renv::install("bioc::GenomeInfoDb")
  #renv::install(paste0("bioc::bambu@", version), type = "binary")
  # https://stackoverflow.com/questions/23922497/create-a-touch-file-on-unix
  write.table(data.frame(), file = outfile, col.names = FALSE)
  return(0)
}

#suppressPackageStartupMessages({
#  library(renv)
#  library(BiocManager)
#})


status <- install_bambu(
  version = snakemake@params[["version"]],
  outfile = snakemake@output[[1]]
)

sessionInfo()

sink()
sink()
