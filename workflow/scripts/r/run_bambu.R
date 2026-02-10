log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

run_bambu <- function(reads, annotations, genome, output_path, ncore, NDR, quantification, lowMemory = TRUE) {
  bambuAnnotations <- prepareAnnotations(annotations)
  bambu_result <- bambu::bambu(
    reads = reads, annotations = bambuAnnotations,
    genome = genome, discovery = TRUE,
    ncore = ncore, lowMemory = lowMemory, NDR = NDR, quant = FALSE
  )
  writeToGTF(bambu_result, file = output_path)
  return(0)
}

suppressPackageStartupMessages(library(bambu))

status <- run_bambu(
  reads = snakemake@input[["reads"]],
  annotations = snakemake@input[["annotations"]],
  genome = snakemake@input[["genome"]],
  output_path = snakemake@output[["output_path"]],
  ncore = snakemake@threads,
  lowMemory = TRUE,
  NDR = NULL
)

sessionInfo()

sink()
sink()
