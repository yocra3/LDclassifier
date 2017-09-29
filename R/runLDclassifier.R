#' @export
#'
#' @param haplos Matrix with the haplotypes (SNPs in rows, samples in columns)
#' @param annot GenomicRanges
#' @param blockSize Numeric with the size of the SNP block
runLDclassifier <- function(haplos, annot, blockSize = 2, mc.cores = 1){

  # Make list of SNP pairs to test
  GRblocks <- GRanges(
    paste0(seqnames(annot)[1], ":",
           start(annot)[1:(length(annot) - blockSize + 1)], "-",
           end(annot)[blockSize:(length(annot))]))
  GRblocks$Ind <- lapply(1:(length(annot) - blockSize + 1),
                         function(x) seq(x, x + blockSize -1, 1))
  overlaps <- findOverlaps(resize(GRblocks, 1e5), GRblocks)
  autoOverlaps <- findOverlaps(GRblocks)
  overlaps <- overlaps[!overlaps %in% autoOverlaps]

  message("Pairs done")

  models <- mclapply(seq_len(length(overlaps)), function(ind){
    bl1 <- from(overlaps)[ind]
    bl2 <- to(overlaps)[ind]
    ind1 <- GRblocks$Ind[[bl1]]
    ind2 <- GRblocks$Ind[[bl2]]
    res <- inversionModel(cbind(
      apply(haplos[, ind1], 1, function(x) paste(x, collapse = "")),
      apply(haplos[, ind2], 1, function(x) paste(x, collapse = ""))
    ))
    res$annot <- c(start = start(GRblocks[bl1]), end = start(GRblocks[bl1]))
    if (res$bic < 50) {
      res$r1 <- NULL
    }
    res
  }, mc.cores = mc.cores)
  models
}
