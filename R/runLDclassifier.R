#' @export
#'
#' @param haplos Matrix with the haplotypes (SNPs in rows, samples in columns)
#' @param annot GenomicRanges
#' @param blockSize Numeric with the size of the SNP block
runLDclassifier <- function(haplos, annot, blockSize = 2, mc.cores = 1){

  # Make list of SNP pairs to test
  GRblocks <- GenomicRanges::GRanges(
    paste0(GenomicRanges::seqnames(annot)[1], ":",
           GenomicRanges::start(annot)[1:(length(annot) - blockSize + 1)], "-",
           GenomicRanges::end(annot)[blockSize:(length(annot))]))
  GRblocks$Ind <- lapply(1:(length(annot) - blockSize + 1),
                         function(x) seq(x, x + blockSize -1, 1))
  overlaps <- GenomicRanges::findOverlaps(GenomicRanges::resize(GRblocks, 1e5), GRblocks)
  autoOverlaps <- GenomicRanges::findOverlaps(GRblocks)
  overlaps <- overlaps[!S4Vectors::`%in%`(overlaps, autoOverlaps)]

  models <- parallel::mclapply(seq_len(length(overlaps)), function(ind){
    bl1 <- S4Vectors::from(overlaps)[ind]
    bl2 <- S4Vectors::to(overlaps)[ind]
    ind1 <- GRblocks$Ind[[bl1]]
    ind2 <- GRblocks$Ind[[bl2]]
    res <- inversionModel(cbind(
      apply(haplos[, ind1], 1, function(x) paste(x, collapse = "")),
      apply(haplos[, ind2], 1, function(x) paste(x, collapse = ""))
    ))
    res$annot <- c(start = GenomicRanges::start(GRblocks[bl1]),
                   end = GenomicRanges::start(GRblocks[bl2]))
    if (res$bic < 10 | res$pval < 0.05) {
      res$r1 <- NULL
    }
    res
  }, mc.cores = mc.cores)
  models
}
