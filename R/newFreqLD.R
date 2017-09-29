
newFreqLD <-function(Resp, Block)
{
  nSNP <- (nchar(Block[1]) - 1)/2

  ## Compute blocks frequency
  freqs <-tapply(Resp, Block, sum)

  # Get all possible left and right levels
  ## Make all combination with 0 and 1
  combs <- expand.grid(lapply(1:nSNP, function(x) c(0, 1)))
  ## Bind to get levels
  leftLevs <- rightLevs <- apply(combs, 1, paste, collapse = "")

  perms <- gtools::permutations(2^(nSNP), 2^(nSNP))

  rightMat <- matrix(rightLevs[perms], ncol = 2^nSNP)
  leftMat <- matrix(leftLevs, nrow = nrow(rightMat), ncol = ncol(rightMat), byrow = TRUE)
  permsmat <- matrix(paste(leftMat, rightMat, sep = "+"), ncol = ncol(rightMat))
  ## Find combination containing more samples
  permsSum <- apply(permsmat, 1, function(x) sum(freqs[x], na.rm = TRUE))
  selAlleles <- permsmat[which.max(permsSum), ]
  props <- freqs[selAlleles]
  props[is.na(props)] <- 0
  names(props)[is.na(names(props))] <- selAlleles[!selAlleles %in% names(props)]
  props <- sort(props, decreasing = TRUE)

  #solve linear algebra to find new frequency values
  AA <- diag(2^nSNP)
  AA[,1] <- -props/max(props)
  AA[1, ] <- 1

  bb <- rep(0, nrow(AA))
  bb[1] <- 1

  res <- qr.solve(AA, bb, tol = 1e-10)
  names(res) <- names(props)

  ## Get list of all combinations
  namesVec <- expand.grid(leftLevs, rightLevs)
  namesVec <- paste(namesVec[, 1], namesVec[, 2], sep = "+")

  res <- c(res, rep(0, length(namesVec) - length(props)))
  res[res < 1e-3] <- 1e-3
  res <- res/sum(res)
  names(res)[(length(props) +1):length(res)] <-
    namesVec[!namesVec %in% names(res)]
  return(res)
}
