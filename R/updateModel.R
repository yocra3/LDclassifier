
updateModel <- function(params, funcProps1, funcProps2){

  newparams <- params

  #compute responsabilities
  newparams$r1 <- params$prob0*params$props1[params$inds]
  newparams$r2 <- (1 - params$prob0)*params$props2[params$inds]

  R1 <- newparams$r1/(newparams$r1 + newparams$r2)
  R2 <- newparams$r2/(newparams$r1 + newparams$r2)

  #Compute the new probability of NO inversion;
  newparams$prob0 <- mean(R1)

  #compute new frequencies for each haplotype in each population (no-inv and inv)
  newparams$props1 <- funcProps1(R1, newparams$inds)
  newparams$props2 <- funcProps2(R2, newparams$inds)
  newparams
}
