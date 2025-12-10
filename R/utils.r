## Utility Functions

## Inverse Multivariate Logit Transform
expitM <- function(x){
  ex <- exp(c(x,0))
  ex/sum(ex)
}
## Multivariate Logit Transform
logitM <- function(p){
  K <- length(p)
  log(p[1:(K-1)]/p[K])
}

## Transformations for Prior Distributions
## This one if prior is given in terms of p.
## e.g. p ~ beta, logitp ~ beta*|p*(1-p)|
logDetJac_logitM <- function(x){
  sum(log(x))
}

## This one if prior is given in terms of y
## e.g. y ~ gamma, logy ~ gamma*|y|
logDetJac_log <- function(x){
  log(x)
}

## Convert parameter list
reList <- function(parsList, parsVec){
  count <- 1
  for( i in 1:length(parsList) ){
    nL <- length(parsList[[i]])
    parsList[[i]] <- parsVec[count:(count + nL - 1)]
    count <- count + nL
  }
  parsList
}

## Convenience function for extracting control values.
extractControls <- function(controlValue, defaultValue)
{
  if(!is.null(controlValue))  
    return(controlValue)
  else 
    return(defaultValue)
}

## Estimate total species per day / bin for test fishery counts:
estimateDailyN <- function(testFishGrp, prob, q){
  if( nrow(testFishGrp) != nrow(prob) )
  Nspp <- 
  for( i in seq_along(prob) ){
    testFishGrp$totalFish*prob[,i]
  }
      
      idx <- which(tfgrp$day == d)
      nd <- 
      Nd <- sum(nd)
      if(Nd == 0) next
      pd <- numeric(nq)
      for( k in 1:nq ) pd[k] <- sum(p[idx,k+(K-nq)]*nd/Nd)              
      prob <- (q*pd)/sum(q*pd)
      negll <- negll - testwgts[d]*dmultinom(testcounts[d,], prob = prob, size = sum(testcounts[d,]), log = TRUE)
}

## Find and return visible modes based on a density function:
findModes <- function(x,...) {
  modes <- NULL
  f <- density(x, ...)  
  for ( i in 2:(length(f$y)-1) ){
    if ( (f$y[i] > f$y[i-1]) & (f$y[i] > f$y[i+1]) ) {
      modes <- c(modes,i)
    }
  }
  if ( length(modes) == 0 ) {
    return(NA)
  }
  return(f$x[modes])
}

## Find and return the global mode based on a density function:
findGlobalMode <- function(x, xrange = NULL, ...) {
  if(!is.null(xrange)){
    x <- x[x <= xrange[2] & x >= xrange[1]]
  }
  if(length(x) < 3) return(NA)
  f <- density(x, ...) 
  xmax <- f$x[which.max(f$y)]  
  return(xmax)
}

