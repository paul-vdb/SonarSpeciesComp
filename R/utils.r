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