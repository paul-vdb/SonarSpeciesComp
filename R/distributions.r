## Distributions to be used for species composition:

## Compute likelihood and posterior probabilities of mixture model.
calcPostProb <- function(x, mu, sigma, prob, wgts = NULL){
  n <- length(x)
  K <- length(mu)
  if(is.null(wgts)) wgts <- rep(1, n)  
  postp <- matrix(0, nrow = n, ncol = K)
  ll <- 0
  
  for( i in 1:n ) {
    if(is.matrix(prob)){ logpostp <- log(prob[i,]) + dnorm(x[i], mu, sigma, log =TRUE)
    }else{ logpostp <- log(prob) + dnorm(x[i], mu, sigma, log =TRUE) }
    C <- max(logpostp)
    logrowsum <- log(sum(exp(logpostp - C))) + C
    postp[i,] <- exp(logpostp - logrowsum)
    ll <- ll + logrowsum*wgts[i]
  }
  list(postp = postp, ll = ll)
}

dTestFishery <- function(x, q, ){
  if(is.matrix(x)){
    ndays <- nrow(x)
    nspp <- ncol(x)
  }else{ 
    ndays <- 1
    nspp <- length(x)
    x <- matrix(x, nrow = 1, ncol = nspp)
  }
  for( d in 1:ndays ){
      idx <- which(tfgrp$day == d)
      nd <- tfgrp$pcount[idx]
      Nd <- sum(nd)
      if(Nd == 0) next
      pd <- numeric(nq)
      for( k in 1:nq ) pd[k] <- sum(p[idx,k+(K-nq)]*nd/Nd)              
      prob <- (q*pd)/sum(q*pd)
      negll <- negll - testwgts[d]*dmultinom(testcounts[d,], prob = prob, size = sum(testcounts[d,]), log = TRUE)
    }
}

## Compute likelihood conditional on posterior probabilities of mixture model.
dnormEM <- function(x, mu, sigma, p, postp, wgts = NULL, log){
    K <- length(mu)
    n <- length(x)
    if(is.null(wgts)) wgts <- rep(1, n)
    logobj <- 0
    for( k in 1:K ){
      logobj <- logobj + sum(wgts*postp[,k]*(logp[grp,k] + dnorm(x, mu[k], sigma[k], log = TRUE)))
    }
    logobj
}

## Negative log likelihood for obtaining confidence intervals:
## To be used with MakeADFun after fitting EM.
## This function will have jacobian options, 
## as it is what would be called by tmbstan.
negLL <- function(pars){
  getAll(pars, pars_fixed, warn = FALSE)
  p <- expitM(logitp)
  logp <- log(p)
  sigma <- exp(logsigma)

  K <- length(sigma)
  n <- length(dataList$x)

  negll <- 0
  for( i in 1:n ){
      logprob <- logp + dnorm(dataList$x[i], mu, sigma, log = TRUE)
      maxlp <- max(logprob)
      negll <- negll - dataList$wgts[i]*(log(sum(exp(logprob - maxlp))) + maxlp)
  }
  negll
}

## Place holders for penalty/prior terms:
prior_beta <- function(beta){
  0
}
prior_sigma0 <- function(sigma0){
  0
}
## This prior includes chinook jack and adult.
prior_sigma <- function(sigma){
  0
}
prior_q <- function(q){
  0
}
## This prior includes chinook jack and adult.
prior_mu <- function(sigma0){
  0
}
