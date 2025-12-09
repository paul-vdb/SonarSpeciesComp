## Distributions to be used for species composition:

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

dresident <- dnorm
dpink <- dnorm
dsockeye <- dnorm
dchinook <- dnorm

dmixture <- function(x, mu, sigma, spp == NULL, log){
  if(is.null(spp)) return(dnorm(x, mu, sigma, log))
  for( i in seq_along() ){
    dresident(x, mu)
  }
}