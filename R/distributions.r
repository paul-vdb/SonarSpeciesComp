#' Calculate the Posterior Probability of Mixture Proportions
#'
#' Used in the Expectation-Maximization Algorithm for fitting the joint mixture model.
#'
#' @param x vector of observations (e.g. fish lengths).
#' @param mu vector of means for each component (e.g. species).
#' @param sigma vector of standard deviations for each componenet (e.g. species).
#' @param prob matrix of mixture weights for each observation and component. Often constant but may vary due to covariates.
#' @wgts vector of likelihood weights for each observation.
#'
#' @return list of posterior probabilities and log likelihood of mixture component.
#'
#'
#' @export
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

#' Calculate the Mixture Proportions for Prediction in model:
#'
#' Used in the Expectation-Maximization Algorithm for fitting the joint mixture model.
#'
#' @param alpha vector of linear predictors for each species proportions.
#' @param alpha vector of linear predictors for jack Chinook as a proportion of all Chinook.
#' @param dataList dataList used in the EM algorithm.
#' @param pAdultChinook proportion of adult Chinook that are small or large. Defaults to fixed = 1 if adults are just one group.
#' @K0 Number of non-Chinook specices
#' @K Number of total categories (species + Chinook types).
#'
#' @return list of posterior probabilities and log likelihood of mixture component.
#'
#'
#' @export
calcProportions <- function(alpha, alphaJackChinook, pAdultChinook, dataList, K0, K){
  ## Set up proportions. alpha parameter for predicting proportions. 
  ## Xalpha is a list of design matrices for each species.
  np <- nrow(dataList$predDF)
  logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
  indx0 <- 1
  for( i in 1:K0 ){
    nc <- ncol(dataList$Xprop[[i]]) 
    indx1 <- indx0 + nc - 1
    logitp[,i] <- as.matrix(dataList$Xprop[[i]]) %*% alpha[indx0:indx1]
    indx0 <- indx1 + 1
  }

  logitpjack <- as.matrix(dataList$XpropChin) %*% alphaJackChinook
  pjack <- 1/(1+exp(-logitpjack))

  p <- matrix(0, nrow = np, ncol = K)
  for( i in 1:np ) {
    p[i, 1:(K0+1)] <- expitM(logitp[i,])  
    p[i, (K0+1):K ] <- p[i, (K0+1)]*c(pjack[i], (1-pjack[i])*pAdultChinook)
  }
  return(p)
}

calcTestFisheryProbs <- function(p, q, dataList, K0, K){
  ndays <- nrow(dataList$testFisheryCounts)
  padultchinook <- 1-rowSums(p[,1:(K0+1)])
  nq <- ncol(dataList$testFisheryCounts)
  for( d in 1:ndays ){
    idx <- which(dataList$predDF$day == d)
    Nd <- sapply(dataList$qSppIndices, FUN = function(x){sum(p[idx,x]*dataList$predDF$SalmonCountTF[idx])})
    Ndc <- sum(padultchinook[idx]*dataList$predDF$SalmonCountTF[idx])
    Nd <- c(Nd, Ndc)
    # Nsalmon <- sum(dataList$lengthData$SalmonCount[idx]/dataList$lengthData$nLengths[idx])
    prob <- (q*Nd)/sum(q*Nd)
  }
  return(prob)
}

## Delete:
#' Calculate the density of the test fishery counts
#'
#' Computes probability of the test fishery counts for some relative catchability value.
#'
#' @param x vector of test fishery counts.
#' @param q vector of relative catchability.
#'
#' @return doesn't return anything...
#'
#'
### @export
# dTestFishery <- function(x, q ){
  # if(is.matrix(x)){
    # ndays <- nrow(x)
    # nspp <- ncol(x)
  # }else{ 
    # ndays <- 1
    # nspp <- length(x)
    # x <- matrix(x, nrow = 1, ncol = nspp)
  # }
  # for( d in 1:ndays ){
      # idx <- which(tfgrp$day == d)
      # nd <- tfgrp$pcount[idx]
      # Nd <- sum(nd)
      # if(Nd == 0) next
      # pd <- numeric(nq)
      # for( k in 1:nq ) pd[k] <- sum(p[idx,k+(K-nq)]*nd/Nd)              
      # prob <- (q*pd)/sum(q*pd)
      # negll <- negll - testwgts[d]*dmultinom(testcounts[d,], prob = prob, size = sum(testcounts[d,]), log = TRUE)
    # }
# }

#' Calculate the normal distribution int he context of the EM algorithm.
#'
#' Used in the Expectation-Maximization Algorithm for fitting the joint mixture model.
#'
#' @param x vector of observations (e.g. fish lengths).
#' @param mu vector of means for each component (e.g. species).
#' @param sigma vector of standard deviations for each componenet (e.g. species).
#' @param prob matrix of mixture weights for each observation and component. Often constant but may vary due to covariates.
#' @wgts vector of likelihood weights for each observation.
#'
#' @return list of posterior probabilities and log likelihood of mixture component.
#'
#'
# dnormEM <- function(x, mu, sigma, p, postp, wgts = NULL, log){
    # K <- length(mu)
    # n <- length(x)
    # if(is.null(wgts)) wgts <- rep(1, n)
    # logobj <- 0
    # for( k in 1:K ){
      # logobj <- logobj + sum(wgts*postp[,k]*(logp[grp,k] + dnorm(x, mu[k], sigma[k], log = TRUE)))
    # }
    # logobj
# }

## Negative log likelihood for obtaining confidence intervals:
## To be used with MakeADFun after fitting EM.
## This function will have jacobian options, 
## as it is what would be called by tmbstan.
# negLL <- function(pars){
  # getAll(pars, pars_fixed, warn = FALSE)
  # p <- expitM(logitp)
  # logp <- log(p)
  # sigma <- exp(logsigma)

  # K <- length(sigma)
  # n <- length(dataList$x)

  # negll <- 0
  # for( i in 1:n ){
      # logprob <- logp + dnorm(dataList$x[i], mu, sigma, log = TRUE)
      # maxlp <- max(logprob)
      # negll <- negll - dataList$wgts[i]*(log(sum(exp(logprob - maxlp))) + maxlp)
  # }
  # negll
# }
