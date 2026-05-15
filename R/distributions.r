#' Calculate the Posterior Probability of Mixture Proportions
#'
#' Used in the Expectation-Maximization Algorithm for fitting the joint mixture model.
#'
#' @param x vector of observations (e.g. fish lengths).
#' @param mu vector of means for each component (e.g. species).
#' @param sigma vector of standard deviations for each componenet (e.g. species).
#' @param prob matrix of mixture weights for each observation and component. Often constant but may vary due to covariates.
#' @param wgts vector of likelihood weights for each observation.
#'
#' @details Computes the posterior probability of each x belonging to component k for a given mean (`mu`), standard deviation (`sigma`), and mixture weight (`prob`). 
#' Allows for likelihood weights to be passed to compute the log-likelihood of the mixture component. 
#'
#' @return list of posterior probabilities and log likelihood of mixture component.
#'
#' @export
calcPostProb <- function(x, mu, sigma, prob, wgts = NULL){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  n <- length(x)
  K <- length(mu)
  if(is.null(wgts)) wgts <- rep(1, n)  
  postp <- matrix(0, nrow = n, ncol = K)
  ll <- 0
  
  for( i in 1:n ) {
    if(is.matrix(prob)){ logpostp <- log(prob[i,]) + dnorm(x[i], mu, sigma, log =TRUE)
    }else{ logpostp <- log(prob) + dnorm(x[i], mu, sigma, log =TRUE) }
    logrowsum <- log_sum_exp(logpostp)
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
#' @param alpha_jackchinook vector of linear predictors for jack Chinook as a proportion of all Chinook.
#' @param pAdultChinook proportion of adult Chinook that are small or large. Defaults to fixed = 1 if adults are just one group.
#' @param Xprop is a list of design matrices with ncol = `length(alpha)`, and the final one ncol = `length(alpha_jackchinook)`.
#' @param K0 Number of non-Chinook species
#' @param K Number of total categories (species + Chinook types).
#'
#' @return list of posterior probabilities and log likelihood of mixture component.
#'
#' @export
calcProportions <- function(alpha, alpha_jackchinook, p_adultchinook, Xprop, K0, K){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  ## Set up proportions. alpha parameter for predicting proportions. 
  ## Xalpha is a list of design matrices for each species.
  np <- nrow(Xprop[[1]])
  logit_p <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
  indx0 <- 1
  for( i in 1:K0 ){
    nc <- ncol(Xprop[[i]]) 
    indx1 <- indx0 + nc - 1
    logit_p[,i] <- as.matrix(Xprop[[i]]) %*% alpha[indx0:indx1]
    indx0 <- indx1 + 1
  }

  logit_pjack <- as.matrix(Xprop[[K0+1]]) %*% alpha_jackchinook
  p_jack <- 1/(1+exp(-logit_pjack))

  p <- matrix(0, nrow = np, ncol = K)
  for( i in 1:np ) {
    p[i, 1:(K0+1)] <- expitM(logit_p[i,])  
    p[i, (K0+1):K ] <- p[i, (K0+1)]*c(p_jack[i], (1-p_jack[i])*p_adultchinook)
  }
  return(p)
}


#' Lognormal distribution (mean and sd on real scale)
#'
#' @param x value on the real scale
#' @param mean_real Mean on real scale of the log normal distribution.
#' @param sd_real Standard deviation on real scale of the log normal distribution.
#' @param log Logical, return log density (TRUE) or on real scale (FALSE).
#'
#' @return Density of log normal distribution for value x.
#'
#' @export
dlognorm_real <- function(x, mean_real, sd_real, log = FALSE){
  sigma <- log(1+(sd_real/mean_real)^2)
  mu <- log(mean_real^2/sqrt(mean_real + sd_real^2))
  if(log) dnorm(log(x), mu, sigma, TRUE) - log(x)
  else dnorm(log(x), mu, sigma, FALSE)/x
}

#' Lognormal distribution (mean and sd on log scale)
#'
#'
#' @param x value on the real scale.
#' @param mean Mean on log scale of the log normal distribution.
#' @param sd Standard deviation on log scale of the log normal distribution.
#' @param log Logical, return log density (TRUE) or on real scale (FALSE).
#'
#' @return Density of log normal distribution for value x.
#'
#' @export
dlognorm <- function(x, mean, sd, log = FALSE){
  if(log) dnorm(log(x), mean, sd, TRUE) - log(x)
  else dnorm(log(x), mean, sd, FALSE)/x
}

#' Binomial distribution allowing for continuous size.
#'
#'
#' @param x Number of successes
#' @param size Number of trials (can be continuous).
#' @param prob Probability of success.
#' @param log Logical, return log density (TRUE) or on real scale (FALSE).
#'
#' @return Density of binomial distribution of x.
#'
#' @export
dbinomc <- function(x, size, prob, log = FALSE){
  logp <- lgamma(size+1) - lgamma(x+1) - lgamma(size - x + 1) + x*log(prob) + (size-x)*log(1-prob)
  if(log) return(logp)
  else return(exp(logp))
}