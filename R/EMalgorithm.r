## EM Algorithm Solver:
library(RTMB)
library(mixtools)
source("utils.r")
source("distributions.r")

predLogitProb <- function(alpha, X){
  
}

EMstep <- function(pars_outer){
  getAll(pars_outer, pars_fixed, dataList, warn = FALSE)
  ll <- 0
  
  ## Parameter Processing
  Kchin <- length(muChin)
  K0 <- length(mu)
  K <- Kchin + K0
  mu <- c(mu, muChin)

  ## Standard deviation incld. observation error
  logsigma_ <- c(logsigma, logsigmaChin)
  sigma <- sqrt(exp(2*logsigma_) + exp(2*logsigma0))

  ## Set up proportions. alpha parameter for predicting proportions. 
  ## Xalpha is a list of design matrices for each species.
  np <- nrow(Xalpha)
  logitp <- matrix(0, nrow = np, ncol = K0+1) ## +1 is Chinook.
  for( i in 1:K0 ){
    logitp[,i] <- Xalpha[,alphaindex[i,1]:alphaindex[i,2]] %*% alpha[alphaindex[i,1]:alphaindex[i,2]]
  }
  logitpChin <- matrix(0, nrow = np, ncol = Kchin)
  ## Default XalphaChin should be intercept.
  for( i in 1:(Kchin-1) ){
    logitpChin[,i] <- XalphaChin[,alphaChinindex[i,1]:alphaChinindex[i,2]] %*% alphaChin[alphaChinindex[i,1]:alphaChinindex[i,2]]
  }
  p <- matrix(0, nrow = np, ncol = K)
  for( i in 1:np ) {
    p[i, 1:(K0+1)] <- expitM(logitp[i,])  
    p[i, (K0+1):K ] <- expitM(logitpChin[i,])*p[i, (K0+1)]
  }
  pobs <- p[grp,]

  ## Adjust lengths for beam spreading:
  ## X is a design matrix that is either distance from shore, or bin id. 
  ## It is expected to contain an intercept term.
  Ladj <- L + beta %*% t(X)

  ## Prior or penalty terms:
  ll <- ll + prior_beta(beta) + prior_sigma0(exp(logsigma0)) + prior_sigma(exp(logsigma)) + prior_mu(mu) ## Add prior to q.
  
  ## Calculate Posterior Probabilities + log likelihood:
  outerLL <- calcPostProb(x = Ladj, mu = mu, sigma = sigma, prob = pobs, wgts = wgts)
  postProb <- outerLL$postp
  ll <- ll + outerLL$ll

  ## Test fishery component:
  if(includeTF){
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


  ## Inner Objective Function
  ## ------------------------------------------
  inner_objective_fn <- function(pars_inner){
    getAll(pars_inner, pars_fixed, warn = FALSE)

    ## Parameter Processing
    # beta <- exp(logbeta)  ## Impact of beam spreading.
    Kchin <- length(muChin)
    K0 <- length(mu)
    K <- Kchin + K0
    mu <- c(mu, muChin)

    ## Standard deviation incld. observation error
    logsigma_ <- c(logsigma, logsigmaChin)
    sigma <- sqrt(exp(2*logsigma_) + exp(2*logsigma0))

    ## Set up proportions
    np <- nrow(logitp)
    p <- matrix(0, nrow = np, ncol = K)
    ## This could be Jack/Adult or Jack/Small Adult/Large Adult. 
    pChin <- expitM(logitpChin)           ## Cooked in assumption that behaviour is same for J + A. pJack < pAdult (can set this via prior).
    for( i in 1:np ) {
      p[i,1:(K0+1)] <- expitM(logitp[i,])  
      p[i, (K0+1):K ] <- pChin*p[i, (K0+1)]
    }
    pobs <- p[grp,]
    logpobs <- log(pobs)

    ## Objective function for length based mixture model, conditional on posterior probs.
    objval <- 0
    for( k in 1:K ){
      objval <- objval - sum(wgts*postProb[,k]*(logpobs[,k] + dnorm(x, mu[k], sigma[k], log = TRUE)))
    }
    objval
  }
  start <- (lapply(pars_outer, RTMB:::getValues))
  F <- MakeTape(inner_objective_fn, start)
  Newton <- F$newton(1:length(F$par()), maxit = 1000 )
  pars_opt <- Newton(numeric(0))
  c(pars_opt, ll)
}


x <- rnormmix(100, c(0.2, 0.1, 0.7), c(0, 5, 10), c(2, 2.5, 3.5)) # simulated data
dataList <- list(x = x, wgts = rep(1, 100))
pars <- list(logitp = logitM(c(0.2, 0.1, 0.7)), mu =  c(0, 5, 10), logsigma = log(c(2, 2, 2)))
pars_fixed <- list()

EM <- MakeTape(EMstep, pars)
## Optimize the tape
EM$simplify()
## Iterate to find fixed point
npar <- length(unlist(pars))
start <- EM$par()
vals <- EM(start)
lli <- -Inf
lli <- c(lli, vals[npar+1])

maxit <- 10000
maxdiff <- 1e-8
iter <- 2
converged <- FALSE
while(iter < maxit & !converged){
  vals <- EM(vals[1:npar])
  lli <- c(lli, vals[npar+1])
  converged <- (lli[iter] - lli[iter-1]) < maxdiff
  iter <- iter + 1
}
cat("number of iterations=", iter, "\n")

pars_fit <- reList(pars, vals)
expitM(pars_fit$logitp)
pars_fit$mu
vals[npar+1]
fitEM <- normalmixEM(x, k = 3)
fitEM$lambda
fitEM$mu
fitEM$loglik
expitM(pars_fit$logitp)
pars_fit$mu

obj <- MakeADFun(negLL, pars)
obj$fn(vals[1:(npar)])
fitEM$loglik

pars_mle <- list(logitp = logitM(fitEM$lambda), mu = fitEM$mu, logsigma = log(fitEM$sigma))
negLL(pars_mle)

n <- length(x)
k <- 3

lambda <- expitM(pars_fit$logitp)
mu <- pars_fit$mu
sigma <- exp(pars_fit$logsigma)

lambda <- fitEM$lambda
mu <- fitEM$mu
sigma <- fitEM$sigma

z <- .C(mixtools:::C_normpost, as.integer(n), as.integer(k),
        as.double(x), as.double(mu), 
        as.double(sigma), as.double(lambda),
        res2 = double(n*k), double(3*k), post = double(n*k),
        loglik = double(1), PACKAGE = "mixtools")
z$loglik
pars_mle <- list(logitp = logitM(fitEM$lambda), mu = fitEM$mu, logsigma = log(fitEM$sigma))
negLL(pars_mle)

ll <- 0
logp <- log(lambda)
for( i in 1:n ){
    logprob <- logp + dnorm(x[i], mu, sigma, log = TRUE)
    maxlp <- max(logprob)
    ll <- ll + (log(sum(exp(logprob - maxlp))) + maxlp)
}
ll
