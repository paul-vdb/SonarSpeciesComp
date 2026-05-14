library(R6)
library(RTMB)

speciesCompModel <- R6Class("model",
  public = list(
    dLogPrior = NULL,
    species_info = NULL,
    test_fishery_counts = NULL
    sonar_lengths = NULL, 
    params_fixed <- NULL,
    params_init <- NULL,
    default_parameters <- NULL,
    initialize = function() {
      self$params_fixed <- list()
      self$params_init <- list()
      defaultPrior(self)
      defaultParams(self)
    },
    setData = function(sonar_lengths, test_fishery_counts, species) {
      self$test_fishery_counts <- test_fishery_counts
      self$sonar_lengths <- sonar_lengths
      self$species_info <- list()
      self$species_info$species <- species
      self$species_info$nspp <- length(species)
    },
    setParameters = function(fixed = list(), inits = list(), delta_limits = c(-1, 1)){
      setParameters(self, inits, fixed, delta_limits)
    }
    model = function(fit = TRUE) {
    },
    setPriors = function(fnlist){
      fnbeta <- extractValue(fnlist[["beta"]], function(x){0})
      fnsd <- extractValue(fnlist[["sd"]], function(x){0})
      include_jac <- any(names(fnlist) == "sd")
      self$prior <- function(x){sum(fnbeta(x$beta)) + sum(fnsd(exp(x$log_sd)) + include_jac*x$log_sd)}
    }
  )
)

defaultParams <- function(self){
  self$default_parameters <- list()
  self$default_parameters$default_q <- c("pink" = 0.25, "sockeye" = 0.5, "jackchinook" = 0.1)
  self$default_parameters$alpha <- c(0,0,0,0)
  self$default_parameters$alpha_jackchinook <- 0
  self$default_parameters$beta <- c(0,0)
  self$default_parameters$mu <- c("smallresident" = 23, "largeresident" = 35, "pink" = 50, "sockeye" = 58, "jackchinook" = 40, "smalladultchinook" = 62, "largeadultchinook" = 80, "adultchinook" = 70)
  self$default_parameters$sigma <- log(c("smallresident" = 2.5, "largeresident" = 2.5, "pink" = 4, "sockeye" = 3.5, "jackchinook" = 2.5, "smalladultchinook" = 4, "largeadultchinook" = 6, "adultchinook" = 7))
  self$default_parameters$sigma0 <- log(4.5)
  self$default_parameters$beta <- c(0,0)
  self$default_parameters$delta_mu <- c("smallresident" = 0, "largeresident" = 0, "pink" = 0, "sockeye" = 0, "jackchinook" = 0, "smalladultchinook" = 0, "largeadultchinook" = 0, "adultchinook" = 0)
}

defaultPrior <- function(self){
  self$dLogPrior <- function(...){0}
}

setParameters = function(self, inits, fixed, delta_limits){
  self$params_fixed <- list()
  self$params_init <- list()
  
  ## Set log catchability:
  out <- extractParams(fixed$q, inits$q, self$default_parameters$q)
  if(!is.null(out$fixed)) self$params_fixed$log_q <- log(out$fixed)
  if(!is.null(out$init)) self$params_init$log_q <- log(out$init)
  
  ## Set alpha
  out <- extractParams(fixed$alpha, inits$alpha, self$default_parameters$alpha)
  self$params_fixed$alpha <- out$fixed
  self$params_init$alpha <- out$init

  ## Set alpha jack chinook
  out <- extractParams(fixed$alpha_jackchinook, inits$alpha_jackchinook, self$default_parameters$alpha_jackchinook)
  self$params_fixed$alpha_jackchinook <- out$fixed
  self$params_init$alpha_jackchinook <- out$init
  
  ## Set beta
  out <- extractParams(fixed$beta, inits$beta, self$default_parameters$beta)
  self$params_fixed$beta <- out$fixed
  self$params_init$beta <- out$init
  
  ## Set mu
  out <- extractParams(fixed$beta, inits$beta, self$default_parameters$beta)
  self$params_fixed$beta <- out$fixed
  self$params_init$beta <- out$init

  ## Set sigma
  out <- extractParams(fixed$sigma, inits$sigma, self$default_parameters$sigma)
  self$params_fixed$log_sigma <- out$fixed
  self$params_init$log_sigma <- out$init
  
  ## Set delta_mu
  out <- extractParams(fixed$delta_mu, inits$delta_mu, self$default_parameters$delta_mu)
  if(!is.null(out$fixed)) self$params_fixed$logit_delta_mu <- logitInterval(out$fixed, delta_limits[1], delta_limits[2])
  if(!is.null(out$init)) self$params_init$logit_delta_mu <- logitInterval(out$init, delta_limits[1], delta_limits[2])
  
  ## Set sigma0
  out <- extractParams(fixed$sigma0, inits$sigma0, self$default_parameters$sigma0)
  if(!is.null(out$fixed)) self$params_fixed$log_sigma0 <- log(out$fixed)
  if(!is.null(out$init)) self$params_init$log_sigma0 <- log(out$init)
}

## ***
## Make function to do fish prediction to pair with test fishery.
## ***

buildEM <- function(self, include_test_fishery = TRUE, adjust_lengths = TRUE){
  
  ## *** data input
  species <- self$species
  species_tf <- colnames(self$testFisheryCounts)
  nspp <- length(species)
  nchin <- grep("chinook", species, values = TRUE)
  nspp_0 <- nspp - nchin

  adjust_lengths <- TRUE
  beta_names <- c("Bin1", "Bin2", "Bin3")
  wgts <- dataList$lengthData$weights
  obs_lengths <- dataList$lengthData$L.cm.adj
  # lower_delta_mu <- -2
  # upper_delta_mu <- 2
  
  EMstep <- function(parsOuter){
    getAll(parsOuter, parsFixed, warn = FALSE)
    ll <- 0
    
    ## Parameter Processing
    ## 1) Mean Fish Size:
    mu <- joinPars(parsOuter$mu, parsFixed$mu, species)
    ## Let the user fit a little shift in the fish size. This will default to 0.
    logit_delta_mu <- joinPars(parsOuter$logit_delta_mu, parsFixed$logit_delta_mu, species) 
    delta_mu <- ilogit_interval(logit_delta_mu, lower_delta_mu, upper_delta_mu)
    mu <- mu + delta_mu

    ## 2) Length standard deviation
    log_sigma <- joinPars(parsOuter$log_sigma, parsFixed$log_sigma, species)  ## *** data input
    ## length measurement error:
    log_sigma0 <- joinPars(parsOuter$log_sigma0, parsFixed$log_sigma0)
    ## Observed standard deviation
    sigma <- sqrt(exp(2*log_sigma) + exp(2*log_sigma0))

    ## 3) Catchability
    log_q <- joinPars(parsOuter$log_q, parsFixed$log_q, species_tf)  ## *** data input
    q <- c(exp(log_q), 1) ## Relative to the last species (chinook).
    nq <- length(q)
    
    ## 4) Proportions relationships
    alpha <- parsOuter$alpha
    alpha_jackchinook <- parsOuter$alpha_jackchinook
    
    ## 5) Beam Spreading Corrections:
    beta <- joinPars(parsOuter$beta, parsFixed$beta, beta_names)
    
    ## Set up proportions. alpha parameter for predicting proportions. 
    ## Xalpha is a list of design matrices for each species.
    p <- calcProportions(alpha = alpha, alpha_jackchinook = alpha_jackchinook, 
      p_adultchinook = p_adultchinook, dataList = dataList, K0 = K0, K = K)

    ## Adjust lengths for beam spreading:
    ## X is a design matrix that is either distance from shore, or bin id. 
    ## It is expected to contain an intercept term.
    if(adjust_lengths){
      L <- obs_lengths - as.matrix(dataList$Xlength) %*% beta
    }else{
      L <- obs_lengths
    }
    
    ## Prior or penalty terms:
    ll <- ll + dprior(sigma = sigma, sigma0 = exp(logsigma0), mu = mu,
                      beta = beta, catchability = q[-nq], alpha = alpha, 
                      alpha_jackchinook = alpha_jackchinook)
                      
    ## Calculate Posterior Probabilities + log likelihood:
    p_obs <- p[dataList$lengthData$grpIndex,]
    outer_ll <- calcPostProb(x = L, mu = muc, sigma = sigma, prob = pobs, wgts = wgts)
    post_prob <- outer_ll$postp
    ll <- ll + outer_ll$ll

    ## Test fishery component:
    if(include_test_fishery){
      ndays <- nrow(dataList$testFisheryCounts)
      padultchinook <- 1-rowSums(p[,1:(K0+1)])  ## Proportion combined adult chinook.
      nq <- ncol(dataList$testFisheryCounts)
      for( d in 1:ndays ){
        idx <- which(dataList$predDF$day == d)
        Nd <- sapply(dataList$qSppIndices, FUN = function(x){sum(p[idx,x]*dataList$predDF$SalmonCountTF[idx])})
        Ndc <- sum(padultchinook[idx]*dataList$predDF$SalmonCountTF[idx])
        Nd <- c(Nd, Ndc)
        prob <- (q*Nd)/sum(q*Nd)
        ll <- ll + dataList$testFisheryWeights[d]*dmultinom(dataList$testFisheryCounts[d,], prob = prob, size = sum(dataList$testFisheryCounts[d,]), log = TRUE)
      }
    }

    ## Inner Objective Function
    ## ------------------------------------------
    inner_objective_fn <- function(parsInner){
      getAll(parsInner, parsFixed, warn = FALSE)

      # Parameter Processing
      Kchin <- length(muChin)
      K0 <- length(mu)
      K <- Kchin + K0
      q <- exp(c(logq, 0))

      ## Allow estimation of small changes from the mean to correct for 'near' shore.
      muDelta <- dataList$dmuLimit[,1] + (dataList$dmuLimit[,2]-dataList$dmuLimit[,1])/(1+exp(-logitdmu))
      mu <- mu + muDelta
      muDeltaChin <- dataList$dmuChinLimit[,1] + (dataList$dmuChinLimit[,2]-dataList$dmuChinLimit[,1])/(1+exp(-logitdmuChin))
      muChin <- muChin + muDeltaChin
      muc <- c(mu, muChin)

      ## Standard deviation incld. observation error
      logsigma_ <- c(logsigma, logsigmaChin)
      sigma <- sqrt(exp(2*logsigma_) + exp(2*logsigma0))

      ## Set up proportions. alpha parameter for predicting proportions. 
      ## Xalpha is a list of design matrices for each species.
      p <- calcProportions(alpha = alpha, alphaJackChinook = alphaJackChinook, 
        pAdultChinook = pAdultChinook, dataList = dataList, K0 = K0, K = K)
      
      # np <- nrow(dataList$predDF)
      # logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
      # indx0 <- 1
      # for( i in 1:K0 ){
        # nc <- ncol(dataList$Xprop[[i]]) 
        # indx1 <- indx0 + nc - 1
        # logitp[,i] <- as.matrix(dataList$Xprop[[i]]) %*% alpha[indx0:indx1]
        # indx0 <- indx1 + 1
      # }

      # logitpjack <- as.matrix(dataList$XpropChin) %*% alphaJackChinook
      # pjack <- 1/(1+exp(-logitpjack))

      # p <- matrix(0, nrow = np, ncol = K)
      # for( i in 1:np ) {
        # p[i, 1:(K0+1)] <- expitM(logitp[i,])  
        # p[i, (K0+1):K ] <- p[i, (K0+1)]*c(pjack[i], (1-pjack[i])*pAdultChinook)
      # }
      logpobs <- log(p[dataList$lengthData$grpIndex,])

      ## Objective function for length based mixture model, conditional on posterior probs.
      objval <- 0
      ## Prior jack chinook based on overal proportion.
      # objval <- objval - dbeta(Njc/(Nkc+Njc), 10, 500, log = TRUE) 

      ## Adjust lengths for beam spreading:
      ## X is a design matrix that is either distance from shore, or bin id. 
      ## It is expected to contain an intercept term.
      if(adjustLengths){
        L <- dataList$lengthData$L.cm.adj - as.matrix(dataList$Xlength) %*% beta
      }else{
        L <- dataList$lengthData$L.cm.adj
      }
      wgts <- dataList$lengthData$weights
      
      ## Prior/penalty terms: Probably should do some sort of transformation:
      objval <- objval - dprior(sigma = exp(logsigma), sigmaChin = exp(logsigmaChin), sigma0 = exp(logsigma0), mu = mu, muChin = muChin, 
                      beta = beta, q = exp(logq), alpha, alphaJackChinook)
      for( k in 1:K ){
        objval <- objval - sum(wgts*postProb[,k]*(logpobs[,k] + dnorm(L, muc[k], sigma[k], log = TRUE)))
      }
      ## Test fishery component:
      if(includeTestFishery){
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
          objval <- objval - dataList$testFisheryWeights[d]*dmultinom(dataList$testFisheryCounts[d,], prob = prob, size = sum(dataList$testFisheryCounts[d,]), log = TRUE)
        }
      }
      objval
    }
    start <- (lapply(parsOuter, RTMB:::getValues))
    F <- MakeTape(inner_objective_fn, start)
    Newton <- F$newton(1:length(F$par()), maxit = 1000)
    pars_opt <- Newton(numeric(0))
    c(pars_opt, ll)
  }
    self$EM <- MakeTape(EMstep, parsInit)
}

## Run the EM algorithm:
runEM <- function(EM, control){
  start <- EM$par()
  npar <- length(start)
  vals <- EM(start)
  if(any(is.nan(vals))) 
    stop("EM Algorithm finding NaN values. Either provide better initial values or more likely choose a different formulation. Occurs often with Test Fishery and Lengths don't agree.")
  lli <- -Inf
  lli <- c(lli, vals[npar+1])

  maxit <- extractControls( control$maxiters, 1000 )
  tol <- extractControls( control$tolerance, 1e-8 )
  relativeDiff <- extractControls( control$relativeDifference, FALSE )
  verbose <- extractControls( control$verbose, FALSE )

  if(verbose) cat("Running EM Algorithm...\n")
  if(verbose) pb <- txtProgressBar(min = 1, max = maxit, initial = 2) 
  iter <- 2
  converged <- FALSE
  diff <- Inf
  vals.prev <- vals
  while(iter < maxit & !converged){
    vals <- tryCatch(EM(vals[1:npar]), warning = function(w){c(NaN, NaN, NaN)})
    if(!is.nan(vals[1])) {
      vals.prev <- vals
    }
    if(is.nan(vals[1])){
      cat("[Warning]  Issues with initial values. Did not converge, returning the last value that successfully fitted.\n")
      vals <- vals.prev
      break;
    }
    lli <- c(lli, vals[npar+1])   
    if(relativeDiff & iter > 2){ diff <- abs((lli[iter] - lli[iter-1])/lli[iter-1])
    }else{ diff <- abs(lli[iter] - lli[iter-1]) }
    converged <- diff < tol
    iter <- iter + 1
    if(verbose) setTxtProgressBar(pb,iter)    
  }
  if(verbose) close(pb)
  if(verbose){
    cat("number of iterations =", iter, "\n")
    if(relativeDiff) cat("Convergence was evaluated based on the relative difference of the log likelihood:", diff, ".\n") 
    else cat("Convergence was evaluated based on the difference of the log likelihood:", diff, ".\n") 
  }
  if(!converged) {
    if(maxit < iter) cat("[Warning]  Maximum iterations of", maxit, "were reached without convergence to a difference of", tol, "between the log likelihood.\n")
    if(maxit >= iter) cat("[Warning]  Failed to converge with a tolerance of", tol, "between the log likelihood.\n")
    if(relativeDiff) cat("If a relative difference of", diff, "seems to be close enough then you may consider ignoring this message.\n")
    else cat("If a difference of", diff, "seems to be close enough then you may consider ignoring this message.\n")
  }
  return(vals)
}

basicMixtureModel <- function(x, K, mu, sigma, prob, mu_fixed = NULL, sigma_fixed = NULL, wgts = NULL, names = NULL, control = list()){

  if(missing(K) & missing(names)) stop("Provide number of components, K")
  if(missing(K)) K <- length(names)
  if(is.null(names)) names <- 1:K

  if(missing(prob)) prob <- rep(1/K, K)
  if(missing(mu)){
    mu <- rnorm(K, mean(x), sd(x))
    names(mu) <- names
    mu <- mu[!names %in% names(mu_fixed)]
    if(length(mu) == 0) mu <- NULL
  }
  if(missing(sigma)){
    sigma <- rep(sd(x)/K, K)
    names(sigma) <- names
    sigma <- sigma[!names %in% names(sigma_fixed)]
    if(length(sigma) == 0) sigma <- NULL
  }
  if(is.null(wgts)) wgts <- rep(1, length(x))
  logit_prob <- fn_null("logitM", prob)
  
  pars_init <- list()
  if(!is.null(mu)) pars_init$mu <- mu
  if(!is.null(sigma)) pars_init$log_sigma <- log_null(sigma)
  pars_init$logit_prob <- logit_prob
  pars_fixed <- list(mu_fixed = mu_fixed, log_sigma_fixed = log_null(sigma_fixed))
  
  EMstep <- function(pars_outer){
    ll <- 0
    n <- length(x)
    mu <- joinPars(pars_outer$mu, pars_fixed$mu, names)
    log_sigma <- joinPars(pars_outer$log_sigma, pars_fixed$log_sigma, names)  ## *** data input
    sigma <- exp(log_sigma)
    ngrp <- length(mu)
    
    logit_prob <- pars_outer$logit_prob
    pmix <- expitM(logit_prob)
    log_pmix <- log(pmix)
    
    post_probs <- matrix(0, nrow = n, ncol = ngrp)
    for( i in 1:n ){
      logfi <- log_pmix + dnorm(x[i], mu, sigma, log = TRUE)
      maxfi <- max(logfi)
      logfisum <- log(sum(exp(logfi - maxfi))) + maxfi
      post_probs[i,] <- exp(logfi - logfisum)
      ll <- ll + logfisum*wgts[i]
    }
    
    ## Inner Objective Function
    ## ------------------------------------------
    innerObjectiveFn <- function(pars_inner){

      mu <- joinPars(pars_inner$mu, pars_fixed$mu, names)
      log_sigma <- joinPars(pars_inner$log_sigma, pars_fixed$log_sigma, names)  ## *** data input
      sigma <- exp(log_sigma)
      ngrp <- length(mu)
      
      logit_prob <- pars_inner$logit_prob
      pmix <- expitM(logit_prob)
      log_pmix <- log(pmix)
      
      ## Objective function for mixture model, conditional on posterior probs.
      objval <- 0
      for( k in 1:ngrp ){
        objval <- objval - sum(wgts*post_probs[,k]*(log_pmix[k] + dnorm(x, mu[k], sigma[k], log = TRUE)))
      }
      objval
    }
    start <- (lapply(pars_outer, RTMB:::getValues))
    F <- MakeTape(innerObjectiveFn, start)
    Newton <- F$newton(1:length(F$par()), maxit = 1000)
    pars_opt <- Newton(numeric(0))
    c(pars_opt, ll)
  }
  EM <- MakeTape(EMstep, pars_init)
  vals <- runEM(EM, control)
  pars_est <- reList(pars_init, vals)
  estimates <- list(mu = joinPars(pars_est$mu, pars_fixed$mu, names))
  estimates$sigma <- joinPars(fn_null('exp', pars_est$log_sigma), fn_null('exp', pars_fixed$log_sigma), names)
  estimates$proportion <- joinPars(fn_null('expitM', pars_est$logit_prob), fn_null('expitM', pars_fixed$logit_prob), names)
  estimates
}


est <- NULL
est2 <- NULL
for( i in 1:500 ){
  N <- c(2000, 1000)
  nsmp <- 25
  p1 <- c(0.3, 0.4, 0.3)
  p2 <- c(0.15, 0.65, 0.2)
  ptrue <- (p1*N[1] + p2*N[2])/sum(N)
  id1 <- sample(3, nsmp, prob = p1, replace = TRUE)
  id2 <- sample(3, nsmp, prob = p2, replace = TRUE)
  x1 <- rnorm(nsmp, c(25, 45, 65)[id1], c(3, 3.7, 2.5)[id1])
  x2 <- rnorm(nsmp, c(25, 45, 65)[id2], c(3, 3.7, 2.5)[id2])
  wgts <- c(rep(N[1]/nsmp, nsmp), rep(N[2]/nsmp, nsmp))
  fit <- basicMixtureModel(x=c(x1, x2), K = 3, prob = c(0.2, 0.4, 0.4), mu = c(20, 40, 60), wgts = wgts)
  fit2 <- basicMixtureModel(x=c(x1, x2), K = 3, prob = c(0.2, 0.4, 0.4), mu = c(20, 40, 60))
  est <- rbind(est, fit$proportion[order(fit$mu)] - ptrue)
  est2 <- rbind(est2, fit2$proportion[order(fit$mu)] - ptrue)
}
boxplot(est, ylim = c(-1, 1))
abline(h = 0, col='red')
boxplot(est2, ylim = c(-1, 1))
abline(h = 0, col = 'red')



#' Run EM Algorithm to fit Joint Species Composition Model
#'
#' Loop through steps of the EM algorithm using the species composition R6 object.
#'
#' @param speciesComp R6 object that contain all the of data and parameters to fit the joint species composition model.
#' @param simulatedData logical (default = FALSE) whether or not to fit to simulated data that is held withing the speciesComp object.
#'
#' @details Run the EM algorithm to fit the joint species composition model. Save the estimated parameters and estimate the daily proportions of each stock
#' and save those within the `speciesComp` R6 object and returned to the user. See \code{speciesComp$estimatedHourlyProportions} to access the stratum hour level 
#' estimated species proportions, \code{speciesComp$estimatedDailyProportions} to access the aggregated estimate of proportions of each stock for the entire day, 
#' within each stratum, and \code{speciesComp$estimatedParameters} to access the estimated and fixed parameters.
#'
#'
#' @export
runEMAlgorithm <- function(speciesComp, simulatedData = FALSE){
  parsInit <- speciesComp$parsInit
  if(!simulatedData){
    dataenv <- local({dataList <- speciesComp$analysisData; includeTestFishery <- speciesComp$includeTestFishery; 
                      parsFixed <- speciesComp$parsFixed; adjustLengths <- speciesComp$adjustLengths; dprior = speciesComp$priorDist;
                      environment()})
  }else{
    dataenv <- local({dataList <- speciesComp$simData; includeTestFishery <- speciesComp$includeTestFishery; 
                    parsFixed <- speciesComp$parsFixed; adjustLengths <- speciesComp$adjustLengths; dprior = speciesComp$priorDist;
                    environment()})
  }
  environment(EMstep) <- dataenv
  EM <- MakeTape(EMstep, parsInit)

  ## Optimize the tape
  EM$simplify()
  ## Iterate to find fixed point
  npar <- length(unlist(speciesComp$parsInit))
  start <- EM$par()
  vals <- EM(start)
  if(any(is.nan(vals))) 
    stop("EM Algorithm finding NaN values. Either provide better initial values or more likely choose a different formulation. Occurs often with Test Fishery and Lengths don't agree.")
  lli <- -Inf
  lli <- c(lli, vals[npar+1])

  maxit <- extractControls( speciesComp$optimControl$maxiters, 1000 )
  tol <- extractControls( speciesComp$optimControl$tolerance, 1e-8 )
  relativeDiff <- extractControls( speciesComp$optimControl$relativeDifference, FALSE )
  verbose <- extractControls( speciesComp$optimControl$verbose, TRUE )

  if(verbose) cat("Running EM Algorithm\n")
  if(verbose) pb <- txtProgressBar(min = 1, max = maxit, initial = 2) 
  iter <- 2
  converged <- FALSE
  diff <- Inf
  vals.prev <- vals
  while(iter < maxit & !converged){
    vals <- tryCatch(EM(vals[1:npar]), warning = function(w){c(NaN, NaN, NaN)})
    if(!is.nan(vals[1])) {
      vals.prev <- vals
    }
    if(is.nan(vals[1])){
      cat("[Warning]  Issues with initial values. Did not converge, returning the last value that successfully fitted.\n")
      vals <- vals.prev
      break;
    }
    lli <- c(lli, vals[npar+1])   
    if(relativeDiff & iter > 2){ diff <- abs((lli[iter] - lli[iter-1])/lli[iter-1])
    }else{ diff <- abs(lli[iter] - lli[iter-1]) }
    converged <- diff < tol
    iter <- iter + 1
    if(verbose) setTxtProgressBar(pb,iter)
  }
  if(verbose) close(pb)
  if(verbose){
    cat("number of iterations =", iter, "\n")
    if(relativeDiff) cat("Convergence was evaluated based on the relative difference of the log likelihood:", diff, ".\n") 
    else cat("Convergence was evaluated based on the difference of the log likelihood:", diff, ".\n") 
  }
  if(!converged) {
    if(maxit < iter) cat("[Warning]  Maximum iterations of", maxit, "were reached without convergence to a difference of", tol, "between the log likelihood.\n")
    if(maxit >= iter) cat("[Warning]  Failed to converge with a tolerance of", tol, "between the log likelihood.\n")
    if(relativeDiff) cat("If a relative difference of", diff, "seems to be close enough then you may consider ignoring this message.\n")
    else cat("If a difference of", diff, "seems to be close enough then you may consider ignoring this message.\n")
  }
  
  parsFit <- reList(parsInit, vals)
  parsFit <- c(parsFit, speciesComp$parsFixed)
  speciesComp$estimatedParameters <- parsFit
  
  ## Set up proportions. alpha parameter for predicting proportions. 
  ## Xalpha is a list of design matrices for each species.
  np <- nrow(speciesComp$analysisData$predDF)
  K0 <- length(parsFit$mu)
  Kchin <- length(parsFit$muChin)
  K <- K0 + Kchin
  logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
  indx0 <- 1
  for( i in 1:length(parsFit$mu) ){
    nc <- ncol(speciesComp$analysisData$Xprop[[i]]) 
    indx1 <- indx0 + nc - 1
    logitp[,i] <- as.matrix(speciesComp$analysisData$Xprop[[i]]) %*% parsFit$alpha[indx0:indx1]
    indx0 <- indx1 + 1
  }
  logitpjack <- as.matrix(speciesComp$analysisData$XpropChin) %*% parsFit$alphaJackChinook
  pjack <- 1/(1+exp(-logitpjack))

  p <- matrix(0, nrow = np, ncol = K)
  N <- matrix(0, nrow = np, ncol = K)
  for( i in 1:np ) {
    p[i, 1:(K0+1)] <- expitM(logitp[i,])  
    p[i, (K0+1):K ] <- p[i, (K0+1)] * c(pjack[i], (1-pjack[i])*parsFit$pAdultChinook) 
    N[i,] <- p[i,] * speciesComp$analysisData$predDF$SalmonCount[i]/speciesComp$analysisData$predDF$MinsCounted[i] # Count or SalmonCount here????? Think about it a bit more...
  }
  ## A lot of work to not have a tidy dependence...
  form <- paste0("cbind(", paste(speciesComp$species, collapse = ","), ") ~ Date + SonarBank + SonarBin + SonarAim + day")
  counts <- speciesComp$analysisData$predDF |> aggregate(cbind(Count = Count/MinsCounted*60, SalmonCount = SalmonCount/MinsCounted*60, NHour = 1) ~ Date + SonarBank + SonarBin + SonarAim + day, sum)
  colnames(N) <- speciesComp$species
  estimatedProp <- cbind(speciesComp$analysisData$predDF, N)
  estimatedProp <- estimatedProp |> aggregate(as.formula(form), sum)
  estimatedProp[, speciesComp$species] <- t(apply(estimatedProp[, speciesComp$species], 1, FUN = function(x){x/sum(x)}))
  estimatedProp <- counts |> merge(estimatedProp)
  speciesComp$estimatedDailyProportions <- estimatedProp
  colnames(p) <- speciesComp$species
  hourlyPredDF <- cbind(speciesComp$analysisData$predDF[, c("Date", "day", "SonarBank", "SonarAim", "SonarBin", "Hour", "HourOrder", "MinsCounted", "Count", "SalmonCount")], p)
  speciesComp$estimatedHourlyProportions <- hourlyPredDF
}

# parsInit <- speciesComp$parsInit
# speciesComp$fitModel()
# dataenv <- local({dataList <- speciesComp$analysisData; includeTestFishery <- speciesComp$includeTestFishery; 
                # parsFixed <- speciesComp$parsFixed; adjustLengths <- speciesComp$adjustLengths; dprior = speciesComp$priorDist;
                # environment()})
# environment(negLogDensity) <- dataenv
# obj <- MakeADFun(negLogDensity, parsInit)
## Find Hessian:
# pars.em <- do.call('c', speciesComp$estimatedParameters[names(speciesComp$parsInit)])
# obj$fn(pars.em)
# sdrep <- sdreport(obj)
# sdrep
## Can ADD ADREPORT to get Standard Error of Estimates or do bootstrapping...

# library(tmbstan)
# fit <- tmbstan(obj, chains = 1, iter = 1000, warmup = 200, init = list(speciesComp$estimatedParameters[names(speciesComp$parsInit)]))

# traceplot(fit, pars = "logq")
# traceplot(fit, pars = "alpha")
# traceplot(fit, pars = "beta")
# fit
# speciesComp$estimatedParameters[names(speciesComp$parsInit)]

# speciesComp$setPriors(priors = list(
  # beta = function(x){sum(dnorm(x, c(-1, -0.8), c(0.1, 0.1), log = TRUE))}, 
  # sigma0 =  function(x){sum(dgamma(x, 1, 0.2, log = TRUE) + log(x))},
  # catchability =  function(x){sum(dnorm(x, c(0.25, 0.1), c(0.3, 0.00001), log = TRUE) + log(x))},
  # alpha =  function(x){sum(dnorm(x, 0, 10, log = TRUE))},  
  # alphaJackChinook = function(x){sum(dnorm(x, 0, 10, log = TRUE))}
  # ))
# speciesComp$fitModel()
