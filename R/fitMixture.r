#' Fit joint species composition model
#'
#' Expectation-Maximization algorithm for fitting the joint hydroacoustic lengths and test fishery model.
#'
#' @param self The R6 `speciesCompModel` that is set up for a model run. 
#'
#' @details Techincal details for this function are provided as a vignette. Computes posterior probabilities of lengths to species
#' given current values of the parameters. Those probabilities are passed to the objective function and then the objective function is maximized to
#' find the new parameters. This process is repeated until the log-likelihood becomes stable and the parameter estimates are then returned along with
#' being being stored in the R6 object $params_estimated.
#'
#' @return List of model estimates (`$params_estimated`), convergence details (`$convergence`), 
#' estimate of daily species passage (`$N_daily`), and daily proportions (`$p_daily`).
#'
#' @export
fit_joint_model <- function(self){
  species <- self$species_info$species
  names_tf <- self$data_info$test_fishery_input
  species_N <- self$species_info$species_predict
  
  nspp <- length(species)
  nchin <- self$species_info$nchinook
  nother <- self$species_info$nother
  p_adultchinook <- 1
  ndays <- self$ndays
  
  adjust_lengths <- self$fit_info$adjust_lengths
  include_test_fishery <- self$fit_info$include_test_fishery

  wgts <- self$data_list$length_data$weights
  obs_lengths <- self$data_list$length_data$L.cm.adj
  lower_delta_mu <- self$data_list$lower_delta_mu
  upper_delta_mu <- self$data_list$upper_delta_mu
  nbeta <- ncol(self$data_list$X_length)
  test_fishery_catch <- self$data_list$test_fishery_catch
  test_fishery_weights <- self$data_info$test_fishery_weights
  
  pars_init <- self$params_init
  pars_fixed <- self$params_fixed
  
  ## Need to ensure that beta is fixed if adjust_lengths = FALSE.
  if(!adjust_lengths){
    pars_fixed$beta <- c(0, 0)
    pars_init$beta <- NULL
  }

  nalpha <- length(c(pars_init$alpha, pars_fixed$alpha))
  nalpha_jackchinook <- length(c(pars_init$alpha_jackchinook, pars_fixed$alpha_jackchinook))
  
  ## pars_outer <- pars_init
  EMstep <- function(pars_outer){
    # Necessary in packages
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")
  
    ll <- 0
    
    ## Parameter Processing
    ## 1) Mean Fish Size:
    mu <- joinPars(pars_outer$mu, pars_fixed$mu, species)
    ## Let the user fit a little shift in the fish size. This will default to 0.
    logit_delta_mu <- joinPars(pars_outer$logit_delta_mu, pars_fixed$logit_delta_mu, species) 
    delta_mu <- ilogitInterval(logit_delta_mu, lower_delta_mu, upper_delta_mu)
    mu <- mu + delta_mu

    ## 2) Length standard deviation
    log_sigma <- joinPars(pars_outer$log_sigma, pars_fixed$log_sigma, species)  ## *** data input
    ## length measurement error:
    log_sigma0 <- joinPars(pars_outer$log_sigma0, pars_fixed$log_sigma0, 1)
    ## Observed standard deviation
    sigma <- sqrt(exp(2*log_sigma) + exp(2*log_sigma0))

    ## 3) Catchability
    log_qinv <- joinPars(pars_outer$log_qinv, pars_fixed$log_qinv, names_tf)  ## *** data input
    qinv <- exp(log_qinv)  ## Relative to the last species (chinook).
    
    ## 4) Proportions relationships
    alpha <- joinPars(pars_outer$alpha, pars_fixed$alpha, 1:nalpha)
    alpha_jackchinook <- joinPars(pars_outer$alpha_jackchinook, pars_fixed$alpha_jackchinook, 1:nalpha_jackchinook)
    
    ## 5) Beam Spreading Corrections:
    beta <- joinPars(pars_outer$beta, pars_fixed$beta, 1:nbeta)
    
    ## Prior Distributions:
    ll <- ll + self$prior_distributions$dlog_qinv(exp(log_qinv)) + self$prior_distributions$dbeta(beta) + self$prior_distributions$dalpha_jackchinook(alpha_jackchinook) +
          self$prior_distributions$dlog_sigma(exp(log_sigma)) + self$prior_distributions$dmu(mu) + self$prior_distributions$dlogit_delta_mu(delta_mu) +
          self$prior_distributions$dlog_sigma0(exp(log_sigma0))

    ll <- ll + self$prior_jacobians$dlog_qinv(exp(log_qinv)) + self$prior_jacobians$dlog_sigma(exp(log_sigma)) + 
               self$prior_jacobians$dlogit_delta_mu(delta_mu, lower = lower_delta_mu, upper = upper_delta_mu) + self$prior_jacobians$dlog_sigma0(exp(log_sigma0))

    ## Set up proportions. alpha parameter for predicting proportions. 
    ## Xalpha is a list of design matrices for each species.
    p <- calcProportions(alpha = alpha, alpha_jackchinook = alpha_jackchinook, 
      p_adultchinook = pars_fixed$proportion_adultchinook, Xprop = self$data_list$X_proportions, 
      K0 = nother, K = nspp)

    ## Adjust lengths for beam spreading:
    ## X is a design matrix that is either distance from shore, or bin id. 
    ## It is expected to contain an intercept term.
    if(adjust_lengths){
      L <- obs_lengths - as.matrix(self$data_list$X_length) %*% beta
    }else{
      L <- obs_lengths
    }
                          
    ## Calculate Posterior Probabilities + log likelihood:
    p_obs <- p[self$data_list$length_data$grpIndex,]
    outer_ll <- calcPostProb(x = L, mu = mu, sigma = sigma, prob = p_obs, wgts = wgts)
    post_prob <- outer_ll$postp
    ll <- ll + outer_ll$ll

    ## Estimate N:
    p_daily <- estimate_daily_proportions(p, pred_df = self$data_list$pred_df, species)
    N_daily <- estimate_daily_salmon(p_daily, total_salmon=self$data_list$total_salmon, species)  ## Combine adult chinook and remove smallresident fish.
    
    if("largeresident" %in% species)  ll <- ll + self$prior_distributions$dlargeresident(N_daily[,1])
    
    ## Test fishery component:
    if(include_test_fishery){
      ## catch[species] ~ dbinom(N[species], q[species]*effort)
      Ntest <- N_daily[cbind(test_fishery_catch$day, test_fishery_catch$N_index)]
      ptest <- 1/qinv[test_fishery_catch$q_index]*test_fishery_catch$effort
      # ll <- ll + test_fishery_weights*sum(dbinomc(x = test_fishery_catch$catch, size = Ntest, prob = ptest, log = TRUE))
      ll <- ll + test_fishery_weights*sum(dpois(x = test_fishery_catch$catch, Ntest*ptest, log = TRUE))
    }
    
    ## Inner Objective Function
    ## ------------------------------------------
    inner_objective_fn <- function(pars_inner){
      # Necessary in packages
      "c" <- ADoverload("c")
      "[<-" <- ADoverload("[<-")

      ## Parameter Processing
      ## 1) Mean Fish Size:
      mu <- joinPars(pars_inner$mu, pars_fixed$mu, species)
      ## Let the user fit a little shift in the fish size. This will default to 0.
      logit_delta_mu <- joinPars(pars_inner$logit_delta_mu, pars_fixed$logit_delta_mu, species) 
      delta_mu <- ilogitInterval(logit_delta_mu, lower_delta_mu, upper_delta_mu)
      mu <- mu + delta_mu

      ## 2) Length standard deviation
      log_sigma <- joinPars(pars_inner$log_sigma, pars_fixed$log_sigma, species)  ## *** data input
      ## length measurement error:
      log_sigma0 <- joinPars(pars_inner$log_sigma0, pars_fixed$log_sigma0, 1)
      ## Observed standard deviation
      sigma <- sqrt(exp(2*log_sigma) + exp(2*log_sigma0))

      ## 3) Catchability
      log_qinv <- joinPars(pars_inner$log_qinv, pars_fixed$log_qinv, names_tf) 
      qinv <- exp(log_qinv)  ## Relative to the last species (chinook).
      
      ## 4) Proportions relationships
      alpha <- joinPars(pars_inner$alpha, pars_fixed$alpha, 1:nalpha)
      alpha_jackchinook <- joinPars(pars_inner$alpha_jackchinook, pars_fixed$alpha_jackchinook, 1:nalpha_jackchinook)
      
      ## 5) Beam Spreading Corrections:
      beta <- joinPars(pars_inner$beta, pars_fixed$beta, 1:nbeta)

      ## Objective function for length based mixture model, conditional on posterior probs.
      objval <- 0
      ## Prior Distributions:
      objval <- objval - self$prior_distributions$dlog_qinv(exp(log_qinv)) - self$prior_distributions$dbeta(beta) - self$prior_distributions$dalpha_jackchinook(alpha_jackchinook) - 
          self$prior_distributions$dlog_sigma(exp(log_sigma)) - self$prior_distributions$dmu(mu) - self$prior_distributions$dlogit_delta_mu(delta_mu) - 
          self$prior_distributions$dlog_sigma0(exp(log_sigma0))
      ## Jacobian transformations for some:
      objval <- objval - self$prior_jacobians$dlog_qinv(exp(log_qinv)) - self$prior_jacobians$dlog_sigma(exp(log_sigma)) - 
                         self$prior_jacobians$dlogit_delta_mu(delta_mu, lower = lower_delta_mu, upper = upper_delta_mu) - self$prior_jacobians$dlog_sigma0(exp(log_sigma0))
            
      ## Set up proportions. alpha parameter for predicting proportions. 
      p <- calcProportions(alpha = alpha, alpha_jackchinook = alpha_jackchinook, 
        p_adultchinook = pars_fixed$proportion_adultchinook, Xprop = self$data_list$X_proportions, 
        K0 = nother, K = nspp)

      ## Adjust lengths for beam spreading:
      if(adjust_lengths){
        L <- obs_lengths - as.matrix(self$data_list$X_length) %*% beta
      }else{
        L <- obs_lengths
      }

      logpobs <- log(p[self$data_list$length_data$grpIndex,])

      ## Prior jack chinook based on overal proportion.
      # objval <- objval - dbeta(Njc/(Nkc+Njc), 10, 500, log = TRUE) 
      
      ## Prior/penalty terms: Probably should do some sort of transformation:
      # objval <- objval - dprior(sigma = exp(logsigma), sigmaChin = exp(logsigmaChin), sigma0 = exp(logsigma0), mu = mu, muChin = muChin, 
                      # beta = beta, q = exp(logq), alpha, alphaJackChinook)
      for( k in 1:nspp ){
        objval <- objval - sum(wgts*post_prob[,k]*(logpobs[,k] + dnorm(L, mu[k], sigma[k], log = TRUE)))
      }
      ## Estimate N:
      p_daily <- estimate_daily_proportions(p, self$data_list$pred_df, species)
      N_daily <- estimate_daily_salmon(p_daily, self$data_list$total_salmon, species)  ## Combine adult chinook and remove smallresident fish.

      if("largeresident" %in% species) objval <- objval - self$prior_distributions$dlargeresident(N_daily[,1])
      
      ## Test fishery component:
      if(include_test_fishery){
        ## catch[species] ~ dbinom(N[species], q[species]*effort)
        Ntest <- N_daily[cbind(test_fishery_catch$day, test_fishery_catch$N_index)]
        ptest <- 1/qinv[test_fishery_catch$q_index]*test_fishery_catch$effort
        # objval <- objval - test_fishery_weights*sum(dbinomc(x = test_fishery_catch$catch, size = Ntest, prob = ptest, log = TRUE))
        objval <- objval - test_fishery_weights*sum(dpois(x = test_fishery_catch$catch, Ntest*ptest, log = TRUE))
      }
      objval
    }
    start <- (lapply(pars_outer, RTMB:::getValues))
    F <- MakeTape(inner_objective_fn, start)
    Newton <- F$newton(1:length(F$par()), maxit = 1000)
    pars_opt <- Newton(numeric(0))
    c(pars_opt, ll)
  }
  EM <- MakeTape(EMstep, pars_init)
  EM$simplify()

  fit <- runEM(EM, control = self$fit_info)  
  
  if(fit$convergence["converged"] > 0){
    cat("For date:", self$est_date, "Model did not successfully converge.\n")
    # cat("Trying to fit a second time.\n")
    # EM(fit$values[-length(fit$values)])
    # fit <- runEM(EM, control = self$fit_info)
  }
  
  pars_est <- reList(pars_init, fit$values)
  estimates <- list(mu = joinPars(pars_est$mu, pars_fixed$mu, species))
  estimates$sigma <- joinPars(fn_null(exp, x = pars_est$log_sigma), fn_null(exp, x = pars_fixed$log_sigma), species)
  estimates$alpha <- joinPars(pars_est$alpha, pars_fixed$alpha, 1:nalpha)
  estimates$alpha_jackchinook <- joinPars(pars_est$alpha_jackchinook, pars_fixed$alpha_jackchinook, 1:nalpha_jackchinook)
  estimates$qinv <- joinPars(fn_null(exp, x = pars_est$log_qinv), fn_null(exp, x = pars_fixed$log_qinv), names_tf)
  estimates$delta_mu <- joinPars(fn_null(ilogitInterval, x=pars_est$logit_delta_mu, lower = lower_delta_mu, upper = upper_delta_mu), fn_null(ilogitInterval, x=pars_fixed$logit_delta_mu, lower = lower_delta_mu, upper = upper_delta_mu), species)
  estimates$mu_adjusted <- estimates$mu + estimates$delta_mu
  estimates$sigma0 <- joinPars(fn_null(exp, x=pars_est$log_sigma0), fn_null(exp, x=pars_fixed$log_sigma0), 1)
  estimates$beta <- joinPars(pars_est$beta, pars_fixed$beta, 1:nbeta)

  ## Set up proportions. alpha parameter for predicting proportions. 
  p <- calcProportions(alpha = pars_est$alpha,alpha_jackchinook =  pars_est$alpha_jackchinook, 
    p_adultchinook = pars_fixed$proportion_adultchinook, Xprop = self$data_list$X_proportions, 
    K0 = nother, K = nspp)
  p_daily <- estimate_daily_proportions(p, self$data_list$pred_df, species)
  N_daily <- estimate_daily_salmon(p_daily, self$data_list$total_salmon, species)  ## Combine adult chinook and remove smallresident fish.
  N_daily <- data.frame(N_daily)
  names(N_daily) <- self$species_info$species_predict
  estimates$N_daily <- cbind(Date = self$data_list$days, N_daily)

  p_daily <- cbind(self$data_list$days, as.data.frame(p_daily))
  names(p_daily) <- c("Date", species)
  estimates$p_daily <- p_daily
  estimates$convergence <- fit$convergence 

  self$params_estimated <- estimates

  estimates
}

#' Run EM Algorithm
#'
#' Iterate through the Expectation-Maximization algorithm until convergence criteria are met.
#'
#' @param EM Taped RTMB model object. See details.
#' @param control List to control convergence criteria of the algorithm. See details.
#'
#' @details The EM object is built by other functions and takes parameter estimates, returns updated parameter estimates and the log-likelihood.
#' Typically we assume that object is taped and built with RTMB but that is not necessary, as long as it does the inner optimization of the objective function.
#' The control list input is
#' \itemize{
#'  \item `maxiters` Maximum iterations before failing (default = 1000).
#'  \item `tolerance` Maximum difference between log-likelihood between each iteration to set convergence (default = 1e-8).
#'  \item `relativeDifference` TRUE then convergence is taken as relative difference between the log-likelihood on each iteration, or FALSE it is total difference (default = FALSE).
#'  \item `verbose` TRUE then print more details upon fitting, or FALSE fit quietly unless there is an error.
#' }
#'
#' @return List of model estimates (`$params_estimated`), convergence details (`$convergence`), 
#' estimate of daily species passage (`$N_daily`), and daily proportions (`$p_daily`).
#'
#' @importFrom RTMB
#' @export
runEM <- function(EM, control){
  start <- EM$par()
  npar <- length(start)
  vals <- EM(start)
  if(any(is.nan(vals))) 
    stop("EM Algorithm finding NaN values. Either provide better initial values or more likely choose a different formulation. Occurs often with Test Fishery and Lengths don't agree.")
  fit_path <- c(start, -Inf)
  fit_path <- rbind(vals)

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
  restarts <- 0
  while(iter < maxit & !converged){
    vals <- tryCatch(EM(vals[1:npar]), warning = function(w){c(NaN, NaN, NaN)})
    if(!is.nan(vals[1])) {
      vals.prev <- vals
    }
    if(is.nan(vals[1])){
      vals <- vals.prev
      vals[npar+1] <- -Inf
      cat("[Warning]  Issues with initial values. Did not converge, returning the last value that successfully fitted.\n")
      break;
    }
    fit_path <- rbind(fit_path, vals)   
    if(relativeDiff & iter > 2){ diff <- abs((fit_path[iter, npar+1] - fit_path[iter-1, npar+1])/fit_path[iter-1, npar+1])
    }else{ diff <- abs(fit_path[iter, npar+1] - fit_path[iter-1, npar+1]) }
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
  return(list(values = vals, convergence = c("converged" = as.numeric(1-converged), "iters" = iter, "log_lik_difference" = as.numeric(diff)), fit_path = fit_path))
}

#' Basic Mixture Model
#'
#' Basic Gaussian mixture model fit using the EM algorithm.
#'
#' @param x data as a vector.
#' @param K number of components, default = NULL and ignored if names are given.
#' @param mu starting values for the mean (default = NULL).
#' @param sigma starting values for the standard deviation (default = NULL).
#' @param prob starting values for mixture proportions (default = NULL).
#' @param mu_fixed named vector with values to hold fixed (default = NULL).
#' @param sigma_fixed named vector with values to hold fixed (default = NULL).
#' @param wgts mixture weights associated with each observation x (default = NULL).
#' @param names Names of each component in the mixture model (default = NULL) Must provide names or K.
#' @param control List for specifying fitting details of EM algorithm (see details).
#'
#' @details The control list input
#' \itemize{
#'  \item `maxiters` Maximum iterations before failing (default = 1000).
#'  \item `tolerance` Maximum difference between log-likelihood between each iteration to set convergence (default = 1e-8).
#'  \item `relativeDifference` TRUE then convergence is taken as relative difference between the log-likelihood on each iteration, or FALSE it is total difference (default = FALSE).
#'  \item `verbose` TRUE then print more details upon fitting, or FALSE fit quietly unless there is an error.
#' }
#'
#' @return List of model estimates.
#'
#' @examples
#'\dontrun{
#' p <- c(0.2, 0.7, 0.1)
#' mu <- c(10, 25, 40)
#' sigma <- c(1.25, 2, 3.3)
#' x <- rnorm(200, mu[sample(3, 200, prob = p, replace = TRUE)], sigma[sample(3, 200, prob = p, replace = TRUE)])
#' fit <- basicMixtureModel(x, K = 3)
#'
#' est <- NULL
#' est2 <- NULL
#' for( i in 1:500 ){
#'    N <- c(2000, 1000)
#'    nsmp <- 25
#'    p1 <- c(0.3, 0.4, 0.3)
#'    p2 <- c(0.15, 0.65, 0.2)
#'    ptrue <- (p1*N[1] + p2*N[2])/sum(N)
#'    id1 <- sample(3, nsmp, prob = p1, replace = TRUE)
#'    id2 <- sample(3, nsmp, prob = p2, replace = TRUE)
#'    x1 <- rnorm(nsmp, c(25, 45, 65)[id1], c(3, 3.7, 2.5)[id1])
#'    x2 <- rnorm(nsmp, c(25, 45, 65)[id2], c(3, 3.7, 2.5)[id2])
#'    wgts <- c(rep(N[1]/nsmp, nsmp), rep(N[2]/nsmp, nsmp))
#'    fit <- basicMixtureModel(x=c(x1, x2), K = 3, prob = c(0.2, 0.4, 0.4), mu = c(20, 40, 60), wgts = wgts)
#'    fit2 <- basicMixtureModel(x=c(x1, x2), K = 3, prob = c(0.2, 0.4, 0.4), mu = c(20, 40, 60))
#'    est <- rbind(est, data.frame(par = c("p1", "p2", "p3"), diff = fit$proportion[order(fit$mu)] - ptrue, method = "weighted"))
#'    est <- rbind(est, data.frame(par = c("p1", "p2", "p3"), diff = fit2$proportion[order(fit2$mu)] - ptrue, method = "unweighted"))
#' }
#' boxplot(diff~method+par, data = est)
#' abline(h = 0, col='red')
#'}
#' @export
basicMixtureModel <- function(x, K = NULL, mu = NULL, sigma = NULL, prob = NULL, mu_fixed = NULL, sigma_fixed = NULL, wgts = NULL, names = NULL, control = list()){

  if(is.null(K) & is.null(names)) stop("Provide number of components, K")
  if(is.null(K)) K <- length(names)
  if(is.null(names)) names <- 1:K

  if(is.null(prob)) prob <- rep(1/K, K)
  if(is.null(mu)){
    mu <- rnorm(K, mean(x), sd(x))
    names(mu) <- names
    mu <- mu[!names %in% names(mu_fixed)]
    if(length(mu) == 0) mu <- NULL
  }
  if(is.null(sigma)){
    sigma <- rep(sd(x)/K, K)
    names(sigma) <- names
    sigma <- sigma[!names %in% names(sigma_fixed)]
    if(length(sigma) == 0) sigma <- NULL
  }
  if(is.null(wgts)) wgts <- rep(1, length(x))
  logit_prob <- fn_null(logitM, x=prob)
  
  pars_init <- list()
  if(!is.null(mu)) pars_init$mu <- mu
  if(!is.null(mu_fixed)) pars_init$mu <- mu[!names(mu) %in% names(mu_fixed)]  
  
  if(!is.null(sigma)) pars_init$log_sigma <- fn_null(log, x=sigma)
  if(!is.null(sigma_fixed)) pars_init$log_sigma <- pars_init$log_sigma[!names(sigma) %in% names(sigma_fixed)]  
  
  pars_init$logit_prob <- logit_prob
  pars_fixed <- list(mu_fixed = mu_fixed, log_sigma_fixed = fn_null(log, x=sigma_fixed))
  
  EMstep <- function(pars_outer){
    # Necessary in packages
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")

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
      # Necessary in packages
      "c" <- ADoverload("c")
      "[<-" <- ADoverload("[<-")
    
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
  EM$simplify()
  fit <- runEM(EM, control)
  pars_est <- reList(pars_init, fit$values)
  estimates <- list(mu = joinPars(pars_est$mu, pars_fixed$mu, names))
  estimates$sigma <- joinPars(fn_null(exp, x = pars_est$log_sigma), fn_null(exp, x = pars_fixed$log_sigma), names)
  estimates$proportion <- joinPars(fn_null(expitM, x=pars_est$logit_prob), fn_null(expitM, x=pars_fixed$logit_prob), names)
  estimates
}
