
#' Joint likelihood function for fitting with RTMB:
#'
#' @param self The R6 `speciesCompModel` that is set up for a model run. 
#'
#' @details Negative log likelihood function for fitting species composition model.
#'
#' @return Function that is the negative log likelihood.
#'
#' @export
make_negll <- function(self){

  species <- self$species_info$species
  
  nspp <- length(species)
  nchin <- self$species_info$nchinook
  nother <- self$species_info$nother
  
  adjust_lengths <- self$fit_info$adjust_lengths
  include_test_fishery <- self$fit_info$include_test_fishery
  proportion_adultchinook <- self$params_fixed[["proportion_adultchinook"]]  

  wgts <- self$data_list$length_data$weights
  ## Scale for natural comparison with test fishery.
  wgts <- wgts/max(wgts)
  obs_lengths <- self$data_list$length_data$L.cm.adj
  lower_delta_mu <- self$data_list$lower_delta_mu
  upper_delta_mu <- self$data_list$upper_delta_mu
  nbeta <- ncol(self$data_list$X_length)
  test_fishery_catch <- self$data_list$test_fishery_catch
  test_fishery_weights <- self$data_info$test_fishery_weights
  X_test_fishery <- self$data_list$X_test_fishery
  names_tf <- colnames(X_test_fishery)
  nobs <- nrow(self$data_list$length_data)

  negll <- function(pars_inner){

    # Necessary in packages
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")

    ## Parameter Processing
    ## 1) Mean Fish Size:
    mu <- pars_inner$mu
    ## Let the user fit a little shift in the fish size. This will default to 0.
    delta_mu <- ilogitInterval(pars_inner$logit_delta_mu, lower_delta_mu, upper_delta_mu)
    mu_adj <- mu + delta_mu
    
    ## 2) Length standard deviation
    log_sigma <- pars_inner$log_sigma
    sigma_spp <- exp(log_sigma)
    ## length measurement error:
    log_sigma0 <- pars_inner$log_sigma0
    sigma0 <- exp(log_sigma0)
    ## Observed standard deviation
    sigma <- sqrt(exp(2*log_sigma) + exp(2*log_sigma0))

    ## 3) Catchability
    log_qinv <- pars_inner$log_qinv
    qinv <- exp(log_qinv)
    qinv_obs <- exp(X_test_fishery %*% log_qinv)[,1]
    
    ## 4) Proportions relationships
    alpha <- pars_inner$alpha
    alpha_jackchinook <- pars_inner$alpha_jackchinook

    ## 5) Beam Spreading Corrections:
    beta <- pars_inner$beta
    
    ## Objective function for length based mixture model, conditional on posterior probs.
    objval <- 0
    ## Prior Distributions:
    objval <- objval - self$prior_distributions$dlog_qinv(qinv) - self$prior_distributions$dbeta(beta) - self$prior_distributions$dalpha_jackchinook(alpha_jackchinook) - 
        self$prior_distributions$dlog_sigma(sigma_spp) - self$prior_distributions$dmu(mu) - self$prior_distributions$dlogit_delta_mu(delta_mu, lower = lower_delta_mu, upper = upper_delta_mu) - 
        self$prior_distributions$dlog_sigma0(sigma0)
          
    ## Set up proportions. alpha parameter for predicting proportions. 
    p <- calcProportions(alpha = alpha, alpha_jackchinook = alpha_jackchinook, 
      p_adultchinook = proportion_adultchinook, Xprop = self$data_list$X_proportions, 
      K0 = nother, K = nspp)

    ## Adjust lengths for beam spreading:
    if(adjust_lengths){
      L <- obs_lengths - as.matrix(self$data_list$X_length) %*% beta
    }else{
      L <- obs_lengths
    }

    logpobs <- log(p[self$data_list$length_data$grpIndex,])

    for( i in 1:nobs ){
      log_pmixi <- (logpobs[i,] + dnorm(L[i], mu_adj, sigma, log = TRUE))
      maxlp <- max(log_pmixi)      
      objval <- objval - wgts[i]*(log(sum(exp(log_pmixi - maxlp))) + maxlp) ## Log sum exponential for numerical robustness.
    }
    ## Estimate N:
    p_daily <- estimate_daily_proportions(p, self$data_list$pred_df, species)
    N_daily <- estimate_daily_salmon(p_daily, self$data_list$total_salmon, species)  ## Combine adult chinook and remove smallresident fish.

    REPORT(p_daily)
    REPORT(N_daily)    
    if("largeresident" %in% species) objval <- objval - self$prior_distributions$dlargeresident(N_daily[,1])
    
    ## Test fishery component:
    if(include_test_fishery){
      Ntest <- N_daily[cbind(test_fishery_catch$day, test_fishery_catch$N_index)]
      ptest <- 1/qinv_obs*test_fishery_catch$effort
      objval <- objval - test_fishery_weights*sum(dpois(x = test_fishery_catch$catch, Ntest*ptest, log = TRUE))
    }
    objval
  }
  return(negll)
}

#' Compute posterior parameter values
#'
#'
#' @param self The R6 `speciesCompModel` that is set up for a model run.
#' @param post Matrix of simualted parameter values that are directly estimated.
#' @param ad_obj AD model object that was used to fit that we can then generate true values from.
#'
#' @return matrix of simulated values transformed to predicted values.
#'
#' @export
compute_param_values <- function(self, post, ad_obj){
  nbeta <- ncol(self$data_list$X_length)
  test_fishery_catch <- self$data_list$test_fishery_catch
  test_fishery_weights <- self$data_info$test_fishery_weights
  X_test_fishery <- self$data_list$X_test_fishery
  names_tf <- colnames(X_test_fishery)
  species <- self$species_info$species
  nspp <- length(species)
  nother <- self$species_info$nother
  proportion_adultchinook <- self$params_fixed[["proportion_adultchinook"]]  

  est_pars <- self$params_estimated

  Nd <- est_pars$N_daily |> subset(select=-Date) |> as.matrix()
  Pd <- est_pars$p_daily |> subset(select=-Date) |> as.matrix()
  Nnames <- paste0("N[",rep(1:nrow(Nd), ncol(Nd)), ",", rep(1:ncol(Nd), each = nrow(Nd)), "]")
  Pnames <- paste0("P[",rep(1:nrow(Pd), ncol(Pd)), ",", rep(1:ncol(Pd), each = nrow(Pd)), "]")
  betanames <- paste0("beta[", 1:length(est_pars$beta), "]")
  qinvnames <- paste0("qinv[", 1:length(est_pars$qinv), "]")
  sigma0names <- "sigma0[1]"
  deltanames <- paste0("delta_mu[", 1:length(est_pars$delta_mu), "]")
  muadjnames <- paste0("mu_adjusted[", 1:length(est_pars$mu_adjusted), "]")
  alphajacknames <- paste0("alpha_jack[", 1:length(est_pars$alpha_jackchinook), "]")
  alphanames <- paste0("alpha[", 1:length(est_pars$alpha), "]")
  munames <- paste0("mu[", 1:length(est_pars$mu), "]")
  sigmanames <- paste0("sigma[", 1:length(est_pars$sigma), "]")

  names <- c(munames, sigmanames, alphanames, alphajacknames, deltanames, muadjnames, sigma0names, betanames, qinvnames, Pnames, Nnames)
  output <- matrix(NA, ncol = length(names), nrow = nrow(post))
  colnames(output) <- names
  
  for( i in 1:nrow(post) ){
    parsi <- ad_obj$env$parList(post[i,])

    output[i, betanames] <- parsi$beta
    output[i, munames] <- parsi$mu
    output[i, sigmanames] <- exp(parsi$log_sigma)
    output[i, deltanames] <- ilogitInterval(parsi$logit_delta_mu, self$data_list$lower_delta_mu, self$data_list$upper_delta_mu)
    output[i, muadjnames] <- parsi$mu + output[i, deltanames]
    output[i, sigma0names] <- exp(parsi$log_sigma0)
    output[i, qinvnames] <- exp(parsi$log_qinv)
    output[i, alphanames] <- parsi$alpha
    output[i, alphajacknames] <- parsi$alpha_jackchinook

    ## Set up proportions. alpha parameter for predicting proportions. 
    p <- calcProportions(alpha = parsi$alpha, alpha_jackchinook = parsi$alpha_jackchinook, 
      p_adultchinook = proportion_adultchinook, Xprop = self$data_list$X_proportions, 
      K0 = nother, K = nspp)
    p_daily <- estimate_daily_proportions(p, self$data_list$pred_df, species)
    N_daily <- estimate_daily_salmon(p_daily, self$data_list$total_salmon, species)  ## Combine adult chinook and remove smallresident fish.
    output[i, Nnames] <- c(N_daily)
    output[i, Pnames] <- c(p_daily)
  }
  output
}

#' Fit joint species composition model as a Bayesian model.
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
#' @import tmbstan
#' @export
run_full_model <- function(self, bayesian = TRUE){
  species <- self$species_info$species
  species_N <- self$species_info$species_predict
  
  nspp <- length(species)
  nchin <- self$species_info$nchinook
  nother <- self$species_info$nother
  p_adultchinook <- 1
  ndays <- self$ndays
  
  adjust_lengths <- self$fit_info$adjust_lengths
  include_test_fishery <- self$fit_info$include_test_fishery

  wgts <- self$data_list$length_data$weights
  ## Scale for natural comparison with test fishery.
  wgts <- wgts/max(wgts)  
  obs_lengths <- self$data_list$length_data$L.cm.adj
  lower_delta_mu <- self$data_list$lower_delta_mu
  upper_delta_mu <- self$data_list$upper_delta_mu
  nbeta <- ncol(self$data_list$X_length)
  test_fishery_catch <- self$data_list$test_fishery_catch
  test_fishery_weights <- self$data_info$test_fishery_weights
  X_test_fishery <- self$data_list$X_test_fishery
  names_tf <- colnames(X_test_fishery)

  pars_init <- self$params_init
  pars_fixed <- self$params_fixed
  proportion_adultchinook <- pars_fixed[["proportion_adultchinook"]]
  pars_fixed[["proportion_adultchinook"]] <- NULL
  
  ## Need to ensure that beta is fixed if adjust_lengths = FALSE.
  if(!adjust_lengths){
    pars_fixed$beta <- c(0, 0)
    pars_init$beta <- NULL
  }

  nalpha <- length(c(pars_init$alpha, pars_fixed$alpha))
  nalpha_jackchinook <- length(c(pars_init$alpha_jackchinook, pars_fixed$alpha_jackchinook))
  nobs <- nrow(self$data_list$length_data)
  
  negll <- make_negll(self)
  
  ## Start with MAP:
  control <- self$fit_info
  control$verbose <- FALSE
  if(self$est_date == max(self$params_estimated$N_daily$Date)){
    est_pars <- self$params_estimated
  }else{
    est_pars <- self$fitModel(control = control)
  }

  ## 1) Mean Fish Size:
  inits <- list()
  inits$mu <- as.numeric(est_pars$mu)
  inits$logit_delta_mu <- logitInterval(as.numeric(est_pars$delta_mu), lower_delta_mu, upper_delta_mu)
  ## 2) Length standard deviation
  inits$log_sigma <- as.numeric(log(est_pars$sigma))  ## *** data input
  ## length measurement error:
  inits$log_sigma0 <- as.numeric(log(est_pars$sigma0))
  ## 3) Catchability
  inits$log_qinv <- as.numeric(log(est_pars$qinv))
  ## 4) Proportions relationships
  inits$alpha <- as.numeric(est_pars$alpha)
  inits$alpha_jackchinook <- as.numeric(est_pars$alpha_jackchinook)
  ## 5) Beam Spreading Corrections:
  inits$beta <- as.numeric(est_pars$beta)

  par_map <- list()
  par_names <- gsub("log_|ilogit_", "", names(pars_fixed))
  for( i in 1:length(pars_fixed) ){
    var_id <- names(pars_fixed)[i]
    pari <- 1:length(inits[[var_id]])
    pari[names(est_pars[[par_names[i]]]) %in% names(pars_fixed[[i]])] <- NA
    par_map[[var_id]] <- as.factor(pari)
  }
  ad_obj <- MakeADFun(negll, inits, map = par_map, silent = TRUE)  
  nsim <- extractControls( control$nsim, 1000 )
  warmup <- extractControls( control$warmup, 500 )
  if(!bayesian){
    fit.mle <- nlminb(ad_obj$par, ad_obj$fn, ad_obj$gr)
    sdrep <- sdreport(ad_obj)
    mle_sum <- data.frame(summary(sdrep, "fixed"))    
    npars <- nrow(mle_sum)
    post <- NULL
    for( i in 1:nsim) post <- rbind(post, mle_sum[,"Estimate"] + rnorm(npars, 0, mle_sum[, "Std..Error"]))
  }else{
    fit <- tmbstan(ad_obj, chains=1, iter = nsim + warmup, warmup = warmup, init = "par")
    post <- as.matrix(fit)
    post <- post[, colnames(post) != "lp__"]
  }
  output <- compute_param_values(self, post, ad_obj)
  self$.posterior_sample <- output

  spp_alpha <- names(self$data_list$X_proportions)[do.call("c", lapply(1:length(self$data_list$X_proportions), FUN = function(x){rep(x, ncol(self$data_list$X_proportions[[x]]))}))]
  alpha_names <- paste0(spp_alpha, "_", do.call("c", lapply(self$data_list$X_proportions, colnames)))
  pnames <- paste0(rep(species, each = nrow(est_pars$p_daily)), "_", rep(est_pars$p_daily$Date, each = length(species)))
  dnames <- paste0(rep(species_N, each = nrow(est_pars$N_daily)), "_", rep(est_pars$N_daily$Date, length(species_N)))

  sum.out <- data.frame(parameter = gsub("\\[.*", "", colnames(output)),
    type = c(species, species, alpha_names, species, species, 1, names(est_pars$beta), names(est_pars$qinv), pnames, dnames),
    param = colnames(output)
  )
  sum.out <- cbind(sum.out, "mean" = apply(output, 2, mean), "sd" = apply(output, 2, sd), t(apply(output, 2, quantile, c(0.025, 0.5, 0.975))))
  self$.posterior_summary <- sum.out
}
