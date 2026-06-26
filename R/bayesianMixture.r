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
fit_posterior <- function(self){
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
  
  ## pars_outer <- pars_init
  negll <- function(pars_inner){
  
    # Necessary in packages
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")

    ## Parameter Processing
    ## 1) Mean Fish Size:
    mu <- pars_inner$mu
    ADREPORT(mu)
    ## Let the user fit a little shift in the fish size. This will default to 0.
    delta_mu <- ilogitInterval(pars_inner$logit_delta_mu, lower_delta_mu, upper_delta_mu)
    REPORT(delta_mu)
    mu_adj <- mu + delta_mu
    ADREPORT(mu_adj)

    ## 2) Length standard deviation
    log_sigma <- pars_inner$log_sigma
    sigma_spp <- exp(log_sigma)
    ADREPORT(sigma_spp)
    ## length measurement error:
    log_sigma0 <- pars_inner$log_sigma0
    sigma0 <- exp(log_sigma0)
    ADREPORT(sigma0)
    ## Observed standard deviation
    sigma <- sqrt(exp(2*log_sigma) + exp(2*log_sigma0))

    ## 3) Catchability
    log_qinv <- pars_inner$log_qinv
    qinv <- exp(log_qinv)
    qinv_obs <- exp(X_test_fishery %*% log_qinv)[,1]
    ADREPORT(qinv)
    
    ## 4) Proportions relationships
    alpha <- pars_inner$alpha
    alpha_jackchinook <- pars_inner$alpha_jackchinook
    ADREPORT(alpha)
    ADREPORT(alpha_jackchinook)
    
    ## 5) Beam Spreading Corrections:
    beta <- pars_inner$beta
    ADREPORT(beta)
    
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
    ADREPORT(p_daily)
    N_daily <- estimate_daily_salmon(p_daily, self$data_list$total_salmon, species)  ## Combine adult chinook and remove smallresident fish.
    ADREPORT(N_daily)
    
    if("largeresident" %in% species) objval <- objval - self$prior_distributions$dlargeresident(N_daily[,1])
    
    ## Test fishery component:
    if(include_test_fishery){
      Ntest <- N_daily[cbind(test_fishery_catch$day, test_fishery_catch$N_index)]
      ptest <- 1/qinv_obs*test_fishery_catch$effort
      objval <- objval - test_fishery_weights*sum(dpois(x = test_fishery_catch$catch, Ntest*ptest, log = TRUE))
    }
    objval
  }

  ## Start with MAP:
  control <- self$fit_info
  control$verbose <- FALSE
  est_pars <- self$fitModel(control = control)

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
  fit.mle <- nlminb(ad_obj$par, ad_obj$fn, ad_obj$gr)
  sdrep <- sdreport(ad_obj)
  mle_sum <- data.frame(summary(sdrep, "fixed"))

  ## Start at MAP to just quickly approx the posterior:
  fit <- tmbstan(ad_obj, chains=1, iter = 2000, warmup = 500, init = "par")

  post_sum <- summary(fit)
  rnames <- rownames(post_sum$summary)
  post_sum <- data.frame(post_sum$summary)
  post_sum$parameters <- rnames
  post_sum <- post_sum |> subset(parameters != "lp__")
  
  which.max(mle_sum$Std..Error - post_sum$sd)
  mle_sum[16,]
  post_sum[16,]
  
  plot(mle_sum$Std..Error, post_sum$sd)
  abline(0,1, col = 'red')
  plot(mle_sum$Estimate, post_sum$mean)
  abline(0,1, col = 'red')
  
  
  plot(fit)
  post <- as.matrix(fit)
  lp <- post[, ncol(post)]
  post <- post[,-ncol(post)]

  i <- 1
  colnames(post) <- gsub("\\[.*", "", colnames(post))
  
  Ndaily <- NULL
  for(i in 1:nrow(post)){
    r <- ad_obj$report(post[i,])
  }
  hist(Ndaily[,2])
  plot(density(Ndaily[,2]))
  abline(v = est_pars$N_daily[3,3], col = 'red')
}
