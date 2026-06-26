validate_estimates = function(self){
  ## Check if delta_mu on edges

  ## Check beta: Make sure bin1 > bin2 or that the line isn't too steep.
  
  ## Check qinv: Multivariate... Maybe need to do quadrature?

  ## pars_outer <- pars_init
  prior_nll <- function(pars_inner){

    # Necessary in packages
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")

    ## Parameter Processing
    ## 1) Mean Fish Size:
    mu <- pars_inner$mu
    ## Let the user fit a little shift in the fish size. This will default to 0.
    delta_mu <- ilogitInterval(pars_inner$logit_delta_mu, lower_delta_mu, upper_delta_mu)

    ## 2) Length standard deviation
    log_sigma <- pars_inner$log_sigma
    sigma_spp <- exp(log_sigma)
    ## length measurement error:
    log_sigma0 <- pars_inner$log_sigma0
    sigma0 <- exp(log_sigma0)

    ## 3) Catchability
    log_qinv <- pars_inner$log_qinv
    qinv <- exp(log_qinv)

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
    objval
  }

  species <- self$species_info$species
  nalpha <- length(obj$params_estimated$alpha)
  nalpha_jackchinook <- length(obj$params_estimated$alpha_jackchinook)
  nbeta <- length(obj$params_estimated$beta)
  lower_delta_mu <- self$data_list$lower_delta_mu
  upper_delta_mu <- self$data_list$upper_delta_mu

  pars <- list(mu = as.numeric(joinPars(self$params_init$mu, self$params_fixed$mu, species)))
  pars$log_sigma <- as.numeric(joinPars(fn_null(exp, x = self$params_init$log_sigma), fn_null(exp, x = self$params_fixed$log_sigma), species))
  pars$alpha <- as.numeric(joinPars(self$params_init$alpha, self$params_fixed$alpha, 1:nalpha))
  pars$alpha_jackchinook <- as.numeric(joinPars(self$params_init$alpha_jackchinook, self$params_fixed$alpha_jackchinook, 1:nalpha_jackchinook))
  pars$log_qinv <- as.numeric(joinPars(self$params_init$log_qinv, self$params_fixed$log_qinv, colnames(self$data_list$X_test_fishery)))
  pars$logit_delta_mu <- as.numeric(joinPars(self$params_init$logit_delta_mu, self$params_fixed$logit_delta_mu, species))
  pars$log_sigma0 <- as.numeric(joinPars(self$params_init$log_sigma0, self$params_fixed$log_sigma0, 1))
  pars$beta <- as.numeric(joinPars(self$params_init$beta, self$params_fixed$beta, 1:nbeta))

  ad_obj_prior <- MakeADFun(prior_nll, pars, silent = TRUE)
  fit_prior <- nlminb(ad_obj_prior$par, ad_obj_prior$fn, ad_obj_prior$gr)
  prior_sds <- sqrt(1/diag(ad_obj_prior$he(fit_prior$par)))

  ## Plot priors:
  prior_means <- ad_obj_prior$env$parList(fit_prior$par)
  prior_lower <- ad_obj_prior$env$parList(fit_prior$par - 5*prior_sds)
  prior_upper <- ad_obj_prior$env$parList(fit_prior$par + 5*prior_sds)
  prior_sds <- ad_obj_prior$env$parList(prior_sds)

  prior_dists <- list()

  dnames <- paste0("d", names(prior_means))
  for( i in 1:length(prior_means) ){
    y <- prior_means[[i]]
    for( j in seq_along(prior_means[[i]])){
      fn <- \(x){
        y[j] <- x
        if(grepl("log_", dnames[i])) y <- exp(y)
        if(grepl("logit_", dnames[i])) y <- ilogitInterval(y, lower_delta_mu, upper_delta_mu)
        if(dnames[i] == "dlogit_delta_mu"){ 
          ans <- self$prior_distributions[[dnames[i]]](y, lower = lower_delta_mu, upper = upper_delta_mu)
        }else{
          ans <- self$prior_distributions[[dnames[i]]](y)
        }
      }
      if(prior_sds[[i]][j] == Inf){
          xx <- seq(-50, 50, length = 300)
          df <- rep(1, 300)
      }else{
        xx <- seq(prior_lower[[i]][j], prior_upper[[i]][j], length = 300)
        df <- exp(do.call('c', lapply(xx, fn)))
      }
      if(grepl("log_", dnames[i])) xx <- exp(xx)
      if(grepl("logit_", dnames[i])) xx <- ilogitInterval(xx, lower_delta_mu[j], upper_delta_mu[j])
      prior_dists[[paste0(names(prior_means)[i], "_", j)]] <- data.frame(x = xx, f = df)
    }
  }
}

# i <- 20
# name <- names(prior_dists)[i]
# name <- gsub("log_|logit_", "", name)
# namei <- gsub("_.*", "", name)
# j <- as.numeric(gsub(".*_", "", name))
# namej <- names(obj$params_estimated[[namei]])[j]
# x <- obj$params_estimated[[namei]][j]
# rangex <- range(prior_dists[[i]]$x)
# if(rangex[1] > x) rangex[1] <- 1.01*x
# if(rangex[2] < x) rangex[2] <- 1.01*x

# plot(prior_dists[[i]], xlim = rangex, xlab = paste(namei, namej), type = 'l', ylab = "Prior Density")
# abline(v = x, col = 'red')
