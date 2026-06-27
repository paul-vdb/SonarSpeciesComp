#' Generate points along the prior distribution provided by the user.
#' 
#' @param self R6 object for sonar species composition analysis.
#' @param parameter Name of the parameter to generate the prior for.
#'
#' @export
prior_assess <- function(self, parameter){
  lower_delta_mu <- self$data_list$lower_delta_mu
  upper_delta_mu <- self$data_list$upper_delta_mu

  name <- gsub("log_|logit_", "", parameter)
  name_d <- paste0("d", parameter)
  if(parameter == "logit_delta_mu"){
    prior <- MakeTape(\(x){-self$prior_distributions[[name_d]](x, lower_delta_mu, upper_delta_mu)}, logitInterval(self$params_estimated[[name]], lower_delta_mu, upper_delta_mu))
  }else{ 
    if(grepl("log_", parameter)){
      prior <- MakeTape(\(x){-self$prior_distributions[[name_d]](exp(x))}, log(self$params_estimated[[name]]))
    }else{
      prior <- MakeTape(\(x){-self$prior_distributions[[name_d]](x)}, self$params_estimated[[name]])
    }
  }
  he_prior <- prior$jacfun()$jacfun()
  prior_newt <- prior$newton(1:length(self$params_estimated[[name]]))
  prior_mean <- prior_newt(numeric(0))
  prior_sd <- sqrt(1/diag(he_prior(prior_mean)))

  priors <- data.frame()  
  for( i in seq_along(prior_mean) ){
    mu <- prior_mean
    fni <- \(x){ 
      mu[i] <- x
      return(prior(mu))
    }
    xx <- seq(prior_mean[i] - 5*prior_sd[i], prior_mean[i] + 5*prior_sd[i], length = 300)
    f <- do.call('c', lapply(xx, fni))
    priors <- rbind(priors, data.frame(parameter = names(self$params_estimated[[name]])[i], x = xx, y = exp(-f)))
  }
  return(priors)
}