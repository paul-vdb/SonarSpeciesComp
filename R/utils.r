
#' Join fixed and estiamted parameter values into a single vector.
#'
#' @param est values that might be estimated.
#' @param fixed values that might be fixed
#' @param species names of the values to join in correct order.
#'
#' @return vector of combined fixed and estimated parameter values.
#'
#' @export
joinPars <- function(est = NULL, fixed = NULL, species){
  if(is.null(est) & !is.null(fixed)) return(fixed)
  nspp <- length(species)
  v <- numeric(nspp)
  idx <- which(species %in% names(fixed))
  if(length(idx) > 0) v[idx] <- fixed
  idy <- which(!species %in% names(fixed))
  if(length(idy) > 0) v[idy[1:length(est)]] <- est
  if(!RTMB:::ad_context()) names(v) <- species[1:length(v)]
  return(v)
}

#' log function with null values.
#'
#' @param x value to log.
#'
#' @return log(x) or NULL if is.null(x).
#'
#' @export
log_null <- function(x){
  if(is.null(x)) return(x)
  return(log(x))
}

#' Call functions with null values.
#'
#' @param f function to do transformation.
#' @param ... additional function information.
#'
#' @details Generalizes R functions to deal with NULL values.
#'
#' @return Function that returns NULL if input is NULL.
#'
#' @export
fn_null <- function(f, ...){
  val <- list()
  val[names(list(...))] <- list(...)
  if(is.null(val$x)) return(val$x)
  fn <- as.call(list(f, quote(...)))
  return(eval(fn))
}

#' Inverse Multivariate Logit Transform
#'
#' @param x vector to transform
#'
#' @details Transform from (-Inf, Inf) -> (0, 1).
#'
#' @return transformed value adding an extra dimension such that sums to 1.
#'
#' @export
expitM <- function(x){
  ex <- exp(c(x,0))
  ex/sum(ex)
}

#' Multivariate Logit Transform
#'
#' @param x vector to transform
#'
#' @details Transform from (0, 1) -> (-Inf, Inf).
#'
#' @return transformed value of reduced dimension.
#'
#' @export
logitM <- function(x){
  K <- length(x)
  log(x[1:(K-1)]/x[K])
}

#' Interval transform logit
#'
#' @param x Value to transform
#' @param lower Lower bound of the value to be transformed
#' @param upper Upper bound of the value to be transformed
#'
#' @details Transform from (lower, upper) -> (-Inf, Inf).
#'
#' @return transformed value
#'
#' @export
logitInterval <- function(x, lower, upper){
  p <- (x-lower)/(upper - lower)
  log(p/(1-p))
}

#' Log Determinant Jacobian of logit interval transformation.
#'
#' @param x Value to transform
#' @param lower Lower bound of the value to be transformed
#' @param upper Upper bound of the value to be transformed
#'
#' @details Compute the jacobian of the transformation to go from (lower, upper) -> (-Inf, Inf).
#'
#' @return log determinant of jacobian.
#'
#' @export
logDetJac_logitInterval <- function(x, lower, upper){
  delta <- upper-lower
  gr <- delta/((lower-x)*(lower+delta-x))
  sum(-log(abs(gr)))
}

#' Inverse interval transform logit
#'
#' @param x Value to transform
#' @param lower Lower bound of the value to be transformed
#' @param upper Upper bound of the value to be transformed
#'
#' @details inverse logit transformation for a value centered at zero to become (0,1). Then scaled to lower and upper.
#'
#' @return value on original scale (lower - upper).
#'
#' @export
ilogitInterval <- function(x, lower, upper){
  p <- 1/(1+exp(-x))
  lower + p*(upper-lower)
}

#' Return jacobian function
#'
#' @param fn Jacobian for the specified transformation.
#' @param includeJacobian Whether to actually use a jacobian transformation.
#'
#' @return A function, either the jacobian originally passed as fn, or a function that only returns 0.
#'
#' @export
addJacobian = function(fn, includeJacobian){
  if(!includeJacobian){
    return(\(...){0})
  }else{
    return(fn)
  }
}

#' Put parameters into original list format.
#'
#' @param parsList original pars list to push values to.
#' @param parsVec vector of par values to push into list.
#'
#' @return updated parsList with values from parsVec.
#'
#' @export
reList <- function(parsList, parsVec){
  count <- 1
  for( i in 1:length(parsList) ){
    nL <- length(parsList[[i]])
    parsList[[i]] <- parsVec[count:(count + nL - 1)]
    count <- count + nL
  }
  parsList
}

#' Extract controls
#'
#' @param controlValue
#' @param defualtValue
#'
#' @return either control value if it is not NULL, or default value.
#'
#' @export
extractControls <- function(controlValue, defaultValue){
  if(!is.null(controlValue))
    return(controlValue)
  else 
    return(defaultValue)
}

#' Extract parameter values
#'
#' @param self R6 speciesCompModel object.
#' @param init_values initial values as a named vector.
#' @param fixed_values values that are to be held fixed in the model.
#' @param default_values values to fill in if init_values are not provided and are not fixed.
#' @param name vector of the names of the vector values that need to be returned.
#' @param transform function to transform the values to the scale needed for estimation.
#' @param ... additional values to pass to the transformation.
#'
#' @return No return, R6 object is updated for both `$params_init` and `$params_fixed`
#'
#' @export
extractParams <- function(self, init_values=NULL, fixed_values=NULL, default_values, name, transform = I, ...){
  id <- names(default_values)
  if(is.null(id)){
    id <- 1:length(default_values)
    names(default_values) <- id
  }
  self$params_fixed[[name]] <- NULL
  self$params_init[[name]] <- NULL
  
  ## Set fixed values, keep NULL if note present:
  self$params_fixed[[name]] <- fn_null(transform, x=fixed_values, ...)

  init <- default_values
  if(!is.null(init_values)){
    if(is.null(names(init))) names(init) <- 1:length(init)
    init[names(init)] <- init
  }
  if(!is.null(fixed_values)){ 
    if(is.null(names(fixed_values))) names(fixed_values) <- 1:length(fixed_values)
    init <- init[!names(init) %in% names(fixed_values)]
    if(length(init) == 0) init <- NULL
  }
  self$params_init[[name]] <- fn_null(transform, x = init, ...)
}

#' Daily proportion of salmon
#'
#' Combine proportions based on estimates for covariates to the daily scale.
#'
#' @param prop Proportion of fish at level of provided covariates, estimated by species composition model.
#' @param pred_df dataframe with information needed to aggregate prop to daily estimates.
#' @param species vector of standard deviations for each componenet (e.g. species).
#'
#' @return Matrix of daily proportions for each species.
#' 
#' @export
estimate_daily_proportions <- function(prop, pred_df, species){
  ndays <- length(unique(pred_df$day))
  nspp <- length(species)
  
  p_daily <- matrix(0, nrow = ndays, ncol = nspp) 
  for( d in 1:ndays ){
    dindx <- which(pred_df$day == d)
    for( i in 1:nspp ) {
      p_daily[d, i] <- sum(pred_df$Count[dindx]*prop[dindx,i])
    }
    p_daily[d, ] <- p_daily[d, ] / sum( p_daily[d, ] )
  }
  return(p_daily)
}

#' Log Sum Exponential
#' 
#' @param x vector on the log scale to sum.
#'
#' @details Find the max value, subtract that on the log scale and exponentiate and sum. Then return back on log scale.
#'
#' @return stable version of \code{sum(x)}.
#' 
#' @export
log_sum_exp(x){
  C <- max(x)
  log(sum(exp(x - C))) + C ## log sum exponential
}

#' Estimate daily salmon
#'
#' Combine proportions and total salmon to estimate daily species escapement
#'
#' @param p_daily Proportion of fish each day, estimated by species composition model.
#' @param total_salmon Daily salmon escapement counted by the hydroacoustic program.
#' @param species vector of standard deviations for each componenet (e.g. species).
#'
#' @details If small resident are included, they are removed when estimating total daily salmon.
#' Large resident fish are included (> 35 cm for updated Mission protocols).
#'
#' @return Matrix for each species, excluding small resident fish, with an additional column added that
#' combines all Chinook salmon.
#' 
#' @export
estimate_daily_salmon = function(p_daily, total_salmon, species){
  if("smalladultchinook" %in% species) {
    chin_indx <- grep("adultchinook", species)
    p_daily[, chin_indx[1]] <- p_daily[,chin_indx[1]] + p_daily[,chin_indx[2]]
    p_daily <- p_daily[,-chin_indx[2], drop = FALSE]
  }
  
  if("smallresident" %in% species){
    correction <- 1 - p_daily[,1]
    p_daily <- p_daily[, -1, drop = FALSE]
  }else{
    correction <- rep(1, nrow(p_daily))
  }
  N_salmon <- matrix(0, nrow = nrow(p_daily), ncol = ncol(p_daily))
  for( d in 1:nrow(p_daily) ){
    N_salmon[d, ] <- as.numeric(p_daily[d,])*total_salmon[d, "count"]/correction[d]
  }
  N_salmon <- cbind(N_salmon, N_salmon[, ncol(N_salmon)-1] + N_salmon[, ncol(N_salmon)])
  return(N_salmon)
}

#' Find modes of a vector.
#'
#' Find the local modes using a kernel density approximation.
#'
#' @param x vector of values to find the modes of.
#' @param ... additional values to pass to the `density` function.
#'
#' @details Any time the density increases and then decreases, a mode is identified.
#'
#' @return Return modes found, or NA if none are found.
#' 
#' @export
findModes <- function(x,...) {
  modes <- NULL
  f <- density(x, ...)  
  for ( i in 2:(length(f$y)-1) ){
    if ( (f$y[i] > f$y[i-1]) & (f$y[i] > f$y[i+1]) ) {
      modes <- c(modes,i)
    }
  }
  if ( length(modes) == 0 ) {
    return(NA)
  }
  return(f$x[modes])
}

#' Find global mode of a vector.
#'
#' Find the global mode using a kernel density approximation.
#'
#' @param x vector of values to find the modes of.
#' @param xrange range of values to look for the global mode in.
#' @param ... additional values to pass to the `density` function.
#'
#' @details Returns largest value found by density approximation.
#'
#' @return Returns scalar.
#' 
#' @export
findGlobalMode <- function(x, xrange = NULL, ...) {
  if(!is.null(xrange)){
    x <- x[x <= xrange[2] & x >= xrange[1]]
  }
  if(length(x) < 3) return(NA)
  f <- density(x, ...) 
  xmax <- f$x[which.max(f$y)]  
  return(xmax)
}

#' Is it a date
#' 
#' Check if the value inherits Date.
#'
#' @param x single value to check if date type.
#'
#' @export
is.Date <- function(x) inherits(x, 'Date')

#' Process Minutes
#'
#' @param time Character in H:M:S format.
#'
#' @return Time in hours.
#'
#' @export
process_minutes <- function(time){
  time <- strsplit(time, ":")
  do.call("c", lapply(time, FUN = function(x){as.numeric(x[1])*60 + as.numeric(x[2]) + as.numeric(x[3])/60}))
}
