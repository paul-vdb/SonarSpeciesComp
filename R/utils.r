## Utility Functions

## Inverse Multivariate Logit Transform
expitM <- function(x){
  ex <- exp(c(x,0))
  ex/sum(ex)
}

## Multivariate Logit Transform
logitM <- function(p){
  K <- length(p)
  log(p[1:(K-1)]/p[K])
}
## Stick Breaking Transformation
stickBreakTransform <- function(p){
  K <- length(p) - 1
  runningSum <- 0
  x <- numeric(K)
  x[1] <- log(p[1]/(1-p[1]))
  for(i in 2:K){
    runningSum <- runningSum + p[i-1]
    pi <- p[i]/(1-runningSum)
    x[i] <- log(pi/(1-pi))
  }
  x
}

## Inverse Stick Breaking Transformation
inverseStickBreakTransform <- function(x){
  K <- length(x) + 1
  p <- numeric(K)
  p[1] <- 1/(1+exp(-x[1]))
  if(K > 2){
    runningSum <- 0
    for( i in 2:(K-1) ){
      runningSum <- runningSum + p[i-1]
      p[i] <- (1-runningSum)*1/(1+exp(-x[i]))
    }
  }
  p[K] <- 1-sum(p[1:(K-1)])
  p
}

## Transformations for Prior Distributions
## This one if prior is given in terms of p.
## e.g. p ~ beta, logitp ~ beta*|p*(1-p)|
logDetJac_logitM <- function(x){
  sum(log(x))
}

## This one if prior is given in terms of y
## e.g. y ~ gamma, logy ~ gamma*|y|
logDetJac_log <- function(x){
  log(x)
}

## Convert parameter list
reList <- function(parsList, parsVec){
  count <- 1
  for( i in 1:length(parsList) ){
    nL <- length(parsList[[i]])
    parsList[[i]] <- parsVec[count:(count + nL - 1)]
    count <- count + nL
  }
  parsList
}

## Convenience function for extracting control values.
extractControls <- function(controlValue, defaultValue)
{
  if(!is.null(controlValue))  
    return(controlValue)
  else 
    return(defaultValue)
}

## Estimate total species per day / bin for test fishery counts:
estimateDailyN <- function(testFishGrp, prob, q){
  if( nrow(testFishGrp) != nrow(prob) )
  Nspp <- 
  for( i in seq_along(prob) ){
    testFishGrp$totalFish*prob[,i]
  }
      
      idx <- which(tfgrp$day == d)
      nd <- 
      Nd <- sum(nd)
      if(Nd == 0) next
      pd <- numeric(nq)
      for( k in 1:nq ) pd[k] <- sum(p[idx,k+(K-nq)]*nd/Nd)              
      prob <- (q*pd)/sum(q*pd)
      negll <- negll - testwgts[d]*dmultinom(testcounts[d,], prob = prob, size = sum(testcounts[d,]), log = TRUE)
}

## Find and return visible modes based on a density function:
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

## Find and return the global mode based on a density function:
findGlobalMode <- function(x, xrange = NULL, ...) {
  if(!is.null(xrange)){
    x <- x[x <= xrange[2] & x >= xrange[1]]
  }
  if(length(x) < 3) return(NA)
  f <- density(x, ...) 
  xmax <- f$x[which.max(f$y)]  
  return(xmax)
}

## Not in Base R? Really...
is.Date <- function(x) inherits(x, 'Date')


## Taken from plotrix for the weighted histogram, reworked for this pacakge.
histplot <- function(x, wgts, range = c(35, 120), delta = 2, ylim = NA, xlim = NA, xlab = NA, ...){
    if (missing(x)) 
        stop("Must pass a value to x.")
    if (missing(wgts)) 
        wgts <- rep(1, length(x))
    extra <- delta - diff(range) %% delta
    if(extra != 0) range <- range + c(-extra/2, extra/2)
    n <- diff(range)/delta
    breaks <- seq(range[1], range[2], by = delta)
    
    counts <- NULL
    for ( bin in 1:(length(breaks)-1) ){
      counts <- c(counts, sum(wgts[x >= breaks[bin] & x < breaks[bin + 1]]))
    }
    density <- counts/sum(counts)
    width <- rep(delta, length(density))
    density <- density/width
    
    if (is.na(ylim)) ylim <- c(0, 1.1 * max(density, na.rm = TRUE))
    if (is.na(xlab)) xlab = "Fish Length (cm)"
    mids <- barplot(density, width =  width, col = "grey", space = 0, 
        ylim = ylim, ylab = "Density", xlab = xlab, ...)
    tickpos <- c(mids - width/2, mids[length(mids)] + width[length(width)]/2)
    axis(1, at = tickpos, labels = signif(breaks, 3))
    return(list(range = range, breaks = breaks, density = density))
}

