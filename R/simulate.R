
## extractControls
## Convenience copy from Nimble to extract values from control input of functions.
extractControls <- function(controlValue, defaultValue)
{
  if(!is.null(controlValue))  
    return(controlValue)
  else 
    return(defaultValue)
}

#' Simulate Lengths and Test Fishery Data:
#'
#' Function to simulate sonar lengths and test fishery data to test hydroacoustic species composition model.
#'
#' @param N Vector containing the total number of fish for each species.
#' @param nLengths Total number of lengths measurements to take.
#' @param sonarBins vector of sonar bins with which counts and lengths are taken. Assumed 1,2,3.
#' @param sonarBanks vector of sonar banks that counts and lengths are taken. Assumed 1 left and 2 right bank.
#' @param mu Vector of mean length for each species.
#' @param sigma Vector of standard deviation of length for each species.
#' @param beta Vector of correction terms for intercept and beam spreading effect.
#' @param sigmaObs Observational error associated with each measured length.0
#' @param q Vector of relative catchability of length equal to length(N) - 1, relative to the last fish.
#' @param control List that allows for controlling additional simulation parameters. See 'details' for a description of all the options.
#'
#'
#' @details Simulations are done such that fish proportions can vary through the day, but measurements are taken according to a stratified design. This is done by
#' simulating the number of fish according to a multinomial proportional to a normal distribution for each fish each hour, and each distance bin from shore, along with 
#' a sit proportion for each bank. Lengths are adjusted by any beam spreading terms provided, and distance from shore is drawn as a uniform distribution for each sonar bin
#' distances, assumed to be 10 m (e.g. sonar bin 1 is 0-10 m).
#'
#' User input options for control are
#' \itemize{
#'  \item{muTime}{Vector containing the mean times that each fish species migrates. Assumed to be 12.}
#'  \item{sigmaTime}{Vector containing the standard devation of times that each fish species migrates. Assumed to be 8. Fish timing is based on a normal distribution so large stanard deviation will tend towards no trend.}
#'  \item{muDist}{Vector containing the mean distance that each fish species migrates. Assumed to be 10.}
#'  \item{sigamDist}{Vector containing the standard deviation that each fish species migrates. Assumed to be 8. Fish distance is based on a normal distribution so large stanard deviation will tend towards no trend.}
#'  \item{bankProp}{Vector containing the proportion of each species on bank 1. Default is 0.5.}
#'  \item{tf_sonarBins}{Vector containing the distance bins that the fish are vulnerable to test fishery catch.}
#'  \item{tf_sonarBanks}{Vector containing the banks of the river where fish are vulnerable to test fishery catch.}
#'  \item{countHours}{The hours that counts are taken at. Defaults to 1:24}
#'  \item{lengthHours}{The hours that lengths are taken at. Defaults to seq(2, 24, 3)}
#'  \item{countDuration}{The number of minutes that fish are counted for. Defaults to 5.}
#'  \item{catchability}{The proportion of total fish, `sum(N)`, that are caught in the test fishery. Defaults to 0.001.}
#'  \item{ncaught}{The number of fish caught by the test fishery. Defaults to NULL but will override the catchability argument if provided.}
#' }
#'
#' @returns A list of all simulated values. 
#' \itemize{
#'  \item{count_data}{data frame containing all the hourly counts.}
#'  \item{length_data{data frame containing all the length measurements.}
#'  \item{test_data}{vector containing the test fishery counts.}
#'  \item{true_proportions}{data frame of the true proportions of each species each hour, bin, and bank, including the total number of fish, N.}
#' }
#'
#' @examples
#' datalist <- simulate(N = c(1000, 50000, 20000), nLengths = 20, 
#'                      sonarBins = c(1,2,3), sonarBanks = c(1,2), mu = c(38, 58, 70), sigma = c(4,3.5,7),
#'                      beta = c(1.5, 0.8), sigmaObs)
#' 
#' @export
simulate <- function(N = c(1000, 5000, 2000),
                     nLengths = 20,
                     sonarBins = c(1,2,3),
                     sonarBanks = c(1,2),
                     mu = c(38, 58, 70), 
                     sigma = c(4, 3.5, 7), 
                     beta = c(0,0), 
                     sigmaObs = 5.5,
                     q = c(1,1),
                     control = list()){

  nspp <- length(N)

  if(length(q) != nspp - 1)
    stop("Please Provide Relative Catchability with length = length of N-1")
  
  ## Adding run timing variability:
  muTime <- extractControls(control$muTime, 12)
  if(length(muTime) == 1){
    muTime = rep(muTime, nspp)
  }
  sigmaTime <- extractControls(control$sigmaTime, 8)
  if(length(sigmaTime) == 1){
    sigmaTime = rep(sigmaTime, nspp)
  }
  ## Adding run distance variability:
  muDist <- extractControls(control$muDist, 10)
  if(length(muDist) == 1){
    muDist = rep(muDist, nspp)
  }
  sigmaDist <- extractControls(control$sigmaDist, 8)
  if(length(sigmaDist) == 1){
    sigmaDist = rep(sigmaDist, nspp)
  }
  ## Proportion of each species on each river bank
  if(length(sonarBanks) == 1){
    bankProp <- rep(1, nspp)
  }else{
    bankProp <- extractControls(control$bankProp, 0.5)
    if(length(bankProp) == 1){
      bankProp <- rep(bankProp, nspp)
    }
  }
  ## Now figure out how to distribute fish across distance from shore.
  distProp <- matrix(0, nrow = length(sonarBins),  ncol = nspp)
  for( i in 1:nspp ){
    distProp[,i] <- dnorm((sonarBins - 1) * 10 + 5, muDist[i], sigmaDist[i])
    distProp[,i] <- distProp[,i]/sum(distProp[,i])
  }

  ## Now figure out what proportion of fish can be caught in the test fishery:
  tfBins <- extractControls(control$tf_sonarBins, sonarBins)
  tfBanks <- extractControls(control$tf_sonarBanks, sonarBanks)

  countHours <- extractControls(control$countHours, 1:24)
  ncount <- length(countHours)
  lengthHours <- extractControls(control$lengthHours, seq(2, 24, 3))
  nlengths <- length(lengthHours)
  countDuration <- extractControls(control$countDuration, 5)

  catchability <- extractControls(control$catchability, 0.001)
  ncaught <- extractControls(control$ncaught, NULL)

  countHours <- sort(countHours)
  prob_id <- expand.grid(Hour = countHours, SonarBin = sonarBins, SonarBank = sonarBanks) 
  
  ## Generate Count Data:
  Nt <- matrix(0,  nrow = nrow(prob_id), ncol = nspp)
  nt <- matrix(0,  nrow = nrow(prob_id), ncol = nspp)
  count_data <- data.frame()
  for( i in 1:nspp ){
    propt <- dnorm(countHours, muTime[i], sigmaTime[i])
    propt <- propt/sum(propt)
    prob_id$pt <- propt[prob_id[, "Hour"]]
    prob_id$pbank <- bankProp[i]*(prob_id$SonarBank == 1) + (1-bankProp[i])*(prob_id$SonarBank == 2) 
    prob_id$pdist <- distProp[prob_id$SonarBin, i]
    prob_id$prob <- (prob_id$pt*prob_id$pbank*prob_id$pdist)/sum(prob_id$pt*prob_id$pbank*prob_id$pdist)
    Nt[,i] <- rmultinom(1, size = N[i], prob = prob_id$prob)
    nt[,i] <- sapply(Nt[,i], FUN = function(x){rbinom(1, size = x, prob = countDuration/60)})
    dati <- data.frame(Species = i, SonarBank = prob_id$SonarBank, SonarBin = prob_id$SonarBin, Hour = prob_id$Hour, count = nt[,i], Duration = countDuration)
    count_data <- rbind(count_data, dati)
  }

  ## Generate Lengths:
  length_data <- data.frame()
  for( i in seq_along(prob_id[,1]) ){
     countsi <- count_data[count_data$SonarBank == prob_id$SonarBank[i] & count_data$SonarBin == prob_id$SonarBin[i],]
     D <- c((prob_id$SonarBin[i]-1)*10, prob_id$SonarBin[i]*10)
    for( j in lengthHours ){
      countsij <- countsi[countsi$Hour == j,]
      countsij <- countsij[order(countsij$Species), ]
      if(sum(countsij$count) == 0) next; ## Just skip lengths if zero fish counted.
      sppij <- sample(1:nspp, nLengths, replace = TRUE, prob = countsij$count/sum(countsij$count))
      length_data <- rbind(length_data, 
        data.frame(SonarBank = prob_id$SonarBank[i], SonarBin = prob_id$SonarBin[i], Hour = j, 
                   length.cm = rnorm(nLengths, mu[sppij], sigma[sppij]), R.m = runif(nLengths, D[1], D[2])) )
    }
  }
  ## Add in the observation error and other information:
  rollAngle <- extractControls(control$rollAngle, 0)
  rollCorrection <- 1/cos(rollAngle*pi/180)
  
  length_data$length.cm <- length_data$length.cm*rollCorrection + (beta[1] + beta[2]*length_data$R.m) + rnorm(nrow(length_data), 0, sigmaObs)
  
  ## Sum over all species:
  count_data_sum <- cbind(prob_id[, c("SonarBank", "SonarBin", "Hour")], data.frame(count = rowSums(nt)))
  count_data_sum$Duration <- countDuration

  ## Return the true proportions:
  true_proportions <- t(apply(Nt, 1, FUN = function(x){x/sum(x)}))
  colnames(true_proportions) <- paste("species", 1:nspp, sep = "_")
  true_proportions <- cbind(prob_id[, c("SonarBank", "SonarBin", "Hour")], true_proportions, data.frame(N = rowSums(Nt)) )
  
  ## Generate some test fishery data.
  q <- c(q, 1)
  r <- numeric(nspp)
  Ntf <- Nt[prob_id$SonarBin %in% tfBins & prob_id$SonarBank %in% tfBanks, ]
  Ntf <- colSums(Ntf)
  ptf <- Ntf/sum(Ntf)
  r <- q*ptf/sum(q*ptf)
  if(is.null(ncaught))
    ncatch <- ceiling(sum(N)*catchability)
  else
    ncatch <- ncaught
  tf_count <- rmultinom(1, size = ncatch, prob = r)[,1]
  names(tf_count) <- paste("species", 1:nspp, sep = "_")
  
  ## Return Data:
  datalist <- list(count_data = count_data_sum, length_data = length_data, test_data = tf_count, true_proportions)
  return(datalist)
}
