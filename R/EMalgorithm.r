## EM Algorithm Solver:
library(RTMB)
library(mixtools)
library(R6)
source("utils.r")
source("distributions.r")
source("processInputs.r")

## Set up an R6 Environment to run the mixture model and keep all the pieces for later.

## Delete me later:
##------------------------------------------------------------------------------------
library(tidyverse)
library(readxl)
setwd("C:/Users/vandambatesp/Documents/GitHub/SonarSpeciesComp/R")
sonarLengths <- read.csv("../example_data/2023MergedLengths.csv")
sonarCounts <- read.csv("../example_data/SONARMERGE export_2023.csv")
testFisheryCounts <- read_excel("../example_data/2023_WhonnockSpeciesComposition.xlsx")
driftTimes <- read_excel("../example_data/Whonnock Drift Times.xlsx")
driftTimes <- driftTimes |> subset(FisheryName == "Area 29 - Whonnock Sockeye Gillnet" & as.numeric(format(TripDate, "%Y")) > 2015)
driftTimes <- driftTimes %>% group_by(TripDate) %>% summarize(setStart = min(FE_NET_START_OUT_DTT), nSets = n(), .groups = "drop")
testFisheryCounts <- testFisheryCounts |> left_join(driftTimes, by = c("TRIP_DTT" = "TripDate"))
estDate <- as.Date("2023-08-15")

othersalmon <- read_excel("../example_data/Whonnock_FLengths_2023-2024_OtherSalmon.xlsx", sheet = "Whonnock 2023", skip = 1)
othersalmon <- othersalmon %>% select(TestFishery = `Test Fishery`, Species, Fork = `Fork (mm)`, POH = `POH (mm)`, WeightLbs = `Weight (lbs)`, Sex, Date)
pssalmon <- read_excel("../example_data/Scale Database Export 2013-2023 Pink and Sockeye AB TF .xlsx")
pssalmon <- pssalmon %>% select(TestFishery = Locality, Species = SpeciesCode, POF, POH, Fork, WeightLbs, Sex, Date = StartDate, Set, CatchPosition)
testFisheryLengths <- pssalmon %>% mutate(Species = ifelse(Species == "S", "Sockeye", "Pink")) %>% 
  bind_rows(othersalmon) %>%
  mutate(FL.cm = Fork/10, POF.cm = POF/10, POH.cm = POH/10) %>% 
  filter(!is.na(Date) & format(Date, "%Y") == "2023")

##------------------------------------------------------------------------------------

speciesComp <- speciesCompSummary$new(species = c("largeresident", "jackchinook", "sockeye", "smalladultchinook", "largeadultchinook"), site = "Mission", estDate = "2023-08-05")
speciesComp$processData(sonarCounts, sonarLengths, testFisheryCounts) ## Process all data:
speciesComp$setDate("2023-08-10", 3)
speciesComp$setSpeciesLengths(testFisheryLengths = testFisheryLengths, ndays = 10)
speciesComp$setModelProportions()
speciesComp$setLengthAdjustment()
speciesComp$setModelParameters()


parsOuter <- speciesComp$parsInit
parsFixed <- speciesComp$parsFixed
dataList <- speciesComp$analysisData
adjustLengths <- speciesComp$adjustLengths
includeTestFishery <- speciesComp$includeTestFishery

EMstep <- function(pars_outer){
  getAll(parsOuter, parsFixed, warn = FALSE)
  ll <- 0
  
  ## Parameter Processing
  Kchin <- length(muChin)
  K0 <- length(mu)
  K <- Kchin + K0
  mu <- c(mu, muChin)
  q <- exp(c(logq, 0))

  ## Standard deviation incld. observation error
  logsigma_ <- c(logsigma, logsigmaChin)
  sigma <- sqrt(exp(2*logsigma_) + exp(2*logsigma0))

  ## Set up proportions. alpha parameter for predicting proportions. 
  ## Xalpha is a list of design matrices for each species.
  np <- nrow(dataList$lengthData)
  logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
  indx0 <- 1
  for( i in 1:K0 ){
    nc <- ncol(dataList$Xprop[[i]]) 
    indx1 <- indx0 + nc - 1
    logitp[,i] <- dataList$Xprop[[i]] %*% alpha[indx0:indx1]
    indx0 <- indx1 + 1
  }

  logitpjack <- dataList$XpropChin %*% alphaJackChinook
  pjack <- 1/(1+exp(-logitpjack))

  p <- matrix(0, nrow = np, ncol = K)
  for( i in 1:np ) {
    p[i, 1:(K0+1)] <- expitM(logitp[i,])  
    p[i, (K0+1):K ] <- p[i, (K0+1)]*c(pjack[i], (1-pjack[i])*pAdultChinook)
  }

  ## Adjust lengths for beam spreading:
  ## X is a design matrix that is either distance from shore, or bin id. 
  ## It is expected to contain an intercept term.
  if(adjustLengths){
    L <- dataList$lengthData$L.cm - dataList$Xlength %*% beta
  }else{
    L <- dataList$lengthData$L.cm
  }
  wgts <- dataList$lengthData$weights
  
  ## Prior or penalty terms:
  # ll <- ll + prior_beta(beta) + prior_sigma0(exp(logsigma0)) + prior_sigma(exp(logsigma)) + prior_mu(mu) ## Add prior to q.
  
  ## Calculate Posterior Probabilities + log likelihood:
  outerLL <- calcPostProb(x = L, mu = mu, sigma = sigma, prob = p, wgts = wgts)
  postProb <- outerLL$postp
  ll <- ll + outerLL$ll

  ## Test fishery component:
  if(includeTestFishery){
    ndays <- nrow(dataList$testFisheryCounts)
    padultchinook <- 1-rowSums(p[,1:(K0+1)])
    nq <- ncol(dataList$testFisheryCounts)
    for( d in 1:ndays ){
      idx <- which(dataList$lengthData$day == d)
      Nd <- sapply(dataList$qSppIndices, FUN = function(x){sum(p[idx,x]*dataList$lengthData$SalmonCount[idx])})
      Ndc <- sum(padultchinook[idx]*dataList$lengthData$SalmonCount[idx])
      Nd <- c(Nd, Ndc)
      if(sum(Nd) == 0) next
      prob <- (q*Nd)/sum(q*Nd)
      ll <- ll + dataList$testFisheryWeights[d]*dmultinom(dataList$testFisheryCounts[d,], prob = prob, size = sum(dataList$testFisheryCounts[d,]), log = TRUE)
    }
  }

  ## Inner Objective Function
  ## ------------------------------------------
  inner_objective_fn <- function(parsInner){
    getAll(parsInner, parsFixed, warn = FALSE)

    ## Parameter Processing
    Kchin <- length(muChin)
    K0 <- length(mu)
    K <- Kchin + K0
    mu <- c(mu, muChin)
    q <- exp(c(logq, 0))

    ## Standard deviation incld. observation error
    logsigma_ <- c(logsigma, logsigmaChin)
    sigma <- sqrt(exp(2*logsigma_) + exp(2*logsigma0))

    ## Set up proportions. alpha parameter for predicting proportions. 
    ## Xalpha is a list of design matrices for each species.
    np <- nrow(dataList$lengthData)
    logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
    indx0 <- 1
    for( i in 1:K0 ){
      nc <- ncol(dataList$Xprop[[i]]) 
      indx1 <- indx0 + nc - 1
      logitp[,i] <- dataList$Xprop[[i]] %*% alpha[indx0:indx1]
      indx0 <- indx1 + 1
    }

    logitpjack <- dataList$XpropChin %*% alphaJackChinook
    pjack <- 1/(1+exp(-logitpjack))

    p <- matrix(0, nrow = np, ncol = K)
    for( i in 1:np ) {
      p[i, 1:(K0+1)] <- expitM(logitp[i,])  
      p[i, (K0+1):K ] <- p[i, (K0+1)]*c(pjack[i], (1-pjack[i])*pAdultChinook)
    }
    logpobs <- log(p)

    ## Adjust lengths for beam spreading:
    ## X is a design matrix that is either distance from shore, or bin id. 
    ## It is expected to contain an intercept term.
    if(adjustLengths){
      L <- dataList$lengthData$L.cm - dataList$Xlength %*% beta
    }else{
      L <- dataList$lengthData$L.cm
    }
    wgts <- dataList$lengthData$weights
    
    ## Prior or penalty terms:
    # ll <- ll + prior_beta(beta) + prior_sigma0(exp(logsigma0)) + prior_sigma(exp(logsigma)) + prior_mu(mu) ## Add prior to q.

    ## Objective function for length based mixture model, conditional on posterior probs.
    objval <- 0
    for( k in 1:K ){
      objval <- objval - sum(wgts*postProb[,k]*(logpobs[,k] + dnorm(x, mu[k], sigma[k], log = TRUE)))
    }
    ## Test fishery component:
    if(includeTestFishery){
      ndays <- nrow(dataList$testFisheryCounts)
      padultchinook <- 1-rowSums(p[,1:(K0+1)])
      nq <- ncol(dataList$testFisheryCounts)
      for( d in 1:ndays ){
        idx <- which(dataList$lengthData$day == d)
        Nd <- sapply(dataList$qSppIndices, FUN = function(x){sum(p[idx,x]*dataList$lengthData$SalmonCount[idx])})
        Ndc <- sum(padultchinook[idx]*dataList$lengthData$SalmonCount[idx])
        Nd <- c(Nd, Ndc)
        if(sum(Nd) == 0) next
        prob <- (q*Nd)/sum(q*Nd)
        objval <- objval - dataList$testFisheryWeights[d]*dmultinom(dataList$testFisheryCounts[d,], prob = prob, size = sum(dataList$testFisheryCounts[d,]), log = TRUE)
      }
    }
    objval
  }
  start <- (lapply(parsOuter, RTMB:::getValues))
  F <- MakeTape(inner_objective_fn, start)
  Newton <- F$newton(1:length(F$par()), maxit = 1000 )
  pars_opt <- Newton(numeric(0))
  c(pars_opt, ll)
}

parsOuter <- speciesComp$parsInit
parsFixed <- speciesComp$parsFixed
dataList <- speciesComp$analysisData
EM <- MakeTape(EMstep, parsOuter)


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
