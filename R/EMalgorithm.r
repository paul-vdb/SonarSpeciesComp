
## Set up an R6 Environment to run the mixture model and keep all the pieces for later.
EMstep <- function(parsOuter){
  getAll(parsOuter, parsFixed, warn = FALSE)
  ll <- 0
 
  ## Parameter Processing
  Kchin <- length(muChin)
  K0 <- length(mu)
  K <- Kchin + K0
  muc <- c(mu, muChin)
  q <- exp(c(logq, 0))

  ## Standard deviation incld. observation error
  logsigma_ <- c(logsigma, logsigmaChin)
  sigma <- sqrt(exp(2*logsigma_) + exp(2*logsigma0))

  ## Set up proportions. alpha parameter for predicting proportions. 
  ## Xalpha is a list of design matrices for each species.
  np <- nrow(dataList$predDF)
  logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
  indx0 <- 1
  for( i in 1:K0 ){
    nc <- ncol(dataList$Xprop[[i]]) 
    indx1 <- indx0 + nc - 1
    logitp[,i] <- as.matrix(dataList$Xprop[[i]]) %*% alpha[indx0:indx1]
    indx0 <- indx1 + 1
  }

  logitpjack <- as.matrix(dataList$XpropChin) %*% alphaJackChinook
  pjack <- 1/(1+exp(-logitpjack))

  p <- matrix(0, nrow = np, ncol = K)
  for( i in 1:np ) {
    p[i, 1:(K0+1)] <- expitM(logitp[i,])  
    p[i, (K0+1):K ] <- p[i, (K0+1)]*c(pjack[i], (1-pjack[i])*pAdultChinook)
  }
  Njc <- sum(p[,K0+1]*dataList$predDF$SalmonCount)
  Nkc <- sum((1-rowSums(p[,1:(K0+1)]))*dataList$predDF$SalmonCount)

  ## Prior jack chinook:
  ll <- ll + dbeta(Njc/(Nkc+Njc), 10, 500, log = TRUE)

  ## Adjust lengths for beam spreading:
  ## X is a design matrix that is either distance from shore, or bin id. 
  ## It is expected to contain an intercept term.
  if(adjustLengths){
    L <- dataList$lengthData$L.cm.adj - as.matrix(dataList$Xlength) %*% beta
  }else{
    L <- dataList$lengthData$L.cm.adj
  }
  wgts <- dataList$lengthData$weights
  
  ## Prior or penalty terms:
  ll <- ll + dprior(sigma = exp(logsigma), sigmaChin = exp(logsigmaChin), sigma0 = exp(logsigma0), mu = mu, muChin = muChin, beta = beta, q = exp(logq), alpha)
  # prior_sigma0(exp(logsigma0)) + prior_beta(beta) # + prior_sigma(exp(logsigma)) + prior_mu(mu) ## Add prior to q.
  
  ## Calculate Posterior Probabilities + log likelihood:
  pobs <- p[dataList$lengthData$grpIndex,]
  outerLL <- calcPostProb(x = L, mu = muc, sigma = sigma, prob = pobs, wgts = wgts)
  postProb <- outerLL$postp
  ll <- ll + outerLL$ll

  ## Test fishery component:
  if(includeTestFishery){
    ndays <- nrow(dataList$testFisheryCounts)
    padultchinook <- 1-rowSums(p[,1:(K0+1)])
    nq <- ncol(dataList$testFisheryCounts)
    for( d in 1:ndays ){
      idx <- which(dataList$predDF$day == d)
      Nd <- sapply(dataList$qSppIndices, FUN = function(x){sum(p[idx,x]*dataList$predDF$SalmonCount[idx])})
      Ndc <- sum(padultchinook[idx]*dataList$predDF$SalmonCount[idx])
      Nd <- c(Nd, Ndc)
      # Nsalmon <- sum(dataList$lengthData$SalmonCount[idx]/dataList$lengthData$nLengths[idx])
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
    muc <- c(mu, muChin)
    q <- exp(c(logq, 0))

    ## Standard deviation incld. observation error
    logsigma_ <- c(logsigma, logsigmaChin)
    sigma <- sqrt(exp(2*logsigma_) + exp(2*logsigma0))

    ## Set up proportions. alpha parameter for predicting proportions. 
    ## Xalpha is a list of design matrices for each species.
    np <- nrow(dataList$predDF)
    logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
    indx0 <- 1
    for( i in 1:K0 ){
      nc <- ncol(dataList$Xprop[[i]]) 
      indx1 <- indx0 + nc - 1
      logitp[,i] <- as.matrix(dataList$Xprop[[i]]) %*% alpha[indx0:indx1]
      indx0 <- indx1 + 1
    }

    logitpjack <- as.matrix(dataList$XpropChin) %*% alphaJackChinook
    pjack <- 1/(1+exp(-logitpjack))

    p <- matrix(0, nrow = np, ncol = K)
    for( i in 1:np ) {
      p[i, 1:(K0+1)] <- expitM(logitp[i,])  
      p[i, (K0+1):K ] <- p[i, (K0+1)]*c(pjack[i], (1-pjack[i])*pAdultChinook)
    }
    logpobs <- log(p[dataList$lengthData$grpIndex,])

    Njc <- sum(p[,K0+1]*dataList$predDF$SalmonCount)
    Nkc <- sum((1-rowSums(p[,1:(K0+1)]))*dataList$predDF$SalmonCount)
    
    ## Objective function for length based mixture model, conditional on posterior probs.
    objval <- 0
    ## Prior jack chinook based on overal proportion.
    objval <- objval - dbeta(Njc/(Nkc+Njc), 10, 500, log = TRUE) 


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
    objval <- objval - dprior(sigma = exp(logsigma), sigmaChin = exp(logsigmaChin), sigma0 = exp(logsigma0), mu = mu, muChin = muChin, beta = beta, q = exp(logq), alpha)
    
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
        Nd <- sapply(dataList$qSppIndices, FUN = function(x){sum(p[idx,x]*dataList$predDF$SalmonCount[idx])})
        Ndc <- sum(padultchinook[idx]*dataList$predDF$SalmonCount[idx])
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
  Newton <- F$newton(1:length(F$par()), maxit = 1000 )
  pars_opt <- Newton(numeric(0))
  c(pars_opt, ll)
}

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
  lli <- -Inf
  lli <- c(lli, vals[npar+1])

  maxit <- speciesComp$optimControl$maxiters
  tol <- speciesComp$optimControl$tolerance
  relativeDiff <- speciesComp$optimControl$relativeDifference
  verbose <- speciesComp$optimControl$verbose

  if(verbose) cat("Running EM Algorithm\n")
  pb <- txtProgressBar(min = 1, max = maxit, initial = 2) 
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
  close(pb)
  if(verbose){ 
    cat("number of iterations =", iter, "\n")
    if(relativeDiff) cat("Convergence was evaluated based on the relative difference of the log likelihood:", diff, ".\n") 
    else cat("Convergence was evaluated based on the difference of the log likelihood:", diff, ".\n") 
  }
  if(!converged) {
    if(maxit < iter) cat("[Warning]  Maximum iterations of", maxit, "were reached without convergence to a difference of", tol, "between the log likelihood.\n")
    if(maxit >= iter) cat("[Warning]  Failed to converge with a tolerance of", tol, "between the log likelihood.\n")
    cat("If a (relative) difference of", diff, "seems to be close enough then you may consider ignoring this message.\n")
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
    p[i, (K0+1):K ] <- p[i, (2+1)]*c(pjack[i], (1-pjack[i])*parsFit$pAdultChinook)    
    N[i,] <- p[i,] * speciesComp$analysisData$predDF$SalmonCount[i]
  }
  ## A lot of work to not have a tidy dependence...
  form <- paste0("cbind(", paste(speciesComp$species, collapse = ","), ") ~ SonarBank + SonarBin + SonarAim + day")
  colnames(N) <- speciesComp$species
  estimatedProp <- cbind(speciesComp$analysisData$predDF, N)
  estimatedProp <- estimatedProp |> aggregate(as.formula(form), sum)
  estimatedProp[, speciesComp$species] <- t(apply(estimatedProp[, speciesComp$species], 1, FUN = function(x){x/sum(x)}))
  speciesComp$estimatedDailyProportions <- estimatedProp
  colnames(p) <- speciesComp$species
  hourlyPredDF <- cbind(speciesComp$analysisData$predDF[, c("Date", "day", "SonarBank", "SonarAim", "SonarBin", "Hour", "HourOrder", "SalmonCount")], p)
  speciesComp$estimatedHourlyProportions <- hourlyPredDF
}
