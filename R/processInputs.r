#' Species Composition Estimation Methods
#'
#' Object to step through all the different setups and steps of estimating species composition on the Fraser River (Qualark and Mission) using Hydroacoustics and the test fishery.
#'
#' @param species vector of species names  e.g. c("largeresident", "jackChinook", "sockeye", "adultchinook")
#' @param estDate Date that you will estimate. Not necessary to set yet.
#' @param ndays Number of days to combine, defaults to 1
#' @param site Choose between "Mission" or "Qualark", the difference being the roll angle.
#'
#' @details This is an R6 object that holds all of the methods to estimate species composition and generate reports. The methods include:
#' \itemize{
#'  \item{setSpecies}{Update the object with species to fit in a model. Input: 'species': options are c("smallresident", "largeresident", "pink", "sockeye", "jackchinook", "adultchinook", "smalladultchinook", "largeadultchinook").}
#'  \item{setDate}{Update the estimation date and number of days to use. Input: 'estDate' and 'ndays'}
#'  \item{controlOptimization}{Control optimization of the EM algorithm. Input: 'tolerance' what difference is okay (default 1e-8), 'maxiters' maximum iterations of EM algorithm (default 1000), 
#'  'relativeDifference' true or false to choose relative vs absolute difference on each iteration of algorithm, 'verbose' true or false to print details of EM algorithm.}
#'  \item{processData}{Core function for adding data to the object. Does the main data processing and adds a few columns needed elsewhere for the entire dataset added to make it fast and easy to do different days. 
#'  Input: 'sonarCounts' and 'sonarLengths' are data frames see \code{?processMissionLengths} or \code{?processQualarkLengths} for details, 
#'  'testFisheryCounts' data frame see \code{?processWhonnockCounts} or \code{?processQualarkTFCounts} for details, 
#'  'includeTangled' Whether or not to include tangled or just gilled counts, and
#'  'dropN' Keep length measurements per hour/bin/bank etc. if > dropN (default 2).}
#'  \item{setSpeciesLengths}{mu = NULL, sigma = NULL, pChin = NULL, testFisheryLengths = NULL, ndays = 6}
#'  \item{setLengthAdjustment}{Choose whether to do length adjustments and the formula to use. Input: 'formula' (default ~ beamWidth.cm) but may consider (~ -1 + SonarBin), 'adjustLengths' (default TRUE).}
#'  \item{setModelParameters}{Core function to set the model before fitting with EM algorithm, especially which values to fix and which to fit. Input: 
#'  'fixedParameters' (default c("mu", "sigma", "muChinook", "sigmaChinook", "pJackChinook")),
#'  'parameterValues' used as either initial values or fixed values, 'adjustLengths' should lengths be adjusted (default NULL), 
#'  'includeTestFishery' Should test fishery data be used (default NULL), 'testFisheryWeights' How much to weight test fishery counts (default 1). NULL implies existing values if already set in object.}
#'  \item{setModelProportions}{Function to set what the relationship with time, date, formula = list()}
#'  \item{fitModel}{Runs EM algorithm based on object values for input. Results are saved as part of object: 'estimatedParameters', 'estimatedHourlyProportions', and 'estimatedDailyProportions'}
#'  \item{simulate}{Simulate the model for the counts on the 'estDate' and 'ndays' with parameter values equal to 'estimatedParameters' within the object. See \code{?simulateGivenData} for details.}
#'  \item{setPriors}{Can set penalties/priors for terms at their real scale. Input: 'priors' List of functions, named after the parameter the prior will set (default list()).}
#' }
#' To get more detailed information about each method, do \code{?method} to see full description.
#'
#' @returns An R6 object with methods and values used for fitting species composition models.
#'
#' @examples
#' speciesComp <- speciesCompSummary$new(species = c("largeresident", "jackchinook", "sockeye", "adultchinook"), site = "Mission", estDate = "2023-08-05")
#' 
#' @export
speciesCompSummary <- R6Class("SpeciesComp",
  public = list(
    site = NULL,
    species = NULL, 
    Species = NULL,
    estDate = NULL,
    ndays = NULL,
    feasibleLengths = c(10, 120),
    includeTestFishery = FALSE,
    adjustLengths = TRUE,
    rollAngle = 0,  ## *** Make sure to set this
    sonarLengths = NULL,
    testFisheryCounts = NULL,
    salmonCounts = NULL,
    speciesPriorLengths = NULL,
    analysisData = list(),
    parameters = list(),
    parsInit = list(),
    parsFixed = list(),
    estimatedParameters = NULL,
    estimatedDailyProportions = NULL,
    estimatedHourlyProportions = NULL,
    simData = NULL,
    priorDist = NULL,
    optimControl = list(tolerance = 1e-8, maxiters = 1000, relativeDifference = FALSE, verbose = FALSE),    
    
    initialize = function(species = c("jackChinook", "sockeye", "adultchinook"), estDate = NULL, ndays = 1, site = NULL) {
      if(missing(site)) stop("Need to initialize with a site = 'Mission' or 'Qualark'")
      if(!site %in% c("Mission", "Qualark")) stop("Need to initialize with a site = 'Mission' or 'Qualark'")
      self$species <- speciesCheck(species)
      self$Species <- speciesLabels(self$species)
      self$estDate <- checkDate(estDate)
      self$ndays <- ndays
      self$site <- site
      if(site == "Mission") self$rollAngle <- 0     
      else self$rollAngle <- 35
      self$priorDist <- function(...){0}  
    },
    ## Set the species names: Convenient and for printing.
    setSpecies = function(species = NULL){
      if(!is.null(species)){
        self$species <- speciesCheck(species)
        self$Species <- speciesLabels(species)
      }
    },
    ## set date:
    setDate = function(estDate = NULL, ndays = NULL){
      if(!is.null(estDate)){
        self$estDate <- checkDate(estDate)
      }
      if(!is.null(ndays)) self$ndays <- ndays
    },
    controlOptimization = function(tolerance = NULL, maxiters = NULL, relativeDifference = NULL, verbose = NULL){
      if(!is.null(tolerance)) self$optimControl$tolerance <- tolerance
      if(!is.null(maxiters)) self$optimControl$maxiters <- maxiters
      if(!is.null(relativeDifference)) self$optimControl$relativeDifference<- relativeDifference
      if(!is.null(verbose)) self$optimControl$verbose <- verbose
    },
    processData = function(sonarCounts, sonarLengths, testFisheryCounts = NULL, includeTangled = TRUE, dropN = 2, salmonPassageTable = NULL, feasibleLengths = NULL){
      if(!is.null(feasibleLengths)) self$feasibleLengths <- feasibleLengths
      if(!is.null(testFisheryCounts))  self$includeTestFishery <- TRUE

      if(self$site == "Qualark") {
        processQualarkLengths(self, sonarcounts, sonarlengths)
      }
      if(self$site == "Mission") {
        processMissionLengths(self, sonarCounts, sonarLengths, dropN)
        processWhonnockCounts(self, testFisheryCounts, tangled = includeTangled)
        if(!is.null(salmonPassageTable)) processMissionSalmonPassage(self, salmonPassageTable)
      }
    },
    setSpeciesLengths = function(mu = NULL, sigma = NULL, pChin = NULL, testFisheryLengths = NULL, ndays = 6){
      setSpeciesLengths(self, mu, sigma, pChin, testFisheryLengths, ndays)
    },
    setLengthAdjustment = function(formula = ~ beamWidth.cm, adjustLengths = TRUE){
      setLengthAdjustment(self, formula, adjustLengths)
    },
    setModelParameters = function(fixedParameters = c("mu", "sigma", "muChinook", "sigmaChinook", "pJackChinook"), 
                                  parameterValues = list(), adjustLengths = NULL, includeTestFishery = NULL, testFisheryWeights = 1){
      setModelParameters(self, fixedParameters, parameterValues, adjustLengths, includeTestFishery, testFisheryWeights)
    },
    setModelProportions = function(formula = list()){
      setModelProportions(self, formula)
    },
    fitModel = function(simulatedData = FALSE){
      runEMAlgorithm(self, simulatedData)
    },
    simulate = function(){
      simulateGivenData(self)
    },
    setPriors = function(priors = list()){
      setPriors(self, priors = priors)
    },
    predictSpeciesComposition = function(){
      predictSpeciesComposition(self)
    },
    plot = function(day = 1, includeProportionLabels = FALSE, SonarBin = NULL, SonarAim = NULL, SonarBank = NULL, ...){
      plotmix(self, day, includeProportionLabels, SonarBin, SonarAim, SonarBank, ...)
    }    
  )
)

#' @export
setPriors <- function(speciesComp, priors = list()){        
        default <- function(x){0}
        dsigma <- extractControls(priors$sigma, default)
        dsigmaChin <- extractControls(priors$sigmaChinook, default)
        dsigma0 <- extractControls(priors$sigma0, default)
        dmu <- extractControls(priors$mu, default)
        dmuChin <- extractControls(priors$muChinook, default)
        dbeta <- extractControls(priors$beta, default)
        dq <- extractControls(priors$catchability, default)
        dalpha <- extractControls(priors$alpha, default)
        dalphaJackChinook <- extractControls(priors$alphaJackChinook, default)
        ## Set priors into species comp object:
        speciesComp$priorDist <- function(sigma, sigmaChin, sigma0, mu, muChin, beta, q, alpha, alphaJackChinook){
          dsigma(sigma) + dsigmaChin(sigmaChin) + dsigma0(sigma0) + dmu(mu) + dmuChin(muChin) + dbeta(beta) + dq(q) + dalpha(alpha) + dalphaJackChinook(alphaJackChinook)
        }
}


speciesCheck <- function(spp){
  allowed <- c("smallresident", "largeresident", "pink", "sockeye", "jackchinook", "adultchinook", "smalladultchinook", "largeadultchinook")
  spp <- tolower(gsub(" ", "", spp))
  if(any(!spp %in% allowed)) stop(paste("Species currently allowed are", allowed))
  if(any(spp == "adultchinook") & any(spp %in% c("smalladultchinook", "largeadultchinook"))){
    spp <- spp[-which(spp == "adultchinook")]
  }
  return(allowed[which(allowed %in% spp)])  ## Order them correctly.
}

speciesLabels <- function(species){
  match <- c("smallresident" = "Small Resident", "largeresident" = "Large Resident", "jackchinook" = "Chinook Jack", 
             "pink" = "Pink", "sockeye" = "Sockeye", "adultchinook" = "Chinook Adult", "smalladultchinook" = "Chinook Small Adult", "largeadultchinook" = "Chinook Large Adult")
  as.character(as.character(match[species]))
}

speciesColours <- function(species){
  cols <- c("smallresident" = "#FFB100", "largeresident" = "#656837",  "jackchinook" = "#A33CC7", "pink" = "#FF8DA1", "sockeye" = "#CD0000", 
    "adultchinook" = "#27408B", "smalladultchinook" = "#27408B", "largeadultchinook" = "#27408B")
  cols[species]
}

speciesOrder <- function(){
  c("smallresident", "largeresident", "jackchinook", "pink", "sockeye", "adultchinook")
}

checkDate <- function(eDate){
  if(is.null(eDate)){return(NULL)}
  if(!is.Date(eDate)){
    eDate <- as.Date(eDate, "%Y-%m-%d")
    if(is.na(eDate)) stop("Date must be submitted as 'Y-m-d', e.g. '2018-08-05'")
  }
  eDate
}

## Core Mission Data Processing function:
processMissionLengths <- function(speciesComp, sonarCounts, sonarLengths, dropN = 2){
  
  year <- as.numeric(format(speciesComp$estDate, "%Y"))
  
  ## Check names in dataframes are present:
  countColNames <- c("Date", "MissionDate", "Hour", "MinsCounted", "SonarCode", "SonarAim", "SonarBin", 
                     "Up.0to5", "Dn.0to5", "Up.5to10", "Dn.5to10")
  lengthColNames <- c("Date", "MissionDate", "Hour", "SonarCode", "SonarAim", "SonarBin", 
                 "R.m", "L.cm", "Dir", "EnterTheta.deg")
                 
  namecheck <- which(!countColNames %in% names(sonarCounts))
  if(length(namecheck)) stop(paste0("We require column names in sonarCounts: ", paste(countColNames[namecheck], collapse = ",")))
  namecheck <- which(!lengthColNames %in% names(sonarLengths))
  if(length(namecheck)) stop(paste0("We require column names in sonarLengths: ", paste(lengthColNames[namecheck], collapse = ",")))

  ## Clean up length data:
  if(!is.Date(sonarLengths$Date)) sonarLengths$Date <- as.Date(sonarLengths$Date, format = "%m/%d/%Y")
  if(!is.Date(sonarLengths$MissionDate)) sonarLengths$MissionDate <- as.Date(sonarLengths$MissionDate, format = "%m/%d/%Y")
  ## Check Mission Date:
  MissionDateCheck <- sonarLengths |> subset(Date == MissionDate & Hour < 5)
  if(nrow(MissionDateCheck) > 0) stop("Invalid Mission Dates Detected in sonarLengths")

  ## Subset Count data to match length data:
  sonarLengths <- sonarLengths |> subset(as.numeric(format(Date, "%Y")) == year)
  sonarLengths <- sonarLengths |> subset(SonarBin %in% c("Bin1", "Bin2", "Bin3") & SonarCode %in% c("A1", "A2") & L.cm >= speciesComp$feasibleLengths[1] & L.cm <= speciesComp$feasibleLengths[2])
  ## Create a look up code for fast merging:
  sonarLengths <- sonarLengths |> within(join_id <- factor(paste(SonarBin, SonarAim, SonarCode, Hour, MissionDate, sep = ""))) |>
                                  within(lookup_code <- as.numeric(join_id))
  sonarLengthsN <- sonarLengths |> aggregate( L.cm ~ lookup_code, FUN = length)
  names(sonarLengthsN)[names(sonarLengthsN) == 'L.cm'] <- 'nLengths'
  sonarLengths <- sonarLengths |> merge(sonarLengthsN, all.x = TRUE, all.y = TRUE)
  
  ## Subset Count data to match length data:
  sonarCounts <- sonarCounts |> subset(SonarBin %in% c("Bin1", "Bin2", "Bin3") & SonarCode %in% c("A1", "A2") & Up.0to5 != "" & MinsCounted > 0)
  suppressWarnings(sonarCounts <- sonarCounts |> within(Up <- ifelse(MinsCounted == 5, as.numeric(Up.0to5), as.numeric(Up.0to5) + as.numeric(Up.5to10))))
  suppressWarnings(sonarCounts <- sonarCounts |> within(Down <- ifelse(MinsCounted == 5, as.numeric(Dn.0to5), as.numeric(Dn.0to5) + as.numeric(Dn.5to10))))
  sonarCounts <- sonarCounts |> within(MinsCounted <- ifelse(is.na(Up.5to10), MinsCounted, 10))
  sonarCounts <- sonarCounts |> within(Count <- Up + Down) |> within(SalmonCount <- Up - Down) |> subset(!is.na(Count))
  
  ## check Date and Mission Date: Format should be m/d/Y.
  if(!is.Date(sonarCounts$Date)) sonarCounts$Date <- as.Date(sonarCounts$Date, format = "%m/%d/%Y")
  if(!is.Date(sonarCounts$MissionDate)) sonarCounts$MissionDate <- as.Date(sonarCounts$MissionDate, format = "%m/%d/%Y")

  ## Check Mission Date:
  MissionDateCheck <- sonarCounts |> subset(Date == MissionDate & Hour < 5)
  if(nrow(MissionDateCheck) > 0) stop("Invalid Mission Dates Detected")

  suppressMessages(
    sonarCounts <- sonarCounts |> within(join_id <- factor(paste(SonarBin, SonarAim, SonarCode, Hour, MissionDate, sep = ""), levels = levels(sonarLengths$join_id))) |>
                   subset(!is.na(join_id)) |> within(lookup_code <- as.numeric(join_id))
  )
  ## Merge count data to lengths to get correction to est proportion with test fishery count and stratum weights:
  sonarLengths <- sonarLengths |> merge(sonarCounts[, c("lookup_code", "Count", "SalmonCount", "MinsCounted")])
  sonarLengths <- sonarLengths |> within(weights <- (Count/MinsCounted*15)/nLengths)  ## Assuming 15 min file.
  
  minhour <- min(sonarCounts$Hour[sonarCounts$MissionDate == sonarCounts$Date], na.rm = TRUE)

  ## Account for "Mission Date" and then just use "Date".
  sonarLengths$HourOrder <- sonarLengths$Hour - minhour
  sonarLengths$HourOrder[sonarLengths$Hour < minhour] <- 24 + sonarLengths$HourOrder[sonarLengths$Hour < minhour]

  sonarLengths$Date <- sonarLengths$MissionDate ## Rename MissionDate to Date for convenience.
  
  ## Add Stratum:
  sonarLengths <- sonarLengths |> within(SonarBank <- ifelse(SonarCode == "A1", "Left Bank", "Right Bank"))  
  sonarLengths <- sonarLengths |> within(stratum <- paste(SonarBank, SonarBin, SonarAim, sep = "_"))
  ## Add beam width:
  sonarLengths <- sonarLengths |> within(beamWidth.cm <- R.m*0.3*pi/180*100)

  ## Will remove single lengths measured in an hour:
  sonarLengths <- sonarLengths |> filter(nLengths > dropN)

  ## I don't like SonarAim being called 3 on right bank if only one aim.
  sonarLengths <- sonarLengths |> within(SonarAimF <- ifelse(SonarBank == "Right Bank", "Aim3", SonarAim))

  speciesComp$sonarLengths <- sonarLengths
}

## Process the Whonnock Test Fishery Counts: 
## Include the test fishery counts in format as per PSC and let them choose to include tangled or only gilled.
processWhonnockCounts <- function(speciesComp, testFisheryCounts, tangled = TRUE){

  if(!is.POSIXct(testFisheryCounts$setStart)){
    stop("We require that 'setStart' be POSIXct.")
  }

  out <- testFisheryCounts[,c("TRIP_DTT", "setStart", "nSets")]
  names(out)[1] <- "Date"
  out$Date <- as.Date(out$Date)
  
  tmp <- data.frame("pink" = testFisheryCounts$`Pink All Gilled` + tangled*testFisheryCounts$`Pink All Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("sockeye" = testFisheryCounts$`Sockeye Adult Gilled` + tangled*testFisheryCounts$`Sockeye Adult Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("jackchinook" = testFisheryCounts$`Chinook Jack Gilled` + tangled*testFisheryCounts$`Chinook Jack Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("adultchinook" = testFisheryCounts$`Chinook Adult Gilled` + tangled*testFisheryCounts$`Chinook Adult Tangled`)
  out <- cbind(out, tmp)

  speciesComp$testFisheryCounts <- out
}

processMissionSalmonPassage = function(speciesComp, salmonPassageTable){
  ## Check names in dataframes are present:
  colNames <- c("SelectedMethod", "LB_BIN1_AIM1", "LB_BIN1_AIM2", "LB_BIN1", "LB_BIN2", 
                     "LB_BIN3", "LB_BIN4", "LB_BIN5", "LB_BIN6", "Offshore", "RB_BIN4", "RB_BIN3",
                     "RB_BIN2", "RB_BIN1", "MissionDate")

  namecheck <- which(!colNames %in% names(salmonPassageTable))
  if(length(namecheck)) stop(paste0("We require column names in salmonPassageTable: ", paste(colNames[namecheck], collapse = ",")))

  namedf <- data.frame(grp = c("LB_BIN1_AIM1", "LB_BIN1_AIM2", "LB_BIN1", "LB_BIN2", "LB_BIN3", "LB_BIN4", "LB_BIN5", 
                               "LB_BIN6", "Offshore", "RB_BIN4", "RB_BIN3", "RB_BIN2", "RB_BIN1"),
                      SonarBank = c(rep("Left Bank",8), "Offshore", rep("Right Bank", 4)),
                      SonarCode = c(rep("A1",8), "mobile", rep("A2", 4)),
                      SonarBin = c("Bin1", "Bin1", "Bin1", "Bin2", "Bin3", "Bin4", "Bin5", "Bin6", NA, "Bin4", "Bin3", "Bin2", "Bin1"),
                      SonarAim = c("Aim1", "Aim2","Aim1+Aim2", "Aim3", "Aim3", "Aim3", "Aim3", "Aim3", NA, "Aim1", "Aim1", "Aim1", "Aim1")
                      )

  ## Clean up length data:
  if(!is.Date(salmonPassageTable$MissionDate)) salmonPassageTable$MissionDate <- as.Date(salmonPassageTable$MissionDate, format = "%m/%d/%Y")
  salmonPassageTable$Date <- salmonPassageTable$MissionDate

  longSalmon <- do.call("rbind", lapply(namedf$grp, FUN = function(x){
    tab = cbind(salmonPassageTable[, c("SelectedMethod", "Date", "MissionDate", x)], "grp" = x); colnames(tab)[4] <-  "count"; tab
    }))
  longSalmon <- longSalmon |> subset(!is.na(count))
  longSalmon <- longSalmon |> merge(namedf)
  speciesComp$salmonCounts <- longSalmon
}

## Add Roxygen
setModelProportions <- function(speciesComp, formula = list()){
  pFormula <- extractControls(formula$all, ~ -1 + stratum:day)

  ## Model Data:
  ndays <- speciesComp$ndays
  dates_ <- seq(speciesComp$estDate - ndays + 1, speciesComp$estDate, 1)
  lengthData <- speciesComp$sonarLengths |> subset(Date %in% dates_) |> within(day <- as.numeric(factor(Date)))
  lengthData <- lengthData |> within(L.cm.adj <- L.cm/cos(speciesComp$rollAngle/180*pi))
  speciesComp$analysisData$lengthData <- lengthData
  
  ## Subset per stratum/Hour to be able to efficiently compute quantities:
  speciesComp$analysisData$lengthData <- speciesComp$analysisData$lengthData |> 
                                          within(grpID <- factor(paste(day, Hour, stratum, sep = "_"))) |>
                                          within(grpIndex <- as.numeric(grpID))

  predDF <- speciesComp$analysisData$lengthData %>% subset(!duplicated(grpIndex))
  predDF <-  predDF[order( predDF$grpIndex),]
  speciesComp$analysisData$predDF <- predDF
  
  ## Potential species models:
  sppNames <- grep("chinook", speciesComp$species, invert = TRUE, value = TRUE)
  speciesComp$analysisData$Xprop <- list()
  for( i in seq_along(sppNames) ){
    sppi <- sppNames[i]
    formulai <- extractControls(formula[[sppi]], pFormula)
    speciesComp$analysisData$Xprop[[sppi]] <- model.matrix(formulai, data =  predDF)
  }
  ## jackChinook model:
  formulajc <- extractControls(formula[["jackchinook"]], ~ 1) ## adult chinook or small/large adults will default to the same, minus input.
  speciesComp$analysisData$XpropChin <- model.matrix(formulajc, data = predDF)

  ## Add in test fishery counts if present:
  if(!is.null(speciesComp$testFisheryCounts)){
    speciesComp$analysisData$testFisheryCounts <- speciesComp$testFisheryCounts |> 
                                                  subset(as.Date(Date) %in% dates_) |> 
                                                  subset(select = c(pink, sockeye, jackchinook, adultchinook))
    if(all(speciesComp$species != "pink")) 
      speciesComp$analysisData$testFisheryCounts <- speciesComp$analysisData$testFisheryCounts |> 
                                                    subset(select = c(-pink))
    if(all(speciesComp$species != "sockeye")) 
      speciesComp$analysisData$testFisheryCounts <- speciesComp$analysisData$testFisheryCounts |> 
                                                    subset(select = c(-sockeye))
    if(!any(grepl("chinook", speciesComp$species))) 
      speciesComp$analysisData$testFisheryCounts <- speciesComp$analysisData$testFisheryCounts |> 
                                                    subset(select = c(-jackchinook, -adultchinook))
  }
  speciesComp$analysisData$testFisheryCounts <- as.matrix(speciesComp$analysisData$testFisheryCounts)
}

setLengthAdjustment <- function(speciesComp, formula = ~ beamWidth.cm, adjustLengths = TRUE){
  ## Beam Width Adjustment or potentially just bins, or nothing (adjustLengths = FALSE)
  speciesComp$analysisData$Xlength <- model.matrix(formula, data = speciesComp$analysisData$lengthData)
  speciesComp$adjustLengths <- adjustLengths
}

setSpeciesLengths <- function(speciesComp, mu = NULL, sigma = NULL, pChin = NULL, testFisheryLengths = NULL, ndays = 7){

  dates_ <- seq(speciesComp$estDate - ndays + 1, speciesComp$estDate, 1)

  ## Default values:
  mu_ <- c("smallresident" = 22.7, "largeresident" = 39, "jackchinook" = 43, "pink" = 50, "sockeye" = 58, "adultchinook" = 75, "smalladultchinook" = 63, "largeadultchinook" = 76)
  sigma_ <- c("smallresident" = 2.5, "largeresident" = 2.5, "jackchinook" = 2.5, "pink" = 3, "sockeye" = 3.5, "adultchinook" = 7.5,"smalladultchinook" = 5, "largeadultchinook" = 6)
  pChin_ <- c("jackchinook" = 0.05, "adultchinook" = 0.95,"smalladultchinook" = 0.17, "largeadultchinook" = 0.78)
  
  ## Subset for current species
  muChin_ <- mu_[grep("chinook", speciesComp$species, value = TRUE)]
  sigmaChin_ <- sigma_[grep("chinook", speciesComp$species, value = TRUE)]
  mu_ <- mu_[grep("chinook", speciesComp$species, value = TRUE, invert = TRUE)]
  sigma_ <- sigma_[grep("chinook", speciesComp$species, value = TRUE, invert = TRUE)]
  pChin_ <- pChin_[grep("chinook", speciesComp$species, value = TRUE)]

  ## Add in user supplied values:
  mu_[grep("chinook", names(mu), value = TRUE, invert = TRUE)] <- mu[grep("chinook", names(mu), value = TRUE, invert = TRUE)]
  sigma_[grep("chinook", names(sigma), value = TRUE, invert = TRUE)] <- sigma[grep("chinook", names(mu), value = TRUE, invert = TRUE)]  
  muChin_[grep("chinook", names(mu), value = TRUE)] <- mu[grep("chinook", names(mu), value = TRUE)]
  sigmaChin_[grep("chinook", names(sigma), value = TRUE)] <- sigma[grep("chinook", names(mu), value = TRUE)]  
  pChin_[grep("chinook", names(pChin), value = TRUE)] <- pChin

  ## Now if test fishery lengths are present use those to infill mu and sigma:
  if(!is.null(testFisheryLengths)){
    if(!all(c("FL.cm", "Date", "Species") %in% names(testFisheryLengths)))
      stop("TestFisheryLengths dataframe must have 'Date', 'FL.cm', and 'Species'.")

    if(!is.POSIXct(testFisheryLengths$Date) & !is.Date(testFisheryLengths$Date))
      stop("Date must be provided in a standard date format")

    testFisheryLengths$Date <- as.Date(testFisheryLengths$Date)
  
    sppNames <- setdiff(grep("resident", speciesComp$species, invert = TRUE, value = TRUE), names(mu))
    sppchin <- grep("chinook", sppNames)
    if(length(sppchin) > 0){
      sppNames <- sppNames[-sppchin]
      sppNames <- c(sppNames, "chinook")
    }
    tfdays <- testFisheryLengths |> subset(Date %in% dates_)
    for(i in seq_along(sppNames) ){
      tfl <- tfdays |> subset(grepl(sppNames[i], tolower(Species)) & !is.na(FL.cm))
      x <- tfl$FL.cm
      if(length(x) < 7) {
        cat("[Warning]  Unable to estimate", sppNames[i] , "lengths from the test fishery data. Using default mean and variance.\n")        
        next
      }
      
      if("chinook" == sppNames[i]){
        nchin <- length(pChin_)
        capture.output(fit <- mixtools::normalmixEM(x, lambda = pChin_[-nchin], mu = muChin_, sigma = sigmaChin_))
        ## Check that we captured jack chinook: 
        if(min(fit$mu) > 50){
           ## Make it more specific about how to actually set alternative.
          if(!"jackChinook" %in% names(mu)) cat("[Warning]  Unable to estimate jack Chinook from the test fishery data. Using default mean and variance.\n")
          if(any(speciesComp$species == "largeadultchinook")){
            capture.output(fit <- mixtools::normalmixEM(x[x > 50], lambda = pChin_[2], mu = muChin_[-1], sigma = sigmaChin_[-1]))
            ord <- order(fit$mu)
            muChin_[c("smalladultchinook", "largeadultchinook")] <- fit$mu[ord]
            sigmaChin_[c("smalladultchinook", "largeadultchinook")] <- fit$sigma[ord]
            pChin_[c("smalladultchinook", "largeadultchinook")] <- (1-pChin_["jackchinook"]) * (fit$lambda[ord])
          }else{
            muChin_["adultchinook"] <- mean(x[x >= 50])
            sigmaChin_["adultchinook"] <- sd(x[x >= 50])
          }
        }else{
          ord <- order(fit$mu)
          muChin_ <- fit$mu[ord]
          sigmaChin_ <- fit$sigma[ord]
          pChin_ <- fit$lambda[ord]
        }
        ## If species are provided overwrite everything you just did.
        muChin_[grep("chinook", names(mu), value = TRUE)] <- mu[grep("chinook", names(mu), value = TRUE)]
        sigmaChin_[grep("chinook", names(sigma), value = TRUE)] <- sigma[grep("chinook", names(sigma), value = TRUE)]
        pChin_[grep("chinook", names(pChin), value = TRUE)] <- pChin
        pChin_ <- pChin_/sum(pChin_)
      }else{
        mu_[sppNames[i]] <- mean(x)
        sigma_[sppNames[i]] <- sd(x)
      }
    }
  }
  
  ## Add the values to the the parameters:
  speciesComp$parameters$muChin <- muChin_
  speciesComp$parameters$sigmaChin <- sigmaChin_
  speciesComp$parameters$pChin <- pChin_
  speciesComp$parameters$mu <- mu_
  speciesComp$parameters$sigma <- sigma_
}

setModelParameters <- function(speciesComp, fixedParameters = c("mu", "sigma", "muChinook", "sigmaChinook", "pJackChinook"), 
                          parameterValues = list(), adjustLengths = NULL, includeTestFishery = NULL, testFisheryWeights = 1){
  ## Another place to set length adjustment:
  if(!is.null(adjustLengths)) speciesComp$adjustLengths <- adjustLengths
  if(!is.null(includeTestFishery)) speciesComp$includeTestFishery <- includeTestFishery

  if(speciesComp$site == "Mission") { 
    q <- c( "sockeye" = 0.33, "pink" = 0.16, "jackchinook" = 0.5)  ## From Yunbo's paper.
  }else{
    q <- c("sockeye" = 0.25, "pink" = 0.10, "jackchinook" = 0.5)
  }
  q <- q[names(q) %in% speciesComp$species]

  speciesComp$analysisData$qSppNames <- names(q)
  speciesComp$analysisData$qSppIndices <- which(speciesComp$species %in% names(q))

  ## Extract all the fixed values present and use the parameter values already present otherwise.
  ## These should be setup in via 'setSpeciesLengths'.
  mu <- as.numeric(extractControls(parameterValues$mu, speciesComp$parameters$mu))
  sigma <- as.numeric(extractControls(parameterValues[["sigma"]], speciesComp$parameters$sigma))
  muChin <- as.numeric(extractControls(parameterValues$muChinook, speciesComp$parameters$muChin))
  sigmaChin <- as.numeric(extractControls(parameterValues$sigmaChinool, speciesComp$parameters$sigmaChin))
  pChin <- extractControls(parameterValues$pChinook, speciesComp$parameters$pChin)
  ## If setLengthAdjustment() has not been called, call it now:
  if(is.null(speciesComp$analysisData$Xlength)) setLengthAdjustment(speciesComp)
  beta <- as.numeric(extractControls(parameterValues$beta, numeric(ncol(speciesComp$analysisData$Xlength))))

  if(length(beta) != ncol(speciesComp$analysisData$Xlength)){
    stop("Provided values for beta parameter do not match the number of columns of the design matrix from 'setLengthAdjustment'")  
  }
  sigma0 <- as.numeric(extractControls(parameterValues$sigma0, 6))
  q <- as.numeric(extractControls(parameterValues$catchability, q))

  ## Make sure all the values are provided:
  if(length(speciesComp$species) != length(c(mu, muChin))){
    stop("Must provide values for mean lengths of all species via `setSpeciesLengths` or via parameterValues input.")
  }
  if(length(speciesComp$species) != length(c(sigma, sigmaChin))){
    stop("Must provide values for standard deviation lengths of all species via `setSpeciesLengths` or via parameterValues input.")
  }
  if(length(grep("chinook", speciesComp$species, value = TRUE)) != length(pChin)){
    stop("Must provide values mean lengths of all species via `setSpeciesLengths` or via parameterValues input.")
  }

  speciesComp$parsInit <- list()
  speciesComp$parsFixed <- list()

  if("mu" %in% fixedParameters){
    speciesComp$parsFixed$mu <- mu
  }else{
    speciesComp$parsInit$mu <- mu
  }
  if("sigma" %in% fixedParameters){
    speciesComp$parsFixed$logsigma <- log(sigma)
  }else{
    speciesComp$parsInit$logsigma <- log(sigma)
  }
  if("muChinook" %in% fixedParameters){
    speciesComp$parsFixed$muChin <- muChin
  }else{
    speciesComp$parsInit$muChin <- muChin
  }
  if("sigmaChinook" %in% fixedParameters){
    speciesComp$parsFixed$logsigmaChin <- log(sigmaChin)
  }else{
    speciesComp$parsInit$logsigmaChin <- log(sigmaChin)
  }
  if("beta" %in% fixedParameters){
    speciesComp$parsFixed$beta <- beta
  }else{
    speciesComp$parsInit$beta <- beta
  }
  if(!speciesComp$adjustLengths){
    speciesComp$parsInit$beta <- NULL
  }
  if("sigma0" %in% fixedParameters){
    speciesComp$parsFixed$logsigma0 <- log(sigma0)
  }else{
    speciesComp$parsInit$logsigma0 <- log(sigma0)
  }
  if("catchability" %in% fixedParameters){
    speciesComp$parsFixed$logq <- log(q)
  }else{
    speciesComp$parsInit$logq <- log(q)
  }
  if(!speciesComp$includeTestFishery){
    speciesComp$parsInit$logq <- NULL
    speciesComp$parsFixed$logq <- log(q)
  }

  nc <- ncol(speciesComp$analysisData$XpropChin)
  logitpChin <- logitM(pChin)
  nr <- length(logitpChin)
  if("pJackChinook" %in% fixedParameters){
    speciesComp$parsFixed$alphaJackChinook <- rep(logitpChin[1], nc)
  }else{
    speciesComp$parsInit$alphaJackChinook <- rep(logitpChin[1], nc)
  }
  ## Fixing probability of adult chinook:
  speciesComp$parsFixed$pAdultChinook <- pChin[-1]/(1-pChin[1])

  ## Initialize proportions:
  probs <- sapply(1:length(mu), FUN = function(x){dnorm(speciesComp$analysisData$lengthData$L.cm, mu[x], sqrt(sigma[x]^2+sigma0^2))})
  probs <- cbind(probs, sapply(speciesComp$analysisData$lengthData$L.cm, FUN = function(x){sum(pChin*dnorm(x, muChin, sqrt(sigmaChin^2+sigma0^2)))}))
  probid <- apply(probs, 1, which.max)
  alpha <- NULL
  for( i in 1:length(mu) ){
    y <- (probid == i)*1
    X <- speciesComp$analysisData$Xprop[[i]][speciesComp$analysisData$lengthData$grpIndex,]
    if(sum(y) > 3){
      fit <- glm(cbind(y, 1)~-1+X, family = "binomial")
      alphai <- coef(fit)
      alphai[is.na(alphai)] <- -2
    }else{
      alphai <- rep(-2, ncol(X))
    }
    alpha <- c(alpha, as.numeric(alphai))
  }
  speciesComp$parsInit$alpha <- alpha
  if(length(testFisheryWeights) == 1) testFisheryWeights <- rep(testFisheryWeights, speciesComp$ndays)
  speciesComp$analysisData$testFisheryWeights <- testFisheryWeights
}

predictSpeciesComposition = function(speciesComp, proportionOffshore = NULL, salmonPassageTable = NULL){
  if(is.null(speciesComp$salmonCounts) ){
    if(is.null(salmonPassageTable) & speciesComp$site == "Mission") stop("A salmonPassageTable needs to be provided to compute species composition at Mission.")
    if(speciesComp$site == "Mission") processMissionSalmonPassage(salmonPassageTable)
  }
  if(!any(speciesComp$estDate %in% speciesComp$estimatedDailyProportions$Date)) stop("Must run 'fitModel' before trying to predict species composition for this date.")
  estimatedP <- speciesComp$estimatedDailyProportions
  estimatedP$SonarCode <- ifelse(estimatedP$SonarBank == "Left Bank", "A1", "A2")
  
  estimatedP$mergecode <- paste(estimatedP$SonarCode, estimatedP$SonarAim, estimatedP$SonarBin)
  salmon$mergecode <- paste(salmon$SonarCode, salmon$SonarAim, salmon$SonarBin)
  keep <- 
  augmentdf <- salmon |> subset(!mergecode %in% estimatedP$mergecode & SonarAim != "Aim1+Aim2" & SonarBank != "Offshore")
  augmentdf <- augmentdf |> aggregate( count ~ Date + SonarBank + SonarCode + SonarAim, FUN = sum)
  if(any("Bin3" %in% estimatedP$SonarBin)){
    augmentdf$SonarBin <- "Bin3"
  }else{
    augmentdf$SonarBin <- "Bin2"
  }
  merge(augmentdf, augmentdf, by = c("Date", "SonarCode", "SonarBin", "SonarAim"), all.x = TRUE, suffixes = c("", "tmp"))
  
  ## Is Bin3 there? use it ow use Bin2.
  ## Duplicate estimatedP for Extra Bins:
  # augmentdf <- salmon |> 
  # if(!any(estimatedp$SonarBin == "Bin3")) {
  # for( i in c("Bin2")
  # estimatedP

  pTF <- speciesComp$testFisheryCounts |> subset(as.Date(Date) %in% estimatedP$Date)
  
  salmon <- speciesComp$salmonCounts |> subset(Date %in% estimatedP$Date)
  salmonspp <- grep("resident", speciesComp$species, value = TRUE, invert = TRUE)
  ## Mobile Prediction:
  pTF <- pTF |> merge(salmon |> subset(SonarBank == "Offshore"), by = "Date")
  
  for(i in salmonspp){

    if(!is.null(proportionOffshore)){
      pTF
    }
  }
}

plotmix <- function(speciesComp, day = 1, includeProportionLabels = FALSE, SonarBin = NULL, SonarAim = NULL, SonarBank = NULL, ...){
  if(!any(speciesComp$estDate %in% speciesComp$estimatedDailyProportions$Date)) 
    stop("Must run 'fitModel' before trying to predict species composition for this date.")

  mu <- speciesComp$estimatedParameters$mu
  sigma <- exp(speciesComp$estimatedParameters$logsigma)
  muChin <- speciesComp$estimatedParameters$muChin
  sigmaChin <- exp(speciesComp$estimatedParameters$logsigmaChin)
  sigma0 <- exp(speciesComp$estimatedParameters$logsigma0)
  beta <- speciesComp$estimatedParameters$beta

  dd <- speciesComp$estDate - day + 1
  stratum <- speciesComp$analysisData$predDF |> subset(Date == dd, select = stratum)
  stratum <- unique(stratum$stratum)

  if(!is.null(SonarBin)) stratum <- grep(SonarBin, stratum, value = TRUE)
  if(!is.null(SonarAim)) stratum <- grep(SonarAim, stratum, value = TRUE)
  if(!is.null(SonarBank)) stratum <- grep(SonarBank, stratum, value = TRUE)

  nst <- length(stratum)
  par(mfrow = c(ceiling(nst/2), 2))
  for( j in seq_along(stratum) ){
    s <- stratum[j]
    keep <- speciesComp$analysisData$lengthData$Date == dd & speciesComp$analysisData$lengthData$stratum == s
    x <- speciesComp$analysisData$lengthData$L.cm.adj[keep] - speciesComp$analysisData$Xlength[keep,] %*% beta
    wgts <- speciesComp$analysisData$lengthData$weights[keep]
    p <- speciesComp$estimatedDailyProportions |> subset(paste(SonarBank, SonarBin, SonarAim, sep = "_") == s & Date == speciesComp$estDate - d + 1)
    if(nrow(p) > 1) stop("I've gotten confused and found more than one proportion that matches with the data.")
    
    range <- histplot(x, wgts, xlab = "Fish Length (cm)", range = c(10, 110), main = paste(dd, gsub("_", " ", s)))
    col <- c("blue", "green")
    xx <- seq(range[1], range[2], by = 0.1)
    f <- numeric(length(xx))
    nonChinook <- grep("chinook", speciesComp$species, invert = TRUE, value = TRUE)
    for( i in seq_along(nonChinook)){
      fi <- p[1,nonChinook[i]]*dnorm(xx, mu[i], sqrt(sigma[i]^2 + sigma0^2))#/(1-pnorm(truncated, mu[i], sqrt(sigma[i]^2 + sigma0^2)))
      f <- f + fi
      lines(xx - range[1], fi, col = speciesColours(nonChinook[i]), lty = 1, lwd = 2) ## Have to adjust to scale of bar plot
      if(includeProportionLabels){
        mtext(paste0(speciesLabels(nonChinook[i]),  ": ", "prop: ", 
              round(p[1,nonChinook[i]], 3), ", mu: ", round(mu[i],1), ", sigma: ", 
              round(sigma[i],1)), side = 3, adj = 1, line = i-3, col = speciesColours(nonChinook[i]),
              cex = 0.75)
      }
    }
    ## Now add Chinook:
    chinook <- grep("chinook", speciesComp$species, invert = FALSE, value = TRUE)
    if(length(chinook) > 0){
    for(i in 1:2 ){
        fi <- p[1,chinook[i]]*dnorm(xx, muChin[i], sqrt(sigmaChin[i]^2 + sigma0^2))
        if(chinook[i] == "smalladultchinook"){
          fi <- fi + p[1,"largeadultchinook"]*dnorm(xx, muChin[i+1], sqrt(sigmaChin[i+1]^2 + sigma0^2))
        }
        f <- f + fi
        lines(xx - range[1], fi, col = speciesColours(chinook[i]), lty = 1, lwd = 2) ## Have to adjust to scale of bar plot
        if(includeProportionLabels){
          mtext(paste0(speciesLabels(chinook[i]),  ": ", "prop: ", 
                round(p[1,chinook[i]], 3), ", mu: ", round(muChin[i],1), ", sigma: ", 
                round(sigmaChin[i],1)), side = 3, adj = 1, line = length(nonChinook)+i-3, col = speciesColours(chinook[i]),
                cex = 0.75)
          if(chinook[i] == "smalladultchinook"){
            mtext(paste0(speciesLabels(chinook[i+1]),  ": ", "prop: ", 
                  round(p[1,chinook[i+1]], 3), ", mu: ", round(muChin[i+1],1), ", sigma: ", 
                  round(sigmaChin[i+1],1)), side = 3, adj = 1, line = length(nonChinook)+i+1-3, col = speciesColours(chinook[i+1]),
                  cex = 0.75)
          }
        }
      }
    }
    lines(xx-range[1], f, col = "black", lwd = 2, lty = 2)
    if(j==1){
      legend("topright", legend = speciesLabels(speciesComp$species), col = speciesColours(speciesComp$species), lty = rep(1, length(speciesComp$species)))
    }
  }
}

calculateResiduals <- function(speciesComp){
  ## Length Residuals:
  pars <- speciesComp$estimatedParameters
  mu <- pars$mu
  muChin <- pars$muChin
  ## Parameter Processing
  Kchin <- length(muChin)
  K0 <- length(mu)
  K <- Kchin + K0
  muc <- c(mu, muChin)
  alpha <- pars$alpha
  beta <- pars$beta
  alphaJackChinook <- pars$alphaJackChinook
  
  ## Standard deviation incld. observation error
  logsigma_ <- c(pars$logsigma, pars$logsigmaChin)
  sigma <- sqrt(exp(2*logsigma_) + exp(2*pars$logsigma0))

  ## Set up proportions. alpha parameter for predicting proportions. 
  ## Xalpha is a list of design matrices for each species.
  np <- nrow(speciesComp$analysisData$predDF)
  logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
  indx0 <- 1
  for( i in 1:K0 ){
    nc <- ncol(speciesComp$analysisData$Xprop[[i]]) 
    indx1 <- indx0 + nc - 1
    logitp[,i] <- as.matrix(speciesComp$analysisData$Xprop[[i]]) %*% alpha[indx0:indx1]
    indx0 <- indx1 + 1
  }

  logitpjack <- as.matrix(speciesComp$analysisData$XpropChin) %*% alphaJackChinook
  pjack <- 1/(1+exp(-logitpjack))

  p <- matrix(0, nrow = np, ncol = K)
  for( i in 1:np ) {
    p[i, 1:(K0+1)] <- expitM(logitp[i,])  
    p[i, (K0+1):K ] <- p[i, (K0+1)]*c(pjack[i], (1-pjack[i])*pars$pAdultChinook)
  }

  if(speciesComp$adjustLengths){
    L <- speciesComp$analysisData$lengthData$L.cm.adj - as.matrix(speciesComp$analysisData$Xlength) %*% beta
  }else{
    L <- speciesComp$analysisData$lengthData$L.cm.adj
  }
  wgts <- speciesComp$analysisData$lengthData$weights
  pobs <- p[speciesComp$analysisData$lengthData$grpIndex,]
  
  id <- sapply(1:length(L), FUN = function(i){which.max(pobs[i,]*dnorm(L[i], muc, sigma))}) ## Definitely not this one.
  residuals <- sapply(1:length(L), FUN = function(i){(L[i]-muc[id[i]])/sigma[id[i]]}) ## Definitely not this one.
    
  
  speciesComp$analysisData$lengthData$residuals <- residuals
  speciesComp$analysisData$lengthData$speciesid <- speciesComp$Species[id]
    
  boxplot(residuals ~ speciesid, data = speciesComp$analysisData$lengthData, xlab = "Species", ylab = "Residuals")
  abline(h = 0, col = 'red')
  qqnorm(residuals)
  qqline(residuals)
  plot(x = 1:length(L), y = residuals, xlab = "Index", ylab = "Residuals", cex = speciesComp$analysisData$lengthData$weights/(median(speciesComp$analysisData$lengthData$weights)))
  abline(h = 0, col = 'red', lty = 1)
  abline(h = c(-1.96, 1.96), col = 'red', lty = 2)

  ## Can print what the 'outlier' values are.
  speciesComp$analysisData$lengthData |> subset(abs(residuals) > 1.96)

}