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
##------------------------------------------------------------------------------------


library(R6)

speciesCompSummary <- R6Class("SpeciesComp",
  public = list(
    site = NULL,
    species = NULL, 
    Species = NULL,
    estDate = NULL,
    ndays = NULL,
    feasibleLengths = c(10, 120),
    includeTestFishery = FALSE,    
    rollAngle = 0,  ## *** Make sure to set this
    sonarLengths = NULL,
    testFisheryCounts = NULL,
    analysisData = NULL,
    parameters = list(),
    optimControl = list(tol = 1e-8, maxiters = 1000, relDiff = FALSE, verbose = FALSE),    
    
    initialize = function(species = c("jackChinook", "sockeye", "chinook"), estDate = NULL, ndays = 1, site = NULL) {
      if(missing(site)) stop("Need to initialize with a site = 'Mission' or 'Qualark'")
      if(!site %in% c("Mission", "Qualark")) stop("Need to initialize with a site = 'Mission' or 'Qualark'")
      self$species <- speciesCheck(species)
      self$Species <- speciesLabels(self$species)
      self$estDate <- checkDate(estDate)
      self$ndays <- ndays
      self$site <- site
      if(site == "Mission") self$rollAngle <- 0
      else self$rollAngle <- 35
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
    controlOptimization = function(tol = NULL, maxiters = NULL, relDiff = NULL, verbose = NULL){
      if(!is.null(tol)) self$optimControl$tol <- tol
      if(!is.null(maxiters)) self$optimControl$maxiters <- maxiters
      if(!is.null(relDiff)) self$optimControl$relDiff <- relDiff
      if(!is.null(verbose)) self$optimControl$verbose <- verbose
    },
    processData = function(sonarCounts, sonarLengths, testFisheryCounts = NULL, includeTangled = TRUE){
      if(!is.null(testFisheryCounts)){
        self$includeTestFishery <- TRUE
      }
      if(self$site == "Qualark") {
        processQualarkLengths(self, sonarcounts, sonarlengths, testfisherycounts)
      }
      if(self$site == "Mission") {
        processMissionLengths(self, sonarCounts, sonarLengths, testFisheryCounts)
        processWhonnockCounts(self, testFisheryCounts, tangled = includeTangled)
      }
    },
    setProportionsModel = function(){  ## Whatever values are in sonarLengths can be included here.
      
    }
  )
)

speciesCheck <- function(spp){
  allowed <- c("smallresident", "largeresident", "jackchinook", "pink", "sockeye", "chinook", "smalladultchinook", "largeadultchinook")
  spp <- tolower(gsub(" ", "", spp))
  if(any(!spp %in% allowed)) stop(paste("Species currently allowed are", allowed))
  if(any(spp == "chinook") & any(spp %in% c("smalladultchinook", "largeadultchinook"))){
    spp <- spp[-which(spp == "chinook")]
  }
  return(spp)
}

speciesLabels <- function(species){
  match <- c("smallresident" = "Small Resident", "largeresident" = "Large Resident", "jackchinook" = "Chinook Jack", 
             "pink" = "Pink", "sockeye" = "Sockeye", "chinook" = "Chinook Adult", "smalladultchinook" = "Chinook Small Adult", "largeadultchinook" = "Chinook Large Adult")
  as.character(match[species])
}

speciesColours <- function(species){
  c("smallresident" = "#FFB100", "largeresident" = "#656837",  "jackchinook" = "#A33CC7", "pink" = "#FF8DA1", "sockeye" = "#CD0000", "chinook" = "#27408B")
}

speciesOrder <- function(){
  c("smallresident", "largeresident", "jackchinook", "pink", "sockeye", "chinook")
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
processMissionLengths <- function(speciesComp, sonarCounts, sonarLengths, testFisheryCounts = NULL){
  
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
  MissionDateCheck <- sonarLengths |> subset(Date == MissionDate & Hour %in% 0:4)
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
  MissionDateCheck <- sonarCounts |> subset(Date == MissionDate & Hour %in% 0:4)
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

  sonarLengths$Date <- sonarLengths$MissionDate
  
  ## Add Stratum:
  sonarLengths <- sonarLengths |> within(stratum <- paste(SonarCode, SonarBin, SonarAim, sep = "_"))

  ## testFisheryCounts
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

  tmp <- data.frame("jackChinook" = testFisheryCounts$`Chinook Jack Gilled` + tangled*testFisheryCounts$`Chinook Jack Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("pink" = testFisheryCounts$`Pink All Gilled` + tangled*testFisheryCounts$`Pink All Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("sockeye" = testFisheryCounts$`Sockeye Adult Gilled` + tangled*testFisheryCounts$`Sockeye Adult Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("chinook" = testFisheryCounts$`Chinook Adult Gilled` + tangled*testFisheryCounts$`Chinook Adult Tangled`)
  out <- cbind(out, tmp)

  speciesComp$testFisheryCounts <- out
}

processQualarkLengths <- function(...){}

modelDataSetup <- function(speciesComp, formula = ~ stratum:day, chinookFormula = ~ 1, beamSpreadingFormula = ~ R.m){
  speciesComp$analysisData <- list()

  ndays <- speciesComp$ndays
  dates_ <- seq(speciesComp$estDate - ndays + 1, speciesComp$estDate, 1)

  ## Model Data:
  lengthData <- speciesComp$sonarLengths |> subset(Date %in% dates_) |> within(day <- as.numeric(factor(Date)))
  speciesComp$analysisData$Xp <- model.matrix(formula, data = lengthData)

  ## Connect Props to Test Fishery:
  speciesComp$analysisData$predictPropTF <- lengthData |> aggregate(update(formula, SalmonCount ~ .), FUN = sum)
  speciesComp$analysisData$Xtf <- model.matrix(formula, data = speciesComp$analysisData$predictPropTF)

  ## Predict chinook props: will assume they don't change per day for now...
  speciesComp$analysisData$Xc <- model.matrix(chinookFormula, data = lengthData)
  speciesComp$analysisData$predictPropTFc <- lengthData |> aggregate(update(chinookFormula, SalmonCount ~ .), FUN = sum)
  speciesComp$analysisData$Xtfc <- model.matrix(chinookFormula, data = speciesComp$analysisData$predictPropTFc)

  ## Beam Width Adjustment - 
  speciesComp$analysisData$Xl <- model.matrix(beamSpreadingFormula, data = lengthData)
}

modelParameterSetup <- function(speciesComp, formula = ~ stratum:day, chinookFormula = ~ 1, beamSpreadingFormula = ~ R.m){
  speciesComp$analysisData <- list()

  ndays <- speciesComp$ndays
  dates_ <- seq(speciesComp$estDate - ndays + 1, speciesComp$estDate, 1)

  ## Model Data:
  lengthData <- speciesComp$sonarLengths |> subset(Date %in% dates_) |> within(day <- as.numeric(factor(Date)))
  speciesComp$analysisData$Xp <- model.matrix(formula, data = lengthData)

  ## Connect Props to Test Fishery:
  speciesComp$analysisData$predictPropTF <- lengthData |> aggregate(update(formula, SalmonCount ~ .), FUN = sum)
  speciesComp$analysisData$Xtf <- model.matrix(formula, data = speciesComp$analysisData$predictPropTF)

  ## Predict chinook props: will assume they don't change per day for now...
  speciesComp$analysisData$Xc <- model.matrix(chinookFormula, data = lengthData)
  speciesComp$analysisData$predictPropTFc <- lengthData |> aggregate(update(chinookFormula, SalmonCount ~ .), FUN = sum)
  speciesComp$analysisData$Xtfc <- model.matrix(chinookFormula, data = speciesComp$analysisData$predictPropTFc)

  ## Beam Width Adjustment - 
  speciesComp$analysisData$Xl <- model.matrix(beamSpreadingFormula, data = lengthData)
}


  n <- nrow(lengths)
  lengths <- lengths %>% filter(!is.na(wgts))
  if(nrow(lengths) != n & verbose) cat("Warning: ", n-nrow(lengths), "observations have been removed due to NA `wgts` values.\n")

  if(any(names(lengths) == "MissionDate")){
    lengths <- lengths %>% select(-Date) %>% rename(Date = MissionDate)
  }
  if(any(names(lengths) == "L.cm.adj")){
    print("Using adjusted lengths `L.cm.adj`.")
    lengths <- lengths %>% select(-L.cm) %>% rename(L.cm = L.cm.adj)
  }

  ## Set up prob strata:
  stratum <- c("day", stratum)
  lengths_df <- lengths %>% 
    mutate(day = as.integer(as.factor(Date)), pcount = Up + Down) %>%
    select(L = L.cm, R = R.m, Date, wgts, pcount, {{stratum}}) %>%
    group_by(across(stratum)) %>%
    mutate(grp = cur_group_id()) %>%
    ungroup() %>%
    mutate(pcount_adj = pcount)

  if(length(stratum) > 1 & subsetTF == TRUE){
    if(qualark){
      lengths_df <- lengths_df %>% mutate(pcount = ifelse(Bank == "Right Bank", pcount, 0))
    }else{
      lengths_df <- lengths_df %>% mutate(pcount = ifelse(SonarBin != "Bin1", pcount, 0))    
    }
  }
  tfgrp <- lengths_df %>% 
            group_by(across(c("grp", stratum))) %>% 
            summarize(pcount = sum(pcount), pcount_adj = sum(pcount_adj), .groups = "drop")
  ndays <- max(lengths_df$day)
  ngrps <- max(lengths_df$grp)

  tfcount_mat <- NULL
  if(includeTF)
    tfcount_mat <- as.matrix(testcounts, nrow = ndays)

  if(length(testwgts) != ndays)
    testwgts <- rep(testwgts[1], ndays)
 
  ## If p isn't provided as an initial value come up with a basic initial value here:
  if(!any(names(initialValues) %in% "p")){
    allpars <- c(fixedValues, initialValues)
    breaks <- c(0, allpars$mu[1:K] + 1.5*sigma[1:K])
    breaks[length(breaks)] <- Inf
    lengths_df <- lengths_df %>% mutate(size_breaks = cut(L, breaks = breaks)) 
    pprior <- lengths_df %>% 
      group_by(grp) %>% mutate(n = n()) %>% 
      group_by(grp, size_breaks) %>% 
      summarize(pprior = n[1]/n(), .groups = "drop") %>%
      pivot_wider(values_from = pprior, names_from = size_breaks)
    pprior[is.na(pprior)] <- 0.00001
    pprior[pprior < 0.00001] <- 0.00001
    pprior <- as.matrix(pprior[,-1])
    pprior <- t(apply(pprior, 1, FUN = function(x){x/sum(x)}))
    initialValues$p <- pprior
  }
  
  input <- processParams(fixedValues, initialValues, K, ngrps)
  pars <- input$pars
  maplist <- input$map

  dataList <- list(L = lengths_df$L, R = lengths_df$R, day = lengths_df$day, grp = lengths_df$grp,
    wgts = lengths_df$wgts, testwgts = testwgts, testcounts = tfcount_mat, tfgrp = tfgrp)


#######################################################

speciesComp <- speciesCompSummary$new(species = c("sockeye", "chinook"), site = "Mission", estDate = "2023-08-05")
speciesComp$rollAngle
speciesComp$estDate
speciesComp$processData(sonarCounts, sonarLengths, testFisheryCounts)
sonarLengths <- speciesComp$sonarLengths
speciesComp$testFisheryCounts


speciesComp$optimControl

lengths, testcounts = NULL, testwgts = 1, fixedValues = list(), initialValues = list(), 
  K = 3, tol = 1e-8, Nmax = 1000, rel = TRUE, stratum = NULL, verbose = FALSE, qualark = FALSE, subsetTF = TRUE
  

ppgilled <- testFisheryCounts$`Pink All Gilled`/(testFisheryCounts$`Pink All Gilled` + testFisheryCounts$`Pink All Tangled`)
psgilled <- testFisheryCounts$`Sockeye Adult Gilled`/(testFisheryCounts$`Sockeye Adult Gilled` + testFisheryCounts$`Sockeye Adult Tangled`)
pcgilled <- testFisheryCounts$`Chinook Adult Gilled`/(testFisheryCounts$`Chinook Adult Gilled` + testFisheryCounts$`Chinook Adult Tangled`)

plot(pcgilled)
plot(testFisheryCounts$TRIP_DTT, ppgilled)
plot(testFisheryCounts$TRIP_DTT, psgilled)
plot(testFisheryCounts$TRIP_DTT, pcgilled)

plot(testFisheryCounts$TRIP_DTT, testFisheryCounts$`Pink All Tangled`)
points(testFisheryCounts$TRIP_DTT, testFisheryCounts$`Pink All Gilled`, col = 'red')

plot(testFisheryCounts$TRIP_DTT, testFisheryCounts$`Sockeye Adult Tangled`)
points(testFisheryCounts$TRIP_DTT, testFisheryCounts$`Sockeye Adult Gilled`, col = 'red')

plot(testFisheryCounts$TRIP_DTT, testFisheryCounts$`Chinook Adult Tangled`)
points(testFisheryCounts$TRIP_DTT, testFisheryCounts$`Chinook Adult Gilled`, col = 'red')
