## Set up an R6 Environment to run the mixture model and keep all the pieces for later.

## Delete me later:
##------------------------------------------------------------------------------------
library(tidyverse)
library(readxl)
sonarLengths <- read.csv("../example_data/2023MergedLengths.csv")
sonarCounts <- read.csv("../example_data/SONARMERGE export_2023.csv")
testFisheryCounts <- read_excel("../example_data/2023_WhonnockSpeciesComposition.xlsx")
estDate <- as.Date("2023-08-15")
##------------------------------------------------------------------------------------

library(R6)

runSpeciesComp <- R6Class("SpeciesComp",
  public = list(
    site = NULL,
    species = NULL, 
    Species = NULL,
    estDate = NULL,
    ndays = NULL,
    feasibleLengths <- c(10, 120)
    includeTestFishery <- FALSE,    
    rollAngle = 0,  ## *** Make sure to set this
    data = list(),
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
    processData = function(sonarCounts, sonarLengths, testFisheryCounts = NULL){
      if(!is.null(testFisheryCounts)){
        self$includeTestFishery <- TRUE
      }
      if(site == "Qualark") self$data <- processQualark(sonarcounts, sonarlengths, testfisherycounts, self$estDate, self$species, self$feasibleLengths)
      if(site == "Mission") self$data <- processMission(sonarCounts, sonarLengths, testFisheryCounts, self$estDate, self$species, self$feasibleLengths)
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
  match <- c("smallresident" = "Small Resident", "largeresident" = "Large Resident", "jackchinook" = "Jack Chinook", "pink" = "Pink", "sockeye" = "Sockeye", "chinook" = "Adult Chinook", "smalladultchinook" = "Small Adult Chinook", "largeadultchinook" = "Large Adult Chinook")
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
processMission <- function(sonarCounts, sonarLengths, testFisheryCounts = NULL, estDate, sppNames = NULL, feasibleLengths){
  year <- as.numeric(format(estDate, "%Y"))
  
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
  sonarLengths <- sonarLengths |> subset(SonarBin %in% c("Bin1", "Bin2", "Bin3") & SonarCode %in% c("A1", "A2") & L.cm >= feasibleLengths[1] & L.cm <= feasibleLengths[2])
  ## Create a look up code for fast merging:
  sonarLengths <- sonarLengths |> within(join_id <- factor(paste(SonarBin, SonarAim, SonarCode, Hour, MissionDate, sep = ""))) |>
                                  within(lookup_code <- as.numeric(join_id))
  sonarLengthsN <- sonarLengths |> aggregate( L.cm ~lookup_code, FUN = length)
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
  sonarLengths <- sonarLengths |> within(stratumWeight <- (Count/MinsCounted*15)/nLengths)  ## Assuming 15 min file.
  
  ## testFisheryCounts
  sonarLengths
  
}

obj <- runSpeciesComp$new(species = c("sockeye", "chinook"), site = "Mission")
obj$rollAngle
obj$estDate
obj$processData(sonarcounts, sonarlengths, testfishcounts, testfishlengths)
obj$optimControl

lengths, testcounts = NULL, testwgts = 1, fixedValues = list(), initialValues = list(), 
  K = 3, tol = 1e-8, Nmax = 1000, rel = TRUE, stratum = NULL, verbose = FALSE, qualark = FALSE, subsetTF = TRUE