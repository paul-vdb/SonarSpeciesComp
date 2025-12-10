## Set up an R6 Environment to run the mixture model and keep all the pieces for later.

library(R6)

runSpeciesComp <- R6Class("SpeciesComp",
  public = list(
    species = NULL, 
    Species = NULL,
    estDate = NULL,
    ndays = NULL,
    # rollAngle = 0,  ## *** Make sure to set this
    data = list(),
    parameters = list(),
    optimControl = list(tol = 1e-8, maxiters = 1000, relDiff = FALSE, verbose = FALSE),    
    initialize = function(species = c("jackChinook", "sockeye", "chinook"), estDate = NULL, ndays = 1) {
      self$species <- speciesCheck(species)
      self$Species <- speciesLabels(self$species)
      self$estDate <- checkDate(estDate)
      self$ndays <- ndays
    },
    ## Set the species names: Convenient and for printing.
    setSpecies = function(species = NULL){
      if(!is.null(species)){
        self$species <- speciesCheck(species)
        self$Species <- speciesPrint(species)
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
    processData = function(sonarcounts, sonarlengths, testfisherylengths, testfisherycounts){
      data <- processData(sonarcounts, sonarlengths, testfisherylengths, testfisherycounts)
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

## Data processing for counts + lengths to add wgts to lengths etc.
processData(sonarcounts, sonarlengths, testfisherylengths = NULL, testfisherycounts = NULL){
  
}

obj <- runSpeciesComp$new()
obj$optimControl

lengths, testcounts = NULL, testwgts = 1, fixedValues = list(), initialValues = list(), 
  K = 3, tol = 1e-8, Nmax = 1000, rel = TRUE, stratum = NULL, verbose = FALSE, qualark = FALSE, subsetTF = TRUE