#' Species Composition Estimation Methods
#'
#' Object to step through all the different setups and steps of estimating species composition on the Fraser River (Qualark and Mission) using Hydroacoustics and the test fishery.
#'
#' @param species vector of species names  e.g. c("largeresident", "jackChinook", "sockeye", "adultchinook").
#' @param date Date that you will estimate. Place holder but does not need to be set.
#' @param ndays Number of days to combine, defaults to 1.
#' @param site Choose between "Mission" or "Qualark", roll angle is set be default based on site.
#' @param roll_angle The roll angle that the sonar is set to (default = NULL).
#' @param feasible_lengths Lower and upper bound of feasible lengths that can be measured on the sonar to remove infeasible (default = c(10, 120)).
#'
#' @details This is an R6 object that holds all of the methods to estimate species composition and generate reports. The methods include:
#' \itemize{
#'  \item{setSpecies}{Update the object with species to fit in a model. Input: 'species': options are c("smallresident", "largeresident", "pink", "sockeye", "jackchinook", "adultchinook", "smalladultchinook", "largeadultchinook").}
#'  \item{setDate}{Update the estimation date and number of days to use. Input: 'date' and 'ndays'}
#'  \item{controlOptimization}{Control optimization of the EM algorithm. Input: list with values `$tolerance` what difference is okay (default 1e-8), `$maxiters` maximum iterations of EM algorithm (default 1000), 
#'  `$relative_difference` true or false to choose relative vs absolute difference on each iteration of algorithm, `$verbose` true or false to print details of EM algorithm.}
#'  \item{processData}{Core function for adding data to the object. Does the main data processing and adds a few columns needed elsewhere for the entire dataset added to make it fast and easy to do different days. 
#'  Input: 'sonar_counts' and 'sonar_lengths' are data frames see \code{?process_mission_lengths} or \code{?process_qualark_lengths} for details, 
#'  'test_fishery_counts' list of dataframes named after the test fishery (e.g. 'whonnock' or 'albion') see \code{?process_whonnock_catch} or \code{?process_albion_catch} for details, 
#'  'include_tangled' Whether or not to include tangled or just gilled counts, 'dropN' Keep length measurements per hour/bin/bank etc. if > dropN (default 2) but not always necessary, 
#'  'salmon_passage table' data frame see \code{?process_mission_salmon_passage} for details, and 'feasible_lengths' set by default but can be changed here for upper and lower bounds of reasonable sonar lengths.}
#'  \item{setSpeciesLengths}{Input: 'mu' = NULL, 'sigma' = NULL, 'proportion_adultchinook' = NULL, 'test_fishery_lengths' = NULL, 'ndays' = 6. See \code{?set_species_lengths} for details.}
#'  \item{setModelParameters}{Core function to set the model before fitting with EM algorithm, especially which values to fix and which to fit. Input: 
#'  'fixed_parameters' - which parameters to hold fixed (default c("mu", "sigma", "proportion_jackchinook")),
#'  'fixed_values' - list fixed values that correspond to fixed_parameters, possible to set just a single value within as fixed (e.g. \code{list(mu = c("jackchinook" = 40))}),
#'  'initial_values' - list provided by the user to improve fitting time for estimated parameters, default are used if not provided.
#'  'formula_proportions' - list of formulas to use for alpha for each species, if \code{list("all" = ~ factor(day))}, then each species is given the same formula, except jack chinook are treated separately.
#'  'formula_lengths' - formula for how to deal with beam spreading (default = \code{~beamWidth.cm}, 'date' - Option to change here if wanted, 'ndays' - Option to change number of days used to estimate, 
#'  'species' - option to choose new species, 'delta_mu_bounds' - choose how much the mean length of each species can vary, 'qinv_names' - how to allow catchability to vary with combinations \code{c("species", "fishery", "net_type")},
#'  'test_fishery_spp' -  combine test fishery and species to inform which species within each test fishery to use e.g. ("sockeye_whonnock", "adultchinook_whonnock", "chinook_albion").}
#'  \item{fitModel}{Runs EM algorithm based on object values for input. Input: 'simulatedDate' - logical if wanting to fit to simulatd data (default = FALSE), 
#'  'control' - list containing controlOptimization arguments as well as '$include_test_fishery' - logical to include or exclude test fishery catch information in the model,
#'  '$adjust_lengths' - logical whether or not to adjust measured lengths for beam spreading, and '$test_fishery_weights' - values to apply for weighting the test fishery. Defaults are 
#'  whatever values are currently stored in the model object, which may be from the previous run. Results are saved as part of object: \code{$params_estimated} and also returns 
#'  the list of estimates.}
#'  \item{simulate}{Simulate the model for the counts on the 'est_date' and 'ndays' with parameter values equal to 'params_estimated' within the object. See \code{?simulateGivenData} for details.}
#'  \item{plotMix}{Plot the fitted mixture model as a histogram. Input: 'day' (default = 1), 'include_proportion_labels' - logical (default = FALSE), '...' - additional plot arguments.}
#'  \item{plotTestFishery}{Plot Pearson residuals for the test fishery data (CPUE - Nq) against each species, test fishery, and net type.}
#'  \item{plotBeamSpreading}{Plot estimated beam spreading effects against sonar range.}
#' }
#'
#' @returns An R6 object with methods and values used for fitting species composition models.
#'
#' @examples
#' self <- speciesCompModel$new(species = c("largeresident", "jackchinook", "sockeye", "adultchinook"), site = "Mission", date = "2023-08-05")
#' 
#' @export
speciesCompModel <- R6Class("SpeciesComp",
  public = list(
    species_info = NULL,
    est_date = NULL,
    ndays = NULL,
    data_list = NULL,
    data_info = NULL,
    fit_info = NULL,
    model = NULL,
    sonar_lengths = NULL,
    test_fishery_catch = NULL,
    salmon_counts = NULL,
    params_estimated = list(),
    params_fixed = list(),
    params_init = list(),
    estimates = list(),
    sim_data = NULL,
    default_parameters = NULL,
    prior_distributions = NULL,
    prior_jacobians = NULL,    
    
    initialize = function(species = c("jackchinook", "sockeye", "adultchinook"), date = NULL, ndays = 1, site = NULL, roll_angle = NULL, feasible_lengths = c(10, 120)) {
      if( missing(site) ) stop("Need to initialize with a site = 'Mission' or 'Qualark'")
      if( !site %in% c("Mission", "Qualark") ) stop("Need to initialize with a site = 'Mission' or 'Qualark'")
      self$species_info <- speciesCheck(species)
      self$est_date <- checkDate(date)
      self$ndays <- ndays

      self$fit_info <- list()
      self$fit_info$include_test_fishery <- FALSE
      self$fit_info$adjust_lengths <- FALSE
      self$fit_info$tolerance <- 1e-8
      self$fit_info$maxiters <- 1000
      self$fit_info$relative_difference <- FALSE
      self$fit_info$verbose <- FALSE

      ## Data specific Information to hold
      self$data_info <- list()
      self$data_info$feasible_lengths <- feasible_lengths
      self$data_info$site <- site
      self$data_info$roll_angle <- roll_angle
      self$data_info$test_fishery_weights <- 1
      
      if(site == "Mission" & is.null(self$data_info$roll_angle)) self$data_info$roll_angle <- 0
      if(site == "Qualark" & is.null(self$data_info$roll_angle)) self$data_info$roll_angle <- 35
      
      if(site == "Mission"){
        self$data_info$test_fishery_input <- c("sockeye_whonnock_vmn", "chinook_whonnock_vmn", "chinook_albion_chin", "chinook_albion_vmn")
        self$data_info$test_fishery_spp <- c("sockeye_whonnock", "chinook_whonnock", "chinook_albion")
      }
      if(site == "Qualark"){
        self$data_info$test_fishery_input <- c("sockeye_qualark", "chinook_qualark")
        self$data_info$test_fishery_spp <- c("sockeye_qualark", "chinook_qualark")
      }
      self$data_info$test_fishery_formula <- c("species", "fishery", "net_type")
      
      self$test_fishery_catch <- list()
      self$default_parameters <- default_params()
      default_prior <- \(...){0}
      self$prior_distributions <- list(dlog_q = default_prior, dbeta = default_prior, dalpha_jackchinook = default_prior, 
                                       dlog_sigma = default_prior, dmu = default_prior, dlogit_delta_mu = default_prior, dlog_sigma0 = default_prior,
                                       dsmallresident = default_prior, dlargeresident = default_prior)
      self$prior_jacobians <- list(dlog_q = default_prior, dbeta = default_prior, dalpha_jackchinook = default_prior, 
                                       dlog_sigma = default_prior, dmu = default_prior, dlogit_delta_mu = default_prior, dlog_sigma0 = default_prior,
                                       dsmallresident = default_prior, dlargeresident = default_prior)                                       
    },
    ## Set the species names: Convenient and for printing:
    setSpecies = function(species = NULL){
      if(!is.null(species)){
        if(!all(species %in% self$species_info$species)){
          cat("Resetting default Parameters.\n")
          self$default_parameters <- default_params()
        }
        self$species_info <- speciesCheck(species)
      }
    },
    ## Set date:
    setDate = function(date = NULL, ndays = NULL){
      if(!is.null(date)) self$est_date <- checkDate(date)
      if(!is.null(ndays)) self$ndays <- ndays
    },
    controlOptimization = function(control = list()){
      self$fit_info$tolerance <-  extractControls(control$tolerance, 1e-8)
      self$fit_info$maxiters <-  extractControls(control$maxiters, 1000)
      self$fit_info$relative_difference<- extractControls(control$relative_difference, 1000)
      self$fit_info$verbose <- extractControls(control$verbose, FALSE)
    },
    processData = function(sonar_counts, sonar_lengths, test_fishery_counts = NULL, include_tangled = TRUE, dropN = 2, salmon_passage_table = NULL, feasible_lengths = NULL){
      if(!is.null(feasible_lengths)) self$data_info$feasible_lengths <- feasible_lengths
      if(!is.null(test_fishery_counts)) self$fit_info$include_test_fishery <- TRUE

      if(self$data_info$site == "Qualark") {
        process_qualark_data(self, net_length = 30)
      }
      if(self$data_info$site == "Mission") {
        process_mission_lengths(self, sonar_counts, sonar_lengths, dropN)
        if("whonnock" %in% names(test_fishery_counts)) process_whonnock_catch(self, test_fishery_counts[["whonnock"]], tangled = include_tangled)
        if("albion" %in% names(test_fishery_counts)) process_albion_catch(self, test_fishery_counts[["albion"]] )        
        if(!is.null(salmon_passage_table)) process_mission_salmon_passage(self, salmon_passage_table)
      }
    },
    setSpeciesLengths = function(mu = NULL, sigma = NULL, proportion_adultchinook = NULL, test_fishery_lengths = NULL, ndays = 6){
      set_species_lengths(self, mu, sigma, proportion_adultchinook, test_fishery_lengths, ndays)
    },
    setModelParameters = function(fixed_parameters = c("mu", "sigma", "proportion_jackchinook"),
                                  fixed_values = list(), initial_values = list(),
                                  formula_proportions = list(),
                                  formula_lengths = ~ beamWidth.cm,
                                  date = NULL, ndays = NULL, species = NULL, delta_mu_bounds = NULL,
                                  qinv_names = NULL,
                                  test_fishery_spp = NULL){
      ## If user wants to update the date and species, they can do that here:
      if(!is.null(qinv_names)) self$data_info$test_fishery_formula <- qinv_names
      if(!is.null(test_fishery_spp)) self$data_info$test_fishery_spp <- test_fishery_spp

      self$setDate(date, ndays)
      self$setSpecies(species)
      set_daily_data(self)

      set_length_adjustment(self, formula_lengths)
      set_model_proportions(self, formula_proportions)
      set_model_parameters(self, fixed_parameters, fixed_values, initial_values, delta_mu_bounds)      
    },
    fitModel = function(simulatedData = FALSE, control = list()){
      self$controlOptimization(control)
      self$fit_info$include_test_fishery <- extractControls(control$include_test_fishery, self$fit_info$include_test_fishery)
      self$data_info$test_fishery_weights <- extractControls(control$test_fishery_weights, self$data_info$test_fishery_weights)
      self$fit_info$adjust_lengths <- extractControls(control$adjust_lengths, self$fit_info$adjust_lengths)
    
      fit_joint_model(self)
    },
    simulate = function(){
      simulateGivenData(self)
    },
    setPriors = function(priors = list(), includeJacobian = TRUE){
      set_priors(self, priors = priors, includeJacobian)
    },
    plotMix = function(day = 1, include_proportion_labels = FALSE, ...){
      plot_mix(self, day = day, include_proportion_labels = include_proportion_labels, ...)
    },
    plotTestFishery = function(){
      plot_test_fishery(self)
    },
    plotBeamSpreading = function(){
      plot_beam_spreading(self)
    }
  )
)

default_params <- function(){
  default_parameters <- list()
  ## Expansion Line Albion: 3049 (2020 Brittany Jenewein)
  default_parameters$qinv <- c("pink_whonnock" = 35000, "sockeye_whonnock" = 12000, "jackchinook_whonnock" = 20000, "adultchinook_whonnock" = 3000,
                               "pink_whonnock_vmn" = 35000, "sockeye_whonnock_vmn" = 12000, "jackchinook_whonnock_vmn" = 20000, "adultchinook_whonnock_vmn" = 3000,   
                              "chinook_whonnock" = 3000, "chinook_whonnock_vmn" = 3000, 
                              "chinook_albion_vmn" = 3049, "chinook_albion_chin" = 6049, "chinook_albion_chum" = 3049,                              
                              "chinook" = 3049, "chinook_albion" = 3049, "chinook_vmn" = 3049, "chinook_chin" = 3049, "chinook_chum" = 3049,
                              "sockeye" = 12000, "sockeye_vmn" = 12000, "pink" = 35000, "pink_vmn" = 35000,
                              "pink_qualark" = 10000, "sockeye_qualark" = 5000, "jackchinook_qualark" = 20000,
                              "adultchinook_qualark" = 5000, "chinook_qualark" = 6000, "coho_albion" = 2548, "chum_albion" = 3201, 
                              "coho_albion_vmn" = 2548, "coho_albion_chin" = 2548, "coho_albion_chum" = 2548,
                              "chum_albion_vmn" = 3201, "chum_albion_chin" = 3201, "chum_albion_chum" = 3201)
  default_parameters$alpha <- c(0,0,0,0)
  default_parameters$alpha_jackchinook <- 0
  default_parameters$beta <- c(0,0)
  default_parameters$mu <- c("smallresident" = 23, "largeresident" = 35, "pink" = 50, "sockeye" = 58, "coho" = 65, "chum" = 70, "jackchinook" = 40, "smalladultchinook" = 62, "largeadultchinook" = 80, "adultchinook" = 70)
  default_parameters$sigma <- c("smallresident" = 2.5, "largeresident" = 2.5, "pink" = 4, "sockeye" = 3.5, "coho" = 3.5, "chum" = 3, "jackchinook" = 2.5, "smalladultchinook" = 4, "largeadultchinook" = 6, "adultchinook" = 7)
  default_parameters$sigma0 <- 4.5
  default_parameters$delta_mu <- c("smallresident" = 0, "largeresident" = 0, "pink" = 0, "sockeye" = 0, "coho" = 0, "chum" = 0, "jackchinook" = 0, "smalladultchinook" = 0, "largeadultchinook" = 0, "adultchinook" = 0)
  default_parameters$proportion_adultchinook <- c("smalladultchinook" = 0.2, "largeadultchinook" = 0.8, "adultchinook" = 1)
  default_parameters
}

#' @export
set_priors <- function(self, priors = list(), includeJacobian = TRUE){        
  default <- function(...){0}
 
  self$prior_distributions <- list()
  self$prior_jacobians <- list()
                                      
  self$prior_distributions$dlog_qinv <- extractControls(priors$qinv, default)
  self$prior_jacobians$dlog_qinv <- addJacobian(\(x){sum(log(abs(x)))}, includeJacobian & !is.null(priors$qinv))                                   
  self$prior_distributions$dlog_sigma <- extractControls(priors$sigma, default)
  self$prior_jacobians$dlog_sigma <- addJacobian(\(x){sum(log(abs(x)))}, includeJacobian & !is.null(priors$sigma))
  self$prior_distributions$dlog_sigma0 <- extractControls(priors$sigma0, default)
  self$prior_jacobians$dlog_sigma0 <- addJacobian(\(x){log(abs(x[1]))}, includeJacobian & !is.null(priors$sigma0))
  self$prior_distributions$dlogit_delta_mu <- extractControls(priors$delta_mu, default)
  self$prior_jacobians$dlogit_delta_mu <- addJacobian(logDetJac_logitInterval, includeJacobian & !is.null(priors$delta_mu))

  self$prior_distributions$dmu <- extractControls(priors$mu, default)
  self$prior_distributions$dbeta <- extractControls(priors$beta, default)
  self$prior_distributions$dalpha_jackchinook <- extractControls(priors$alpha_jackchinook, default)
  self$prior_distributions$dlargeresident <- extractControls(priors$N_largeresident, default)
  self$prior_distributions$dsmallresident <- extractControls(priors$N_smallresident, default)
}

#' @export
speciesCheck <- function(spp){
  allowed <- c("smallresident", "largeresident", "pink", "sockeye", "coho", "chum", "jackchinook", "adultchinook", "smalladultchinook", "largeadultchinook")
  spp <- tolower(gsub(" ", "", spp))
  if(any(!spp %in% allowed)) stop(paste("Species currently allowed are", allowed))
  if(any(spp == "adultchinook") & any(spp %in% c("smalladultchinook", "largeadultchinook"))){
    spp <- spp[-which(spp == "adultchinook")]
  }
  spp_info <- list()
  spp_info$species <- allowed[which(allowed %in% spp)]
  spp_info$nspecies <- length(spp_info$species)
  spp_info$species_other <- grep("chinook", spp_info$species, invert = TRUE, value = TRUE)
  spp_info$nother <- length(spp_info$species_other)
  spp_info$species_chinook <- grep("chinook", spp_info$species, invert = FALSE, value = TRUE)
  spp_info$nchinook <- length(spp_info$species_chinook)
  spp_info$labels <- speciesLabels(spp_info$species)
  spp_info$chinook_indx <- grep("chinook", spp_info$species)
  spp_info$adultchinook_indx <- grep("adultchinook", spp_info$species)
  spp_info$other_indx <- grep("chinook", spp_info$species, invert = TRUE)

  species_N <- spp_info$species[spp_info$species != "smallresident"]
  species_N <- grep("adultchinook", species_N, invert = TRUE, value = TRUE)
  if(any(grepl("adultchinook", spp_info$species))) species_N <- c(species_N, "adultchinook", "chinook")
  spp_info$species_predict <- species_N
  
  return(spp_info)
}

#' Check value is a in correct Date format
#'
#' @param e_date Value to check.
#'
#' @return A date in Y-m-d format.
#'
#' @export
checkDate <- function(e_date){
  if(is.null(e_date)){return(NULL)}
  if(!is.Date(e_date)){
    e_date <- as.Date(e_date, "%Y-%m-%d")
    if(is.na(e_date)) stop("Date must be submitted as 'Y-m-d', e.g. '2018-08-05'")
  }
  e_date
}

#' Process Mission Sonar Lengths
#'
#' @param self R6 speciesCompModel Object.
#' @param sonar_counts Sonar Lengths in standard Mission format.
#' @param sonar_lengths Sonar Lengths in standard Mission Format.
#' @param dropN Minimum sample size in a stratum.
#'
#' @return No object returned, appends sonar_lengths to self.
#'
#' @export
process_mission_lengths <- function(self, sonar_counts, sonar_lengths, dropN = 2){
  
  year <- as.numeric(format(self$est_date, "%Y"))
  
  ## Check names in dataframes are present:
  count_col_names <- c("Date", "MissionDate", "Hour", "MinsCounted", "SonarCode", "SonarAim", "SonarBin", 
                     "Up.0to5", "Dn.0to5", "Up.5to10", "Dn.5to10")
  length_col_names <- c("Date", "MissionDate", "Hour", "SonarCode", "SonarAim", "SonarBin", 
                 "R.m", "L.cm", "Dir", "EnterTheta.deg")
                 
  namecheck <- which(!count_col_names %in% names(sonar_counts))
  if(length(namecheck)) stop(paste0("We require column names in sonarCounts: ", paste(count_col_names[namecheck], collapse = ",")))
  namecheck <- which(!length_col_names %in% names(sonar_lengths))
  if(length(namecheck)) stop(paste0("We require column names in sonarLengths: ", paste(length_col_names[namecheck], collapse = ",")))

  ## Clean up length data:
  if(!is.Date(sonar_lengths$Date)) sonar_lengths$Date <- as.Date(sonar_lengths$Date, format = "%m/%d/%Y")
  if(!is.Date(sonar_lengths$MissionDate)) sonar_lengths$MissionDate <- as.Date(sonar_lengths$MissionDate, format = "%m/%d/%Y")
  sonar_lengths$Hour <-   as.numeric(sonar_lengths$Hour)
  ## Check Mission Date:
  mission_date_check <- sonar_lengths |> subset(Date == MissionDate & Hour < 5)
  if(nrow(mission_date_check) > 0) stop("Invalid Mission dates for sonar lengths detected where 'Date == MissionDate' at Hours 0-4.")

  ## Subset Count data to match length data:
  sonar_lengths <- sonar_lengths |> subset(as.numeric(format(Date, "%Y")) == year)
  sonar_lengths <- sonar_lengths |> subset(SonarBin %in% c("Bin1", "Bin2", "Bin3") & SonarCode %in% c("A1", "A2") & L.cm >= self$data_info$feasible_lengths[1] & L.cm <= self$data_info$feasible_lengths[2])
  ## Create a look up code for fast merging:
  sonar_lengths <- sonar_lengths |> within(join_id <- factor(paste(SonarBin, SonarAim, SonarCode, Hour, MissionDate, sep = ""))) |>
                                  within(lookup_code <- as.numeric(join_id))
  sonar_lengthsN <- sonar_lengths |> aggregate( L.cm ~ lookup_code, FUN = length)
  names(sonar_lengthsN)[names(sonar_lengthsN) == 'L.cm'] <- 'nLengths'
  sonar_lengths <- sonar_lengths |> merge(sonar_lengthsN, all.x = TRUE, all.y = TRUE)
  
  ## Subset Count data to match length data:
  sonar_counts <- sonar_counts |> subset(SonarBin %in% c("Bin1", "Bin2", "Bin3") & SonarCode %in% c("A1", "A2") & Up.0to5 != "" & MinsCounted > 0)
  suppressWarnings(sonar_counts <- sonar_counts |> within(Up <- ifelse(MinsCounted == 5, as.numeric(Up.0to5), as.numeric(Up.0to5) + as.numeric(Up.5to10))))
  suppressWarnings(sonar_counts <- sonar_counts |> within(Down <- ifelse(MinsCounted == 5, as.numeric(Dn.0to5), as.numeric(Dn.0to5) + as.numeric(Dn.5to10))))
  sonar_counts <- sonar_counts |> within(MinsCounted <- ifelse(is.na(Up.5to10), MinsCounted, 10))
  sonar_counts <- sonar_counts |> within(Count <- Up + Down) |> within(SalmonCount <- Up - Down) |> subset(!is.na(Count))
  
  ## check Date and Mission Date: Format should be m/d/Y.
  if(!is.Date(sonar_counts$Date)) sonar_counts$Date <- as.Date(sonar_counts$Date, format = "%m/%d/%Y")
  if(!is.Date(sonar_counts$MissionDate)) sonar_counts$MissionDate <- as.Date(sonar_counts$MissionDate, format = "%m/%d/%Y")

  ## Check Mission Date:
  sonar_counts$Hour <-   as.numeric(sonar_counts$Hour)
  mission_date_check <- sonar_counts |> subset(Date == MissionDate & Hour < 5)
  if(nrow(mission_date_check) > 0) stop("Invalid Mission dates for sonar counts detected where 'Date == MissionDate' at Hours 0-4.")

  suppressMessages(
    sonar_counts <- sonar_counts |> within(join_id <- factor(paste(SonarBin, SonarAim, SonarCode, Hour, MissionDate, sep = ""), levels = levels(sonar_lengths$join_id))) |>
                   subset(!is.na(join_id)) |> within(lookup_code <- as.numeric(join_id))
  )
  ## Merge count data to lengths to get correction to est proportion with test fishery count and stratum weights:
  sonar_lengths <- sonar_lengths |> merge(sonar_counts[, c("lookup_code", "Count", "SalmonCount", "MinsCounted")])
  sonar_lengths <- sonar_lengths |> within(weights <- (Count/MinsCounted*5)/nLengths)  ## Assuming 5 min file taken.
  
  minhour <- min(sonar_counts$Hour[sonar_counts$MissionDate == sonar_counts$Date], na.rm = TRUE)

  ## Account for "Mission Date" and then just use "Date".
  sonar_lengths$HourOrder <- sonar_lengths$Hour - minhour
  sonar_lengths$HourOrder[sonar_lengths$Hour < minhour] <- 24 + sonar_lengths$HourOrder[sonar_lengths$Hour < minhour]

  sonar_lengths$Date <- sonar_lengths$MissionDate ## Rename MissionDate to Date for convenience.
  
  ## Add Stratum:
  sonar_lengths <- sonar_lengths |> within(SonarBank <- ifelse(SonarCode == "A1", "Left Bank", "Right Bank"))  
  sonar_lengths <- sonar_lengths |> within(stratum <- paste(SonarBank, SonarBin, SonarAim, sep = "_"))
  ## Add beam width:
  sonar_lengths <- sonar_lengths |> within(beamWidth.cm <- R.m*0.3*pi/180*100)

  ## Will remove single lengths measured in an hour:
  sonar_lengths <- sonar_lengths |> filter(nLengths > dropN)

  ## I don't like SonarAim being called 3 on right bank if only one aim.
  sonar_lengths <- sonar_lengths |> within(SonarAimF <- ifelse(SonarBank == "Right Bank", "Aim3", SonarAim))

  self$sonar_lengths <- sonar_lengths
}


#' Process Whonnock Catch
#'
#' @param self R6 speciesCompModel Object.
#' @param test_fishery_counts Data frame in standard format for Whonnock catch.
#' @param tangle Logical whether to include or exclude tangled counts.
#'
#' @return No object returned, appends test fishery catch to \code{self$test_fishery_catch[["whonnock"]]}
#'
#' @export
process_whonnock_catch <- function(self, test_fishery_counts, tangled = TRUE){

  ## Process times:
  test_fishery_counts <- test_fishery_counts |> 
                         within(start_out <- process_minutes(start_out)) |>
                         within(full_out <- process_minutes(full_out)) |>
                         within(start_in <- process_minutes(start_in)) |>
                         within(full_in <- process_minutes(full_in))  

  ## Soak Time and Effort:
  ## *** pvdb Don't like correction factor. Discounting the impact of in-time makes the CPUE relationship linear...
  test_fishery_counts <- test_fishery_counts |> 
                         within(soak_time <- (full_out-start_out)/2 + (start_in - full_out) + (full_in-start_in)^0.5/2) |>
                         within(effort <- soak_time*net_length/1000)
  
  out <- test_fishery_counts[,c("TRIP_DTT", "effort", "net_length")]
  names(out)[1] <- "Date"
  out$Date <- as.Date(out$Date)
  
  tmp <- data.frame("pink" = test_fishery_counts$`Pink All Gilled` + tangled*test_fishery_counts$`Pink All Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("sockeye" = test_fishery_counts$`Sockeye Adult Gilled` + tangled*test_fishery_counts$`Sockeye Adult Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("jackchinook" = test_fishery_counts$`Chinook Jack Gilled` + tangled*test_fishery_counts$`Chinook Jack Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("adultchinook" = test_fishery_counts$`Chinook Adult Gilled` + tangled*test_fishery_counts$`Chinook Adult Tangled`)
  out <- cbind(out, tmp)

  out <- out |> within(chinook <- adultchinook + jackchinook)
  
  out <- out |> aggregate(cbind(pink, sockeye, jackchinook, adultchinook, chinook, effort) ~ Date, sum)

  out$net_type <- "vmn"

  self$test_fishery_catch[["whonnock"]] <- out
}

#' @export
process_mission_salmon_passage = function(self, salmon_passage_table){
  ## Check names in dataframes are present:
  col_names <- c("SelectedMethod", "LB_BIN1_AIM1", "LB_BIN1_AIM2", "LB_BIN1", "LB_BIN2", 
                     "LB_BIN3", "LB_BIN4", "LB_BIN5", "LB_BIN6", "Offshore", "RB_BIN4", "RB_BIN3",
                     "RB_BIN2", "RB_BIN1", "MissionDate")

  namecheck <- which(!grep("BIN4|BIN5|BIN6", col_names, invert = TRUE, value = TRUE) %in% names(salmon_passage_table))
  if(length(namecheck)) stop(paste0("We require column names in salmonPassageTable: ", paste(col_names[namecheck], collapse = ",")))

  namedf <- data.frame(grp = c("LB_BIN1_AIM1", "LB_BIN1_AIM2", "LB_BIN1", "LB_BIN2", "LB_BIN3", "LB_BIN4", "LB_BIN5", 
                               "LB_BIN6", "Offshore", "RB_BIN4", "RB_BIN3", "RB_BIN2", "RB_BIN1"),
                      SonarBank = c(rep("Left Bank",8), "Offshore", rep("Right Bank", 4)),
                      SonarCode = c(rep("A1",8), "mobile", rep("A2", 4)),
                      SonarBin = c("Bin1", "Bin1", "Bin1", "Bin2", "Bin3", "Bin4", "Bin5", "Bin6", NA, "Bin4", "Bin3", "Bin2", "Bin1"),
                      SonarAim = c("Aim1", "Aim2","Aim1+Aim2", "Aim3", "Aim3", "Aim3", "Aim3", "Aim3", NA, "Aim1", "Aim1", "Aim1", "Aim1")
                      )

  ## Clean up length data:
  if(!is.Date(salmon_passage_table$MissionDate)) salmon_passage_table$MissionDate <- as.Date(salmon_passage_table$MissionDate, format = "%m/%d/%Y")
  salmon_passage_table$Date <- salmon_passage_table$MissionDate

  long_salmon <- do.call("rbind", lapply(namedf$grp[namedf$grp %in% names(salmon_passage_table)], FUN = function(x){
    tab = cbind(salmon_passage_table[, c("SelectedMethod", "Date", "MissionDate", x)], "grp" = x); colnames(tab)[4] <-  "count"; tab
    }))
  long_salmon <- long_salmon |> subset(!is.na(count))
  long_salmon <- long_salmon |> merge(namedf)
  self$salmon_counts <- long_salmon
}

#' @export
set_daily_data <- function(self){

  self$data_list <- list()
  ## Model Lengths Data:
  ndays <- self$ndays
  dates_ <- seq(self$est_date - ndays + 1, self$est_date, 1)
  length_data <- self$sonar_lengths |> subset(Date %in% dates_) |> within(day <- as.numeric(factor(Date)))
  length_data <- length_data |> within(L.cm.adj <- L.cm/cos(self$data_info$roll_angle/180*pi))
  
  length_data <- length_data |> 
                    within(grpID <- factor(paste(day, Hour, stratum, sep = "_"))) |>
                    within(grpIndex <- as.numeric(grpID))
  
  self$data_list$length_data <- length_data

  ## Test Fishery Data:
  ## Add in test fishery counts if present:
  if(length(self$test_fishery_catch) > 0){
    tf_names <- names(self$test_fishery_catch)
    catch <- NULL
    cols <- c("Date", "pink", "sockeye", "jackchinook", "adultchinook", "chinook", "effort", "net_type")
    spp_tf <- c("pink", "sockeye", "jackchinook", "adultchinook", "chinook")
    for( i in 1:length(self$test_fishery_catch) ){
      cols_i <- cols[cols %in% names( self$test_fishery_catch[[i]] )]
      catch_i <- self$test_fishery_catch[[i]] |> 
                                    subset(as.Date(Date) %in% dates_) |> 
                                    subset(select = cols_i)
      spp_i <- spp_tf[spp_tf %in% names(catch_i)]
      for( j in 1:nrow(catch_i) ){
        catch_ij <- catch_i[j,]
        catch <- rbind(catch, data.frame(Date = catch_ij$Date, fishery = tf_names[i], net_type = catch_ij$net_type, species = spp_i, catch = as.numeric(catch_ij[,spp_i]), effort = catch_ij$effort))
      }
    }
    
    ## par_name defined by user:
    if(length(self$data_info$test_fishery_formula) > 1){ catch$par_name <- do.call(paste, c(catch[, self$data_info$test_fishery_formula], sep = "_"))
    }else{catch$par_name <- c(catch[, self$data_info$test_fishery_formula]) }
    # catch <- catch |> within(par_name <- paste(species, fishery, sep = "_"))
    catch <- catch |> within(day <- as.numeric(factor(Date)))

    ## Setup names of test fishery information:
    catch <- catch |> subset(paste(species, fishery, sep = "_") %in% self$data_info$test_fishery_spp)
    self$data_info$test_fishery_input <- unique(catch$par_name)
    
    # catch <- catch |> subset(par_name %in% self$data_info$test_fishery_input)
    catch <- catch |> within(q_index <- sapply(par_name, FUN = function(x) { which(x == self$data_info$test_fishery_input) }) )
    catch <- catch |> within(N_index <- sapply(species, FUN = function(x) { which(x == self$species_info$species_predict)}) )
    self$data_list$test_fishery_catch <- catch
  }
  ## Set prediction data for estimating prop to counts:
  pred_df <- self$data_list$length_data |> subset(!duplicated(grpIndex))
  pred_df <-  pred_df[order( pred_df$grpIndex),]
  self$data_list$pred_df <- pred_df  
  self$data_list$days <- dates_

  ## Estimate total salmon:
  salmon_counts <- self$salmon_counts |> subset(Date %in% dates_)
  self$data_list$total_salmon <- salmon_counts |> aggregate(count~Date, sum)
}

## Add Roxygen
#' @export
set_model_proportions <- function(self, formula = list()){
  if(self$ndays == 1)
    p_formula <- extractControls(formula$all, ~ -1 + stratum)
  else
    p_formula <- extractControls(formula$all, ~ -1 + stratum:factor(day))

  ## Model Data:
  ndays <- self$ndays
  dates_ <- seq(self$est_date - ndays + 1, self$est_date, 1)
  if( !identical(dates_, self$data_list$days) ) set_daily_data(self)
    
  ## Potential species models:
  species <- self$species_info$species
  species_other <- self$species_info$species_other
  ni <- 0
  self$data_list$alpha_formula <- list()
  self$data_list$X_proportions <- list()
  for( i in seq_along(species_other) ){
    sppi <- species_other[i]
    formulai <- extractControls(formula[[sppi]], p_formula)
    self$data_list$alpha_formula[[sppi]] <- formulai
    if(!all(all.vars(formulai) %in% names(self$data_list$pred_df))) 
      stop("Variables provided in formula are missing from the sonar lengths data. See 'obj$analysisData$lengthData' for potential variables.\n")
    if(any(all.vars(formulai) %in% "day") & ndays == 1)
      stop("Exclude 'day' from formula if using only one day of data.\n")      
    self$data_list$X_proportions[[sppi]] <- model.matrix(formulai, data = self$data_list$pred_df)
  }

  ## jack chinook model:
  formula_jc <- extractControls(formula[["jackchinook"]], ~ 1) ## adult chinook or small/large adults will default to the same, minus input.
  # if("all" %in% names(formula)) formula_jc <- p_formula
  self$data_list$X_proportions[["jackchinook"]]  <- model.matrix(formula_jc, data = self$data_list$pred_df)
}

set_length_adjustment <- function(self, formula = NULL){
  formula <- extractControls(formula, ~ beamWidth.cm)
  ## Beam Width Adjustment or potentially just bins, or nothing (adjustLengths = FALSE)
  self$data_list$X_length <- model.matrix(formula, data = self$data_list$length_data)
  self$data_list$length_adjust_formula <- formula
  self$fit_info$adjust_lengths <- TRUE
}

#' @export
set_species_lengths <- function(self, mu = NULL, sigma = NULL, proportions_chinook = NULL, test_fishery_lengths = NULL, ndays = 7){

  dates_ <- seq(self$est_date - ndays + 1, self$est_date, 1)

  ## Default values:
  mu_ <-  self$default_parameters$mu[self$species_info$species]
  sigma_ <- self$default_parameters$sigma[self$species_info$species]
  proportions_chinook_ <- c("jackchinook" = 0.05, "adultchinook" = 0.95,"smalladultchinook" = 0.17, "largeadultchinook" = 0.78)
  proportions_chinook_ <- proportions_chinook_[self$species_info$species_chinook]
  
  if(!all(names(mu) %in% self$species_info$species)) stop("Provided species mean size 'mu' that is not a current species or is spelled wrong.")
  if(!all(names(sigma) %in% self$species_info$species)) stop("Provided species standard deviation 'sigma' that is not a current species or is spelled wrong.")
  
  ## Add in user supplied values:
  mu_[names(mu)] <- mu
  sigma_[names(sigma)] <- sigma 
  proportions_chinook_[names(proportions_chinook)] <- proportions_chinook

  ## Now if test fishery lengths are present use those to infill mu and sigma:
  if(!is.null(test_fishery_lengths)){
    if(!all(c("FL.cm", "Date", "Species") %in% names(test_fishery_lengths)))
      stop("test_fishery_lengths dataframe must have 'Date', 'FL.cm', and 'Species'.")

    if(!is.POSIXct(test_fishery_lengths$Date) & !is.Date(test_fishery_lengths$Date))
      stop("Date must be provided in a standard date format")

    test_fishery_lengths$Date <- as.Date(test_fishery_lengths$Date)
  
    ## Do non-chinook species first:
    spp_names_mu <- setdiff(grep("resident", self$species_info$species_other, invert = TRUE, value = TRUE), names(mu))
    spp_names_sigma <- setdiff(grep("resident", self$species_info$species_other, invert = TRUE, value = TRUE), names(sigma))
    spp_names <- unique(c(spp_names_mu, spp_names_sigma))
    tfdays <- test_fishery_lengths |> subset(Date %in% dates_)
    
    for(i in seq_along(spp_names) ){
      tfl <- tfdays |> subset(grepl(spp_names[i], tolower(Species)) & !is.na(FL.cm))
      x <- tfl$FL.cm
      if(length(x) < 5) {
        cat("[Warning]  Unable to estimate", spp_names[i] , "lengths from the test fishery data. Using default mean and variance.\n")        
        next
      }
      ## If mu not provided, compute here
      if(!spp_names[i] %in% names(mu)){
        mu_[spp_names[i]] <- mean(x) 
      }
      ## If sigma not provided, compute here
      if(!spp_names[i] %in% names(sigma)){
        sigma_[spp_names[i]] <- sd(x)
      }
    }

    tfl <- tfdays |> subset(grepl("chinook", tolower(Species)) & !is.na(FL.cm))
    x <- tfl$FL.cm
    chinook_names <- self$species_info$species_chinook
    nchinook <- self$species_info$nchinook
    mu_chin <- self$default_parameters$mu[chinook_names]
    sigma_chin <- self$default_parameters$sigma[chinook_names]
    
    chinook_list <- fitChinookLengths(x, chinook_names, mu_chin, sigma_chin, proportions_chin = proportions_chinook_, 
      mu_fixed = mu, sigma_fixed = sigma)
    mu_[chinook_names] <- chinook_list$mu
    sigma_[chinook_names] <- chinook_list$sigma
    proportions_chinook_[chinook_names] <- chinook_list$p
  }
  ## Add the values to the the parameters:
  self$default_parameters$mu <- mu_
  self$default_parameters$sigma <- sigma_
  self$default_parameters$alpha_jackchinook <- log(proportions_chinook_[1]/(1-proportions_chinook_[1]))
  self$default_parameters$proportion_adultchinook <- chinook_list$p[2:nchinook]/(1-chinook_list$p[1])
}

#' @export
fitChinookLengths <- function(x, chinook_names, mu_chin, sigma_chin, proportions_chin, mu_fixed = NULL, sigma_fixed = NULL){
    nchinook <- length(chinook_names)

    ## See what values are fixed:
    mu_jack_fixed <- "jackchinook" %in% mu_fixed
    sigma_jack_fixed <- "jackchinook" %in% sigma_fixed
    
    mu_adult_fixed <- any(grepl("adultchinook", mu_fixed))
    sigma_adult_fixed <- any(grepl("adultchinook", sigma_fixed))
        
    ## Scenario 1: All fixed values:
    if(sigma_adult_fixed & mu_adult_fixed & mu_jack_fixed & sigma_jack_fixed){
      return(list(mu_chin = mu_chin, sigma_chin = sigma_chin, proportions_chin = proportions_chin))
    }

    if(any(sigma_chin < 1)) sigma_chin[sigma_chin < 1] <- 2.5 ## Update sigma_chin just in case it was sent too small.
    
    mu_est <- numeric(nchinook)
    sigma_est <- numeric(nchinook)
    p_est <- numeric(nchinook)
    
    xjack <- x[x<=50]
    xadult <- x[x > 50]
    jack_success <- FALSE
    adult_success <- FALSE
    fit <- NULL
    
    mu_fixed <- mu_fixed[grep("chinook", names(mu_fixed), value = TRUE)]
    sigma_fixed <- sigma_fixed[grep("chinook", names(sigma_fixed), value = TRUE)]
    
    if(length(sigma_fixed) == 0) sigma_fixed <- NULL
    if(length(mu_fixed) == 0) mu_fixed <- NULL
    
    if(length(xjack) > 3 & length(xadult) > 3){
      fit_all <- basicMixtureModel(x, K = nchinook, mu = mu_chin, sigma = sigma_chin, prob = proportions_chin, 
                               mu_fixed = mu_fixed, sigma_fixed = sigma_fixed, names = chinook_names, control = list(verbose = FALSE))
      if(min(fit_all$mu) < 48){
        jack_success <- TRUE
        adult_success <- TRUE
        ord <- order(fit_all$mu)
        mu_est <- fit_all$mu[ord]
        sigma_est <- fit_all$sigma[ord]
        p_est <- fit_all$proportion[ord]
      }else{
        cat("[Warning]  Chinook mixture for jack Chinook did not fit.\n")
      }
      if(nchinook == 3){
        if(fit_all$mu[3] - fit_all$mu[2] < 3) cat("[Warning]  Chinook mixture means for small and large adults are very close. Suggest choosing one size class.\n")
      }
    }

    if(length(xadult) >= 7 & nchinook == 3 & !jack_success){
      if(!"jackchinook" %in% names(mu_fixed)) mu_fixed <- c(mu_chin[1], mu_fixed)
      if(!"jackchinook" %in% names(sigma_fixed)) sigma_fixed <- c(sigma_chin[1], sigma_fixed)
      fit_adult <- basicMixtureModel(x, K = nchinook, mu = mu_chin[-1], sigma = sigma_chin[-1], prob = proportions_chin, 
                               mu_fixed = mu_fixed, sigma_fixed = sigma_fixed, names = chinook_names, control = list(verbose = FALSE))
      ord <- order(fit_adult$mu)
      mu_est <- fit_adult$mu[ord]
      sigma_est <- fit_adult$sigma[ord]
      p_est <- fit_adult$proportion[ord]
      adult_sucess <- TRUE
      if(mu_est[2] < 50){
        adult_sucess <- FALSE  ## Can't allow adult chinook to be size of jack.
        cat("[Warning]  Chinook mixture means for small adults look to be fitting the size of jacks. Suggest choosing one size class.\n")
      }
      if(abs(diff(fit_adult$mu[2:3])) < 3) cat("[Warning]  Chinook mixture means for small and large adults are very close. Suggest choosing one size class.\n")
    }
    
    if(!jack_success){
      if(length(xjack) > 3) {
        mu_est[1] <- mean(xjack)
        sigma_est[1] <- sd(xjack)
        p_est <- proportions_chin[1]
      }else{
        mu_est[1] <- mu_chin[1]
        sigma_est[1] <- sigma_chin[1]
        p_est[1] <- proportions_chin[1]
      }
    }
    if(!adult_success){
      if(length(xadult) > 3) {
        mu_est[2:nchinook] <- mean(xadult)
        sigma_est[2:nchinook] <- sd(xadult)
        p_est[2:nchinook] <- proportions_chin[2:nchinook]
      }else{
        mu_est[2:nchinook] <- mu_chin[2:nchinook]
        sigma_est[2:nchinook] <- sigma_chin[2:nchinook]
        p_est[2:nchinook] <- proportions_chin[2:nchinook]
      }
    }
    names(p_est) <- chinook_names
    names(sigma_est) <- chinook_names
    names(mu_est) <- chinook_names
    
    p_est <- p_est/sum(p_est)
    return(list(mu = mu_est, sigma = sigma_est, p = p_est))
}

#' @export
set_model_parameters <- function(self, fixed_parameters = c("mu", "sigma", "proportion_jackchinook"),
                                  fixed_values = list(), initial_values = list(),
                                  delta_mu_bounds = NULL){
  ## Set which counts apply to the test fishery analysis:
  ## *****
  ## self$analysisData$predDF$SalmonCountTF <- ifelse(self$analysisData$predDF$SonarBin %in% testFisheryBins, self$analysisData$predDF$SalmonCount, 0)
  species <- self$species_info$species

  for(i in seq_along(fixed_parameters)){
    if(!fixed_parameters[i] %in% names(fixed_values)){
      fixed_values[[fixed_parameters[i]]] <- self$default_parameters[[fixed_parameters[i]]]
      if(fixed_parameters[i] %in% c("mu", "delta_mu", "sigma")) fixed_values[[fixed_parameters[i]]] <- fixed_values[[fixed_parameters[i]]][species]
      
    }
  }

  self$params_fixed <- list()
  self$params_init <- list()
  
  if(is.null(delta_mu_bounds)) delta_mu_limits = cbind(rep(-1, self$species_info$nspecies), rep(1, self$species_info$nspecies))
  if(length(delta_mu_bounds$lower) == 1) delta_mu_limits = cbind(rep(delta_mu_bounds$lower[1], self$species_info$nspecies), rep(delta_mu_bounds$upper[1], self$species_info$nspecies))
  self$data_list$lower_delta_mu <- delta_mu_limits[,1]
  self$data_list$upper_delta_mu <- delta_mu_limits[,2]

  ## Set log catchability:
  extractParams(self, initial_values$qinv, fixed_values$qinv, self$default_parameters$qinv[self$data_info$test_fishery_input], name = "log_qinv", transform = log)  
  ## Set alpha
  nalpha <- do.call(sum, lapply(self$data_list$X_proportions[self$species_info$species_other], ncol))
  self$default_parameters$alpha <- numeric(nalpha)
  extractParams(self, initial_values$alpha, fixed_values$alpha, self$default_parameters$alpha, name = "alpha")  
  ## Set alpha jack chinook
  self$default_parameters$alpha_jackchinook <- rep(self$default_parameters$alpha_jackchinook[1], ncol(self$data_list$X_proportions[["jackchinook"]]))
  extractParams(self, initial_values$alpha_jackchinook, fixed_values$alpha_jackchinook, self$default_parameters$alpha_jackchinook, name = "alpha_jackchinook")  
  ## Set beta
  nbeta <- ncol(self$data_list$X_length)
  self$default_parameters$beta <- rep(0, nbeta)
  extractParams(self, initial_values$beta, fixed_values$beta, self$default_parameters$beta, name = "beta")
  ## Set mu
  extractParams(self, initial_values$mu, fixed_values$mu, self$default_parameters$mu[self$species_info$species], name = "mu")  
  ## Set sigma
  extractParams(self, initial_values$sigma, fixed_values$sigma, self$default_parameters$sigma[self$species_info$species], name = "log_sigma", transform = log)  
  ## Set delta_mu
  extractParams(self, initial_values$delta_mu, fixed_values$delta_mu, self$default_parameters$delta_mu[self$species_info$species], name = "logit_delta_mu", 
                transform = logitInterval, lower = delta_mu_limits[,1], upper = delta_mu_limits[,2])
  ## Set sigma0
  extractParams(self, initial_values$sigma0, fixed_values$sigma0, self$default_parameters$sigma0, name = "log_sigma0", transform = log)

  ## proportion adult chinook:
  self$params_fixed$proportion_adultchinook <- self$default_parameters$proportion_adultchinook[grep("adult",self$species_info$species_chinook, value=TRUE)]
  if("proportion_adultchinook" %in% names(fixed_values)) self$params_fixed$proportion_adultchinook <- fixed_values$proportion_adultchinook
}

#' @export
process_albion_catch <- function(self, albion_catch){
  ## Process Date:
  if(!is.POSIXct(albion_catch$CRPT_DTT)) stop("CRPT_DTT must be of type POSIXct.")
  if(!is.POSIXct(albion_catch$NET_START_OUT)) stop("CRPT_DTT must be of type POSIXct.")
  if(!is.POSIXct(albion_catch$NET_FULL_OUT)) stop("CRPT_DTT must be of type POSIXct.")
  if(!is.POSIXct(albion_catch$NET_START_IN)) stop("CRPT_DTT must be of type POSIXct.")
  if(!is.POSIXct(albion_catch$NET_FULL_IN)) stop("CRPT_DTT must be of type POSIXct.")
  
  if(is.null(albion_catch$NET_CONFIG)) stop("Require a column named NET_CONFIG to state Chinook or Chum  net.")
  
  albion_catch <- albion_catch |> subset( format(CRPT_DTT, "%Y") == format(self$est_date, "%Y"))
  
  ## Net lengths in fathoms:
  ## Chum net = 150, Chinook net = 200 fathoms
  albion_catch <- albion_catch |> within(net_length <- ifelse(grepl("chum", NET_CONFIG, ignore.case = TRUE), 150, 200))
  albion_total <- albion_catch |> aggregate(CATCH_QTY~CRPT_DTT+NET_CONFIG+SET_NO, sum)

  albion_set <- albion_catch |> by( ~ CRPT_DTT+NET_CONFIG+SET_NO, 
    function(x){data.frame(CRPT_DTT = x$CRPT_DTT[1], NET_CONFIG = x$NET_CONFIG[1], SET_NO = x$SET_NO[1], net_length = x$net_length[1],
                NET_START_IN = min(x$NET_START_IN), NET_FULL_IN = max(x$NET_FULL_IN), 
                NET_START_OUT = min(x$NET_START_OUT), NET_FULL_OUT = max(x$NET_FULL_OUT))})
  albion_set <- do.call("rbind", albion_set)

  albion_set <- albion_set |> 
    within(time_in <- as.numeric(difftime(NET_FULL_IN,  NET_START_IN, units = "mins"))) |> 
    within(time_out <- as.numeric(difftime(NET_FULL_OUT, NET_START_OUT, units = "mins"))) |>
    within(time_full <- as.numeric(difftime(NET_START_IN, NET_FULL_OUT, units = "mins"))) |>
    within(time <- time_in^0.5/2 + time_out/2 + time_full)
  albion_set <- albion_set |> within(effort <- net_length*time/1000)
  albion_set <- albion_set |> aggregate(cbind(effort, time, net_length) ~ CRPT_DTT + NET_CONFIG + SET_NO, sum)

  albion_total <- albion_total |> merge(albion_set, by = c("CRPT_DTT", "NET_CONFIG", "SET_NO"))
  names(albion_total)[names(albion_total) == "CATCH_QTY"] <- "chinook"  
  albion_total <- albion_total |> within(cpue <- chinook/effort) 
  albion_total$Date <- as.Date(albion_total$CRPT_DTT)
  albion_total <- albion_total |> within(net_type <- ifelse(grepl("SP CHINOOK", NET_CONFIG, ignore.case=TRUE), "chin", NET_CONFIG))
  albion_total <- albion_total |> within(net_type <- ifelse(grepl("VMN", NET_CONFIG, ignore.case=TRUE), "vmn", net_type))
  albion_total <- albion_total |> within(net_type <- ifelse(grepl("SP Chum", NET_CONFIG, ignore.case=TRUE), "chum", net_type))
  self$test_fishery_catch[["albion"]] <- albion_total
}

# self <- speciesComp$new(
            # species = c("largeresident", "jackchinook", "sockeye", "adultchinook"), 
            # site = "Mission", date = "2025-08-08",  ndays = 3)
# self$processData(sonar_counts = sonarCounts, sonar_lengths = sonarLengths, test_fishery_counts = list(whonnock = testFisheryCounts, albion = albion), 
  # include_tangled = TRUE, dropN = 2, salmon_passage_table = salmonCounts, feasible_lengths = c(10,120))

# test_fishery_lengths <- testFisheryLengths
# sonar_counts <- sonarCounts
# sonar_lengths <- sonarLengths
# test_fishery_counts <- testFisheryCounts
# include_tangled <- TRUE
# dropN <- 2
# salmon_passage_table <- salmonCounts
# feasible_lengths <- c(10,120)

# self$setSpeciesLengths(test_fishery_lengths = testFisheryLengths, ndays = 6)

# self$setModelParameters(fixed_parameters = c("mu", "sigma", "proportion_jackchinook", "beta"),
                                  # fixed_values = list(beta = c(0, 1)), initial_values = list(),
                                  # formula_proportions = list("all" = ~factor(day)),
                                  # delta_mu_bounds = list(lower = -2.5, upper = 2.5),
                                  # formula_lengths = ~ beamWidth.cm, adjust_lengths = TRUE)

# pars_est <- fitJointModel(self)

# (1/pars_est$q)
# pars_est$N_daily
# pars_est$mu
# pars_est$mu_adj
# pars_est$beta

# self$data_list$test_fishery_catch
# self$data_list$total_salmon

# self$params_init
# self$params_fixed

# self$data_list$X_length
# self$data_list$X_proportions
# self$species_info$species_predict

# test_fishery_catch <- data.frame(Date = rep(seq(self$est_date - self$ndays + 1, self$est_date, 1), each = 3), 
                                 # fishery = rep(c("whonnock", "whonnock", "albion"), 3), 
                                 # species = c("sockeye", "chinook", "chinook"),
                                 # catch = c(25,6,12,32,8,15,42,9,7),
                                 # effort = rep(1, 9))
# test_fishery_catch$day <- as.numeric(factor(test_fishery_catch$Date))
# test_fishery_catch$q_index <- sapply(paste(test_fishery_catch$species,test_fishery_catch$fishery, sep = "_"), FUN = function(x) { which(x == self$data_info$test_fishery_input)})
# test_fishery_catch$N_index <- sapply(test_fishery_catch$species, FUN = function(x) { which(x == self$species_info$species_predict)})




# calculateResiduals <- function(self){
  # Length Residuals:
  # pars <- self$estimatedParameters
  # mu <- pars$mu
  # muChin <- pars$muChin
  # Parameter Processing
  # Kchin <- length(muChin)
  # K0 <- length(mu)
  # K <- Kchin + K0
  # muc <- c(mu, muChin)
  # alpha <- pars$alpha
  # beta <- pars$beta
  # alphaJackChinook <- pars$alphaJackChinook
  
  # Standard deviation incld. observation error
  # logsigma_ <- c(pars$logsigma, pars$logsigmaChin)
  # sigma <- sqrt(exp(2*logsigma_) + exp(2*pars$logsigma0))

  # Set up proportions. alpha parameter for predicting proportions. 
  # Xalpha is a list of design matrices for each species.
  # np <- nrow(self$analysisData$predDF)
  # logitp <- matrix(0, nrow = np, ncol = K0) ## +1 is Chinook.
  # indx0 <- 1
  # for( i in 1:K0 ){
    # nc <- ncol(self$analysisData$Xprop[[i]]) 
    # indx1 <- indx0 + nc - 1
    # logitp[,i] <- as.matrix(self$analysisData$Xprop[[i]]) %*% alpha[indx0:indx1]
    # indx0 <- indx1 + 1
  # }

  # logitpjack <- as.matrix(self$analysisData$XpropChin) %*% alphaJackChinook
  # pjack <- 1/(1+exp(-logitpjack))

  # p <- matrix(0, nrow = np, ncol = K)
  # for( i in 1:np ) {
    # p[i, 1:(K0+1)] <- expitM(logitp[i,])  
    # p[i, (K0+1):K ] <- p[i, (K0+1)]*c(pjack[i], (1-pjack[i])*pars$pAdultChinook)
  # }

  # if(self$adjustLengths){
    # L <- self$analysisData$lengthData$L.cm.adj - as.matrix(self$analysisData$Xlength) %*% beta
  # }else{
    # L <- self$analysisData$lengthData$L.cm.adj
  # }
  # wgts <- self$analysisData$lengthData$weights
  # pobs <- p[self$analysisData$lengthData$grpIndex,]
  
  # id <- sapply(1:length(L), FUN = function(i){which.max(pobs[i,]*dnorm(L[i], muc, sigma))}) ## Definitely not this one.
  # residuals <- sapply(1:length(L), FUN = function(i){(L[i]-muc[id[i]])/sigma[id[i]]}) ## Definitely not this one.
    
  
  # self$analysisData$lengthData$residuals <- residuals
  # self$analysisData$lengthData$speciesid <- self$Species[id]
    
  # boxplot(residuals ~ speciesid, data = self$analysisData$lengthData, xlab = "Species", ylab = "Residuals")
  # abline(h = 0, col = 'red')
  # qqnorm(residuals)
  # qqline(residuals)
  # plot(x = 1:length(L), y = residuals, xlab = "Index", ylab = "Residuals", cex = self$analysisData$lengthData$weights/(median(self$analysisData$lengthData$weights)))
  # abline(h = 0, col = 'red', lty = 1)
  # abline(h = c(-1.96, 1.96), col = 'red', lty = 2)

  # Can print what the 'outlier' values are.
  # self$analysisData$lengthData |> subset(abs(residuals) > 1.96)

# }