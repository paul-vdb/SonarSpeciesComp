#' Base class for species composition estimation methods
#'
#' @description An R6 base class that is used to step through all the different setups and steps of estimating species composition on the Fraser River (Qualark and Mission) using Hydroacoustics and the test fishery.
#'
#' @field species_info List with information about species naming and ordering.
#' @field est_date Date for the day of analysis.
#' @field ndays Number of days to join.
#' @field data_list List of all data used for analysis for the ndays of interest.
#' @field data_info List of information specific to site and analysis (e.g. test fishery information).
#' @field fit_info List of optimiziation and fitting information for fitting the model.
#' @field sonar_lengths Dataframe of sonar lengths passed for entire data processed.
#' @field test_fishery_catch List of test fishery catch for each fishery.
#' @field salmon_counts Dataframe daily salmon counts estimate.
#' @field params_estimated List of all estimated parameters after fitting model.
#' @field params_fixed List all values to hold fixed when fitting the model.
#' @field params_init List all initial values to use when fitting the model.
#' @field sim_data List of daily simulated data to match data_list.
#' @field default_parameters List of default parameter values to use when inits are not provided.
#' @field prior_distributions List of prior distributions for all parameters.
#'
#' @importFrom R6 R6Class
#' @export
speciesCompModel <- R6::R6Class("SpeciesCompModel",
  public = list(
    # --- Fields ---  
    species_info = NULL,
    est_date = NULL,
    ndays = NULL,
    data_list = NULL,
    data_info = NULL,
    fit_info = NULL,
    sonar_lengths = NULL,
    test_fishery_catch = NULL,
    salmon_counts = NULL,
    params_estimated = list(),
    params_fixed = list(),
    params_init = list(),
    sim_data = NULL,
    default_parameters = NULL,
    prior_distributions = NULL,
    .prior_dens = NULL,
    .posterior_dens = NULL,
    
    #' @description Initialize the R6 object
    #' @param species vector of species names  e.g. c("largeresident", "jackChinook", "sockeye", "adultchinook").
    #' @param date Date that you will estimate. Place holder but does not need to be set.
    #' @param ndays Number of days to combine, defaults to 1.
    #' @param site Choose between "Mission" or "Qualark", roll angle is set be default based on site.
    #' @param roll_angle The roll angle that the sonar is set to (default = NULL).
    #' @param feasible_lengths Lower and upper bound of feasible lengths that can be measured on the sonar to remove infeasible (default = c(10, 120)).    
    initialize = function(species = c("jackchinook", "sockeye", "adultchinook"), date = NULL, ndays = 1, site = NULL, roll_angle = NULL, feasible_lengths = c(10, 120)) {
      if( missing(site) ) stop("Need to initialize with a site = 'Mission' or 'Qualark'")
      if( !site %in% c("Mission", "Qualark") ) stop("Need to initialize with a site = 'Mission' or 'Qualark'")
      self$species_info <- speciesCheck(species)
      if(is.null(date)) date <- Sys.Date()
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
        self$data_info$test_fishery_spp <- c("sockeye_whonnock", "chinook_whonnock", "chinook_albion")
      }
      if(site == "Qualark"){
        self$data_info$test_fishery_spp <- c("sockeye_qualark", "chinook_qualark")
      }
      self$data_info$test_fishery_formula <- ~ 0 + species:fishery
      
      self$test_fishery_catch <- list()
      self$default_parameters <- default_params()
      default_prior <- \(...){0}
      self$prior_distributions <- list(dlog_qinv = default_prior, dbeta = default_prior, dalpha_jackchinook = default_prior, 
                                       dlog_sigma = default_prior, dmu = default_prior, dlogit_delta_mu = default_prior, dlog_sigma0 = default_prior,
                                       dsmallresident = default_prior, dlargeresident = default_prior)
    },
    #' @description Update the R6 object with species names to fit in a model.
    #' @param species Character vector with options are "smallresident", "largeresident", "pink", "sockeye", "coho", "chum", "jackchinook", "adultchinook", "smalladultchinook", "largeadultchinook"
    setSpecies = function(species = NULL){
      if(!is.null(species)){
        if(!all(species %in% self$species_info$species)){
          cat("Resetting default Parameters.\n")
          self$default_parameters <- default_params()
        }
        self$species_info <- speciesCheck(species)
      }
    },
    #' @description Update the estimation date and number of days to use. 
    #' @param date Date to change the analysis date to.
    #' @param ndays Number of days to combine.
    setDate = function(date = NULL, ndays = NULL){
      if(!is.null(date)) self$est_date <- checkDate(date)
      if(!is.null(ndays)) self$ndays <- ndays
    },
    #' @description Control optimization of the EM algorithm. 
    #' @param control List with values `$tolerance` what difference is okay (default 1e-8), `$maxiters` maximum iterations of EM algorithm (default 1000), `$relative_difference` true or false to choose relative vs absolute difference on each iteration of algorithm, `$verbose` true or false to print details of EM algorithm.
    controlOptimization = function(control = list()){
      self$fit_info$tolerance <-  extractControls(control$tolerance, 1e-8)
      self$fit_info$maxiters <-  extractControls(control$maxiters, 1000)
      self$fit_info$relative_difference<- extractControls(control$relative_difference, 1000)
      self$fit_info$verbose <- extractControls(control$verbose, FALSE)
    },
    #' @description Function to process input data for analysis. Does the main data processing and adds a few columns needed elsewhere for the entire dataset added to make it fast and easy to do different days.
    #' @param sonar_counts Data frame of sonar counts in standard format.
    #' @param sonar_lengths Data frame of sonar lengths in standard format.
    #' @param test_fishery_counts List of dataframes named after the test fishery (e.g. 'whonnock' or 'albion')
    #' @param include_tangled Logical whether to use the tangled fish in the gilnet. 
    #' @param dropN Minimum number of observations to keep length measurements per hour/bin/bank etc. if > dropN (default 2) but not always necessary. 
    #' @param salmon_passage_table Data frame that has the total number of salmon predicted each day at the site.
    #' @param feasible_lengths Minimum and maximum lengths that are feasible. Set by default but can be changed here.
    processData = function(sonar_counts, sonar_lengths, test_fishery_counts = NULL, include_tangled = TRUE, dropN = 2, salmon_passage_table = NULL, feasible_lengths = NULL){
      if(!is.null(feasible_lengths)) self$data_info$feasible_lengths <- feasible_lengths
      if(!is.null(test_fishery_counts)) self$fit_info$include_test_fishery <- TRUE

      if(self$data_info$site == "Qualark") {
        process_qualark_data(self, net_length = 30)
      }
      if(self$data_info$site == "Mission") {
        process_mission_lengths(self, sonar_counts, sonar_lengths, dropN)
        if("whonnock" %in% names(test_fishery_counts)) process_mission_catch(self, test_fishery_counts[["whonnock"]], name = "whonnock", tangled = include_tangled)
        if("albion" %in% names(test_fishery_counts)) process_albion_catch(self, test_fishery_counts[["albion"]] )
        if("brownsvillebar" %in% names(test_fishery_counts)) process_mission_catch(self, test_fishery_counts[["brownsvillebar"]], name = "brownsvillebar", tangled = include_tangled)        
        if(!is.null(salmon_passage_table)){
          if(any("TotalSalmon Official" == colnames(salmon_passage_table))) {
            salmon_passage_table$Date <- as.Date(salmon_passage_table$MissionDate)
            suppressWarnings(newtab <- data.frame(Date = salmon_passage_table$Date, count = as.numeric(salmon_passage_table$`TotalSalmon Official`)) |> subset(!is.na(count)))
            self$salmon_counts <- newtab
          }else{ 
            process_mission_salmon_passage(self, salmon_passage_table)
          }
        }
      }
    },
    #' @description Set species length parameter information 
    #' @param mu Named vector the user wants to set as default values \code{c("largeresident" = 38)}.
    #' @param sigma Named vector the user wants to set as default values e.g. \code{c("largeresident" = 2.5)}.
    #' @param proportion_adultchinook Named vector that sets default proportions of Chinook e.g. \code{c("jackchinook" = 0.05, "adultchinook" = 0.95)}.
    #' @param test_fishery_lengths Data frame that contains test fishery lengths for use of parameter estimation default values.
    #' @param ndays Number of prior days to use from test fishery lengths (default = 6).
    #' @param plot Logical to plot the estimated mean and variance set by the function.
    #' @param verbose Logical to print warning messages or not.
    setSpeciesLengths = function(mu = NULL, sigma = NULL, proportion_adultchinook = NULL, test_fishery_lengths = NULL, ndays = 6, plot = FALSE, verbose = FALSE){
      set_species_lengths(self, mu, sigma, proportion_adultchinook, test_fishery_lengths, ndays, verbose = verbose)
      if(plot) plot_test_fishery_lengths(self, test_fishery_lengths = test_fishery_lengths, ndays = ndays)
    },
    #' @description Set the model parameters before fitting with EM algorithm, including which values to fix and which to fit.
    #' @param fixed_parameters Which parameters to hold fixed (default \code{c("mu", "sigma", "proportion_jackchinook")}).
    #' @param fixed_values List provided by the user to hold particular values fixed, if nothing provided but parameters are held fixed then default values will be used.
    #' @param initial_values List provided by the user to improve fitting time for estimated parameters, default are used if not provided.
    #' @param formula_proportions list of formulas to use for alpha for each species, if \code{list("all" = ~ factor(day))}, then each species is given the same formula, except jack chinook are treated separately.
    #' @param formula_lengths Formula for how to deal with beam spreading (default = \code{~beamWidth.cm}, 
    #' @param date Date (optional) to change analysis date.
    #' @param ndays Number of days to combine (optional).
    #' @param delta_mu_bounds Matrix or vector to schoose how much the mean length of each species can vary. Matrix must have nrows of number of species. Otherwise a vector is length 2 and declares lower and upper for all.
    #' @param expansion_formula Formula to set the test fishery expansion line relationships (defualt = ~ species:fishery).
    #' @param test_fishery_spp Vector of joined species and test fishery names to used e.g. ("sockeye_whonnock", "adultchinook_whonnock", "chinook_albion").
    setModelParameters = function(fixed_parameters = c("mu", "sigma", "proportion_jackchinook"),
                                  fixed_values = list(), initial_values = list(),
                                  formula_proportions = list(),
                                  formula_lengths = ~ beamWidth.cm,
                                  date = NULL, ndays = NULL, delta_mu_bounds = NULL,
                                  expansion_formula = NULL,
                                  test_fishery_spp = NULL){
      ## If user wants to update the date and species, they can do that here:
      if(!is.null(expansion_formula)) self$data_info$test_fishery_formula <- expansion_formula
      if(!is.null(test_fishery_spp)) self$data_info$test_fishery_spp <- test_fishery_spp

      self$setDate(date, ndays)
      set_daily_data(self)

      ## Set test fishery formula
      self$data_list$X_test_fishery <- model.matrix(self$data_info$test_fishery_formula, data = self$data_list$test_fishery_catch)
      ## Remove columns that are all zero (e.g. sockeye has 1 net type).
      csum <- colSums(abs(self$data_list$X_test_fishery)) 
      self$data_list$X_test_fishery <- self$data_list$X_test_fishery[, csum > 0, drop = FALSE]
      colnames(self$data_list$X_test_fishery) <- gsub("species|net_type|fishery", "", colnames(self$data_list$X_test_fishery))
      colnames(self$data_list$X_test_fishery) <- gsub(":", "_", colnames(self$data_list$X_test_fishery))

      set_length_adjustment(self, formula_lengths)
      set_model_proportions(self, formula_proportions)
      set_model_parameters(self, fixed_parameters, fixed_values, initial_values, delta_mu_bounds)      
    },
    #' @description EM algorithm for fitting the joint species composition model.
    #' @param control List containing \code{controlOptimization} arguments as well as '$include_test_fishery' - logical to include or exclude test fishery catch information in the model, '$adjust_lengths' logical include beam spreading adjusment, and '$test_fishery_weights' the model weights for the test fishery.
    #' @return List of parameter estimates, which is also stored in the \code{params_estimated} field of the R6 object.
    fitModel = function(control = list()){
      self$controlOptimization(control)
      self$fit_info$include_test_fishery <- extractControls(control$include_test_fishery, self$fit_info$include_test_fishery)
      self$data_info$test_fishery_weights <- extractControls(control$test_fishery_weights, self$data_info$test_fishery_weights)
      self$fit_info$adjust_lengths <- extractControls(control$adjust_lengths, self$fit_info$adjust_lengths)
    
      fit_joint_model(self)
    },
    #' @description Simulate the model for the counts on the 'est_date' and 'ndays' with parameter values equal to 'params_estimated' within the object.
    simulate = function(){
      simulate_given_data(self)
    },
    #' @description Define Prior distributions (or penalty functions) for fitting the joint species composition model.
    #' @param priors List of prior distributions. Defaults to \code{\(...){0}}.
    #' @param includeJacobian Logical of whether or not to include jacobian transformation.
    setPriors = function(priors = list(), includeJacobian = TRUE){
      set_priors(self, priors = priors, includeJacobian)
    },
    #' @description Plot the fitted mixture model against the historgram of the length data.
    #' @param day to plot (day = 1 is `est_date`) (default = 1).
    #' @param ... Additional plot arguments.
    plotMix = function(day = 1, ...){
      plot_mix(self, day = day, ...)
    },
    #' @description Plot Pearson residuals for the test fishery data (CPUE - Nq) against each species, test fishery, and net type.
    #' @param includePrior Should the prior distribution be plotted alongside the estimates.
    plotTestFishery = function(includePrior = TRUE){
      plot_test_fishery(self, includePrior = includePrior)
    },
    #' @description Plot estimated beam spreading effects against sonar range.
    #' @param includePrior Should the prior distribution be plotted alongside the estimates.
    plotBeamSpreading = function(includePrior = TRUE){
      plot_beam_spreading(self, includePrior = includePrior)
    },
    #' @description Plot test fishery lengths against the assumed mean and standard deviation to be used in the model.
    #' @param test_fishery_lengths Length data frame that contains Date, FL.cm, and Species.
    #' @param ndays Number of days to join for plotting test fishery data.
    plotTestFisheryLengths = function(test_fishery_lengths = NULL, ndays = 6){
      plot_test_fishery_lengths(self, test_fishery_lengths = test_fishery_lengths, ndays = ndays)
    },
    priorCheck = function(parameter){
      if(parameter %in% c("sigma", "mu", "sigma0", "qinv")) parameter <- paste0("log_", parameter)
      if(parameter == "delta_mu") parameter <- paste0("logit_delta_mu")
      plot_prior(self, parameter)
    }
  )
)

#' @export
default_params <- function(){

  .default_parameters <- list()
  if(!is.null(.GlobalEnv$default_parameters)){
    .default_parameters <- .GlobalEnv$default_parameters
  }
  if(is.null(.default_parameters$alpha)) .default_parameters$alpha <- c(0,0,0,0)
  if(is.null(.default_parameters$alpha_jackchinook)) .default_parameters$alpha_jackchinook <- 0
  if(is.null(.default_parameters$beta)) .default_parameters$beta <- c(0, 0)
  if(is.null(.default_parameters$mu)) .default_parameters$mu <- c("smallresident" = 23, "largeresident" = 35, "pink" = 50, "sockeye" = 58, "coho" = 65, "chum" = 70, "jackchinook" = 40, "smalladultchinook" = 62, "largeadultchinook" = 80, "adultchinook" = 70)
  if(is.null(.default_parameters$sigma)) .default_parameters$sigma <- c("smallresident" = 2.5, "largeresident" = 2.5, "pink" = 4, "sockeye" = 3.5, "coho" = 3.5, "chum" = 3, "jackchinook" = 2.5, "smalladultchinook" = 4, "largeadultchinook" = 6, "adultchinook" = 7)
  if(is.null(.default_parameters$sigma0)) .default_parameters$sigma0 <- 4.5
  if(is.null(.default_parameters$delta_mu)) .default_parameters$delta_mu <- c("smallresident" = 0, "largeresident" = 0, "pink" = 0, "sockeye" = 0, "coho" = 0, "chum" = 0, "jackchinook" = 0, "smalladultchinook" = 0, "largeadultchinook" = 0, "adultchinook" = 0)
  if(is.null(.default_parameters$proportion_adultchinook)) .default_parameters$proportion_adultchinook <- c("smalladultchinook" = 0.2, "largeadultchinook" = 0.8, "adultchinook" = 1)
  .default_parameters
}

#' @export
set_default_expansion <- function(tfnames){
  ## Expansion Line Albion: 3049 (2020 Brittany Jenewein)
  init <- numeric(length(tfnames))
  names(init) <- tfnames
  init[grep("sockeye", tfnames)] <- 10000
  init[grep("chinook", tfnames)] <- 5000
  init[grep("adultchinook", tfnames)] <- 4500
  init[grep("jackchinook", tfnames)] <- 15000
  init[grep("coho", tfnames)] <- 2548
  init[grep("chum", tfnames)] <- 3201
  init[grep("pink", tfnames)] <- 20000

  init
}

#' @export
set_priors <- function(self, priors = list(), includeJacobian = TRUE){        
  ## Make Tape and use Atomic here:
  self$prior_distributions <- list()
    
  self$prior_distributions$dlog_qinv <- addPrior(priors$qinv, \(x){log(abs(x))}, includeJacobian & !is.null(priors$qinv))
  self$prior_distributions$dlog_sigma <- addPrior(priors$sigma, \(x){log(abs(x))}, includeJacobian & !is.null(priors$sigma))
  self$prior_distributions$dlog_sigma0 <- addPrior(priors$sigma0,\(x){log(abs(x))}, includeJacobian & !is.null(priors$sigma0))
  self$prior_distributions$dlogit_delta_mu <- addPrior(priors$delta_mu, logDetJac_logitInterval, includeJacobian & !is.null(priors$delta_mu))
  self$prior_distributions$dmu <- addPrior(priors$mu, NULL, FALSE)
  self$prior_distributions$dbeta <- addPrior(priors$beta, NULL, FALSE)
  self$prior_distributions$dalpha_jackchinook <- addPrior(priors$alpha_jackchinook, NULL, FALSE)
  self$prior_distributions$dlargeresident <- addPrior(priors$N_largeresident, NULL, FALSE)
  self$prior_distributions$dsmallresident <- addPrior(priors$N_smallresident, NULL, FALSE)
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
  sonar_lengths <- sonar_lengths |> subset(nLengths > dropN)

  ## I don't like SonarAim being called 3 on right bank if only one aim.
  sonar_lengths <- sonar_lengths |> within(SonarAimF <- ifelse(SonarBank == "Right Bank", "Aim3", SonarAim))

  self$sonar_lengths <- sonar_lengths
}


#' Process PSC Catch
#'
#' @param self R6 speciesCompModel Object.
#' @param test_fishery_counts Data frame in standard format for Whonnock catch.
#' @param name Test fishery name (default = "whonnock").
#' @param tangle Logical whether to include or exclude tangled counts.
#'
#' @details input data frame should have columns start_out, full_out, start_in, full_in along with TRIP_DTT, and net_length. The other 
#' required columns are generally `Sockeye Adult Gilled` and `Sockeye Adult Tangled`.
#'
#' @return No object returned, appends test fishery catch to \code{self$test_fishery_catch[[name]]}
#'
#' @export
process_mission_catch <- function(self, test_fishery_counts, name = "whonnock", tangled = TRUE){

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

  out <- test_fishery_counts[,c("TRIP_DTT", "FE_SET_NO", "effort", "net_length")]
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

  tmp <- data.frame("coho" = test_fishery_counts$`Coho All Gilled` + tangled*test_fishery_counts$`Coho All Tangled`)
  out <- cbind(out, tmp)

  tmp <- data.frame("chum" = test_fishery_counts$`Chum All Gilled` + tangled*test_fishery_counts$`Chum All Tangled`)
  out <- cbind(out, tmp)

  out <- out |> within(chinook <- adultchinook + jackchinook)
  
  out <- out |> aggregate(cbind(pink, sockeye, coho, chum, jackchinook, adultchinook, chinook, effort) ~ FE_SET_NO + Date, sum)

  out$net_type <- "vmn"
  
  if(name == "brownsvillebar") out$Date <- out$Date + 1 ## Add a day to the dataset to match it roughly to Mission.

  self$test_fishery_catch[[name]] <- out
}

#' Process Mission Salmon Passage
#'
#' @param self R6 speciesCompModel Object.
#' @param salmon_passage_table Data frame in standard format for Mission count based on Right and Left Bank.
#'
#' @return No object returned, appends salmon counts \code{self$salmon_counts}
#'
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

#' Set Daily Data
#'
#' @param self R6 speciesCompModel Object.
#'
#' @details Process daily data for the model to prepare for analysis.
#'
#' @return No object returned, appends a data list with all the daily data needed for the model \code{self$data_list}
#'
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
      if(nrow(catch_i) == 0) next;
      spp_i <- spp_tf[spp_tf %in% names(catch_i)]
      for( j in 1:nrow(catch_i) ){
        catch_ij <- catch_i[j,]
        catch <- rbind(catch, data.frame(Date = catch_ij$Date, fishery = tf_names[i], net_type = catch_ij$net_type, species = spp_i, catch = as.numeric(catch_ij[,spp_i]), effort = catch_ij$effort))
      }
    }
    catch <- catch |> within(day <- as.numeric(factor(Date)))
    if(any(catch$effort == 0)){
      cat("[Warning]  Removing test fishery observations from due to zero effort (recorded) for those observations.\n")
      catch <- catch |> subset(effort > 0)
    }

    ## Setup names of test fishery information:
    catch <- catch |> subset(paste(species, fishery, sep = "_") %in% self$data_info$test_fishery_spp)
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
  if(self$data_info$site == "Mission" & "grp" %in% names(salmon_counts)) salmon_counts <- salmon_counts |> subset(!grepl("AIM", grp))
  if(any(duplicated(salmon_counts$Date))){ self$data_list$total_salmon <- salmon_counts |> aggregate(count~Date, sum)
  }else{ self$data_list$total_salmon <- salmon_counts }
}

#' Set Model Proportions
#'
#' @param self R6 speciesCompModel Object.
#' @param formula List of formulas for proportions of different species.
#'
#' @details If formula in list \code{list("all" = ~factor(day))} then this is applied to all species except Chinook. Default formula is currently
#' \code{~-1+stratum:factor(day)}.
#'
#' @return No object returned, appends a data list with all the daily data needed for the model \code{self$data_list}
#'
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
  self$data_list$alpha_formula[["jackchinook"]] <- formula_jc
  # if("all" %in% names(formula)) formula_jc <- p_formula
  self$data_list$X_proportions[["jackchinook"]]  <- model.matrix(formula_jc, data = self$data_list$pred_df)
}

#' Set Length Adjustment
#'
#' @param self R6 speciesCompModel Object.
#' @param formula Formula for how to model beam spreading.
#'
#' @details Default formula is \code{~ beamWidth.cm} will use mean correction of L.cm - beta_0 - beta_1 * beamWidth.cm. Other suggestion
#' is to use SonarBin. Any columns available within \code{self$data_list$length_data} is available here.
#'
#' @return No object returned, appends a data list with all the daily data needed for the model \code{self$data_list$X_length}
#'
#' @export
set_length_adjustment <- function(self, formula = NULL){
  formula <- extractControls(formula, ~ beamWidth.cm)
  ## Beam Width Adjustment or potentially just bins, or nothing (adjustLengths = FALSE)
  self$data_list$X_length <- model.matrix(formula, data = self$data_list$length_data)
  self$data_list$length_adjust_formula <- formula
  self$fit_info$adjust_lengths <- TRUE
}

#' Set Species Length
#'
#' @param self R6 speciesCompModel Object.
#' @param mu Named vector of lengths to assume in the model e.g. \code{c("largeresident" = 35)}.
#' @param sigma Named vector of lengths to assume in the model e.g. \code{c("largeresident" = 2.5)}.
#' @param proportions_chinook Proportion of Chinook that are different age/size classes ("jackchinook", "adultchinook").
#' @param test_fishery_lengths Data frame with columns ("FL.cm", "Date", "Species").
#' @param ndays Number of days to include for estimating expected mean and variance of lengths for each species (including current date).
#'
#' @details If a mu value is not provided and test_fishery_lengths are, then we will try and estimate mean and variance based on the 
#' number of days provided (default = 7). If not enough catches are available, or test fishery data is not provided, then default values are used.
#' For Chinook, if proportions_chinook is not provided, we attempt to estimate this through a mixture model.
#'
#' @return No object returned, appends a intial values to model object `self`, \code{self$default_parameters}.
#'
#' @export
set_species_lengths <- function(self, mu = NULL, sigma = NULL, proportions_chinook = NULL, test_fishery_lengths = NULL, ndays = 7, verbose = FALSE){

  dates_ <- seq(self$est_date - ndays + 1, self$est_date, 1)

  ## Default values:
  mu_ <-  self$default_parameters$mu[self$species_info$species]
  sigma_ <- self$default_parameters$sigma[self$species_info$species]
  proportions_chinook_ <- c("jackchinook" = 0.05, "adultchinook" = 0.95,"smalladultchinook" = 0.17, "largeadultchinook" = 0.78)
  proportions_chinook_ <- proportions_chinook_[self$species_info$species_chinook]
  
  if(!all(names(mu) %in% self$species_info$species)){
    missing <- setdiff(names(mu), self$species_info$species)
    mu <- mu[names(mu) %in% self$species_info$species]
    if(verbose) cat("[Warning]  Ignoring user supplied mu values:", missing, ". Either missing species or spelling wrong.\n")
  }
  if(!all(names(sigma) %in% self$species_info$species)){
    missing <- setdiff(names(sigma), self$species_info$species)
    sigma <- sigma[names(sigma) %in% self$species_info$species]
    if(verbose) cat("[Warning]  Ignoring user supplied sigma values:", missing, ". Either missing species or spelling wrong.\n")
  }
  
  ## Add in user supplied values:
  mu_[names(mu)] <- mu
  sigma_[names(sigma)] <- sigma 
  proportions_chinook_[names(proportions_chinook)] <- proportions_chinook

  ## Now if test fishery lengths are present use those to infill mu and sigma:
  if(!is.null(test_fishery_lengths)){
    if(!all(c("FL.cm", "Date", "Species") %in% names(test_fishery_lengths)))
      stop("test_fishery_lengths dataframe must have 'Date', 'FL.cm', and 'Species'.")

    if(!inherits(test_fishery_lengths$Date, "POSIXct") & !is.Date(test_fishery_lengths$Date))
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
        if(verbose) cat("[Warning]  Unable to estimate", spp_names[i] , "lengths from the test fishery data. Using default mean and variance.\n")        
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
      mu_fixed = mu, sigma_fixed = sigma, verbose = verbose)
    mu_[chinook_names] <- chinook_list$mu
    sigma_[chinook_names] <- chinook_list$sigma
    proportions_chinook_[chinook_names] <- chinook_list$p
  }
  ## Add the values to the the parameters:
  self$default_parameters$mu <- mu_
  self$default_parameters$sigma <- sigma_
  self$default_parameters$alpha_jackchinook <- log(proportions_chinook_[1]/(1-proportions_chinook_[1]))
  self$default_parameters$proportion_adultchinook <- proportions_chinook_[self$species_info$species_chinook[-1]]/(1-proportions_chinook_[self$species_info$species_chinook[1]])
}

#' Fit Chinook Lengths
#'
#' @param x Vector of Chinook length values.
#' @param mu_chin Named vector with mean chinook values assumed.
#' @param sigma Named vector with standard deviation of chinook values assumed.
#' @param proportions_chin Named vector of expected proportion of each Chinook age/size class.
#' @param mu_fixed Named vector of mean lengths that were provided by the user to hold fixed.
#' @param sigma_fixed Named vector of standard deviation of lengths that were provided by the user to hold fixed.
#'
#' @return List with mean, standard deviation, and proportions for Chinook.
#'
#' @export
fitChinookLengths <- function(x, chinook_names, mu_chin, sigma_chin, proportions_chin, mu_fixed = NULL, sigma_fixed = NULL, verbose = FALSE){
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
    mu_chin[names(mu_fixed)] <- mu_fixed
    sigma_chin[names(sigma_fixed)] <- sigma_fixed
    
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
        if(verbose) cat("[Warning]  Chinook mixture for jack Chinook did not fit.\n")
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
        if(verbose) cat("[Warning]  Chinook mixture means for small adults look to be fitting the size of jacks. Suggest choosing one size class.\n")
      }
      if(abs(diff(fit_adult$mu[2:3])) < 3) cat("[Warning]  Chinook mixture means for small and large adults are very close. Suggest choosing one size class.\n")
    }
    
    if(!jack_success){
      if(length(xjack) > 3) {
        mu_est[1] <- mean(xjack)
        sigma_est[1] <- sd(xjack)
        p_est <- proportions_chin[1]
      }else{
        if(verbose) cat("[Warning]  Unable to estimate jack chinook lengths from the test fishery data. Using default mean and variance.\n")        
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
        if(verbose) cat("[Warning]  Unable to estimate adult chinook lengths from the test fishery data. Using default mean and variance.\n")        
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

#' Set Model Parameters
#'
#' @param self R6 speciesCompModel Object.
#' @param fixed_parameters Character vector with entire parameters to hold fixed.
#' @param fixed_values List with names equal to parameters and named vectors providing fixed values.
#' @param initial_values List with names equal to parameters and named vectors providing initial values.
#' @param delta_mu_bounds Matrix (first column lower bound, second upper bound) or vector (lower bound, upper bound) (Defaults to (-1,1)).
#'
#' @details Set initial values and fixed values for running the species composition model. Intended to be updated before \code{self$fitModel}.
#'
#' @return No object returned, appends \code{self$params_fixed} and \code{self$params_init}
#'
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
  tf_names <- colnames(self$data_list$X_test_fishery)
  qdefault <- set_default_expansion(tf_names)
  extractParams(self, initial_values[["qinv"]], fixed_values[["qinv"]], qdefault, name = "log_qinv", transform = log)  
  ## Set alpha
  nalpha <- do.call(sum, lapply(self$data_list$X_proportions[self$species_info$species_other], ncol))
  self$default_parameters$alpha <- numeric(nalpha)
  extractParams(self, initial_values[["alpha"]], fixed_values[["alpha"]], self$default_parameters$alpha, name = "alpha")  
  ## Set alpha jack chinook
  self$default_parameters[["alpha_jackchinook"]] <- rep(self$default_parameters$alpha_jackchinook[1], ncol(self$data_list$X_proportions[["jackchinook"]]))
  extractParams(self, initial_values$alpha_jackchinook, fixed_values$alpha_jackchinook, self$default_parameters$alpha_jackchinook, name = "alpha_jackchinook")  
  ## Set beta
  nbeta <- ncol(self$data_list$X_length)
  if(length(self$default_parameters$beta) != nbeta) self$default_parameters$beta <- rep(0, nbeta)
  extractParams(self, initial_values[["beta"]], fixed_values[["beta"]], self$default_parameters[["beta"]], name = "beta")
  ## Set mu
  extractParams(self, initial_values[["mu"]], fixed_values[["mu"]], self$default_parameters$mu[self$species_info$species], name = "mu")  
  ## Set sigma
  extractParams(self, initial_values[["sigma"]], fixed_values[["sigma"]], self$default_parameters$sigma[self$species_info$species], name = "log_sigma", transform = log)  
  ## Set delta_mu
  extractParams(self, initial_values[["delta_mu"]], fixed_values[["delta_mu"]], self$default_parameters$delta_mu[self$species_info$species], name = "logit_delta_mu", 
                transform = logitInterval, lower = delta_mu_limits[,1], upper = delta_mu_limits[,2])
  ## Set sigma0
  extractParams(self, initial_values[["sigma0"]], fixed_values[["sigma0"]], self$default_parameters[["sigma0"]], name = "log_sigma0", transform = log)

  ## proportion adult chinook:
  self$params_fixed$proportion_adultchinook <- self$default_parameters$proportion_adultchinook[grep("adult",self$species_info$species_chinook, value=TRUE)]
  if("proportion_adultchinook" %in% names(fixed_values)) self$params_fixed$proportion_adultchinook <- fixed_values$proportion_adultchinook
}

#' Process Albion Catch
#'
#' @param self R6 speciesCompModel Object.
#' @param albion_catch Data frame with columns "CRPT_DTT", "NET_START_OUT", "NET_FULL_OUT", "NET_START_IN", "NET_FULL_IN", and "NET_CONFIG".
#'
#' @return No object returned, appends test fishery catch to `self`, \code{self$test_fishery_catch[["albion"]]}.
#'
#' @export
process_albion_catch <- function(self, albion_catch){
  ## Process Date:
  if(!inherits(albion_catch$CRPT_DTT, "POSIXct")) stop("CRPT_DTT must be of type POSIXct.")
  if(!inherits(albion_catch$NET_START_OUT, "POSIXct")) stop("CRPT_DTT must be of type POSIXct.")
  if(!inherits(albion_catch$NET_FULL_OUT, "POSIXct")) stop("CRPT_DTT must be of type POSIXct.")
  if(!inherits(albion_catch$NET_START_IN, "POSIXct")) stop("CRPT_DTT must be of type POSIXct.")
  if(!inherits(albion_catch$NET_FULL_IN, "POSIXct")) stop("CRPT_DTT must be of type POSIXct.")
  
  if(is.null(albion_catch$NET_CONFIG)) stop("Require a column named NET_CONFIG to state Chinook or Chum  net.")
  
  albion_catch <- albion_catch |> subset( format(CRPT_DTT, "%Y") == format(self$est_date, "%Y"))
  
  ## Net lengths in fathoms:
  ## Chum net = 150, Chinook net = 200 fathoms
  # albion_catch <- albion_catch |> within(net_length <- ifelse(grepl("chum", NET_CONFIG, ignore.case = TRUE), 150, 200))
  albion_catch <- albion_catch |> within(net_length <- 200)
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
