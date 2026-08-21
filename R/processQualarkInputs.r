
#' Process Qualark Sonar Lengths
#'
#' @param self R6 speciesCompModel Object.
#' @param sonar_counts Sonar Lengths in standard Mission format.
#' @param sonar_lengths Sonar Lengths in standard Mission Format.
#' @param dropN Minimum sample size in a stratum.
#'
#' @return No object returned, appends sonar_lengths to self.
#'
#' @export
process_qualark_lengths <- function(self, sonar_counts, sonar_lengths, dropN = 2){
  
  year <- as.numeric(format(self$est_date, "%Y"))
  
  ## Check names in dataframes are present:
  count_col_names <- c("Date", "Count Hour", "Duration", "SonarBin", "SonarBank",
                     "Up", "Down", "Net Up", "Salmon/Hour")
  length_col_names <- c("file_date", "file_time_stamp", "frequency_setting", "length_cm", 
                 "range_m", "motion", "theta_deg", "bank")
                 
  namecheck <- which(!count_col_names %in% names(sonar_counts))
  if(length(namecheck)) stop(paste0("We require column names in sonarCounts: ", paste(count_col_names[namecheck], collapse = ",")))
  namecheck <- which(!length_col_names %in% names(sonar_lengths))
  if(length(namecheck)) stop(paste0("We require column names in sonarLengths: ", paste(length_col_names[namecheck], collapse = ",")))

  ## Clean up length data:
  names(sonar_lengths)[names(sonar_lengths) == "file_date"] <- "Date"
  names(sonar_lengths)[names(sonar_lengths) == "length_cm"] <- "L.cm"
  names(sonar_lengths)[names(sonar_lengths) == "range_m"] <- "R.m"
  names(sonar_lengths)[names(sonar_lengths) == "bank"] <- "SonarBank"

  if(!is.Date(sonar_lengths$Date)) sonar_lengths$Date <- as.Date(sonar_lengths$Date, format = "%Y-%m-%d")
  sonar_lengths$Hour <-   as.numeric(substr(sonar_lengths$`file_time_stamp`, 1, 2))
  
  sonar_lengths$SonarBin <- ifelse(grepl("high", sonar_lengths$frequency_setting, ignore.case = TRUE), "Bin1", "Bin2")

  ## Subset Count data to match length data:
  sonar_lengths <- sonar_lengths |> subset(as.numeric(format(Date, "%Y")) == year)
  sonar_lengths <- sonar_lengths |> subset(L.cm >= self$data_info$feasible_lengths[1] & L.cm <= self$data_info$feasible_lengths[2])
  
  ## Create a look up code for fast merging:
  sonar_lengths <- sonar_lengths |> within(join_id <- factor(paste(SonarBank, SonarBin, Hour, Date, sep = ""))) |>
                                  within(lookup_code <- as.numeric(join_id))
  sonar_lengthsN <- sonar_lengths |> aggregate( L.cm ~ lookup_code, FUN = length)
  names(sonar_lengthsN)[names(sonar_lengthsN) == 'L.cm'] <- 'nLengths'
  sonar_lengths <- sonar_lengths |> merge(sonar_lengthsN, all.x = TRUE, all.y = TRUE)
  
  ## Subset Count data to match length data:
  names(sonar_counts)[names(sonar_counts) == "Salmon/Hour"] <- "SalmonFlux"
  names(sonar_counts)[names(sonar_counts) == "Count Hour"] <- "Hour"
  names(sonar_counts)[names(sonar_counts) == "Duration"] <- "MinsCounted"

  sonar_counts <- sonar_counts |> within(Count <- Up + Down) |> within(SalmonCount <- Up - Down) |> subset(!is.na(Count))
  
  ## check Date and Mission Date: Format should be m/d/Y.
  if(!is.Date(sonar_counts$Date)) sonar_counts$Date <- as.Date(sonar_counts$Date, format = "%Y-%m-%d")

  suppressMessages(
    sonar_counts <- sonar_counts |> within(join_id <- factor(paste(SonarBank, SonarBin, Hour, Date, sep = "")), levels = levels(sonar_lengths$join_id)) |>
                   subset(!is.na(join_id)) |> within(lookup_code <- as.numeric(join_id))
  )
  ## Merge count data to lengths to get correction to est proportion with test fishery count and stratum weights:
  sonar_lengths <- sonar_lengths |> merge(sonar_counts[, c("lookup_code", "Count", "SalmonCount", "MinsCounted")])
  sonar_lengths <- sonar_lengths |> within(weights <- Count/nLengths)  ## Assuming 5 min file taken.
  
  ## Add beam width:
  sonar_lengths <- sonar_lengths |> within(beamWidth.cm <- R.m*0.3*pi/180*100)

  ## Will remove single lengths measured in an hour:
  sonar_lengths <- sonar_lengths |> subset(nLengths > dropN)
  
  ## Sonar Lengths _ Adjust for roll angle:
  sonar_lengths$L.cm <- sonar_lengths$L.cm/cos(35*pi/180) ## Hard coded but the roll angle does sometimes change...**********************
  
  self$sonar_lengths <- sonar_lengths
  
  salmon_counts <- sonar_counts |> aggregate(SalmonFlux ~ Date + SonarBin, sum)  
  salmon_counts <- salmon_counts |> reshape(direction = "wide", idvar = "Date", timevar = "SonarBin")
  salmon_counts <- salmon_counts |> within(nearshore <- SalmonFlux.Bin1 + SalmonFlux.Bin2) |>
                                    within(offshore <- SalmonFlux.Bin3) |>
                                    within(count <- nearshore + offshore) 
  salmon_counts <- salmon_counts |> subset(select = c("Date", "count", "nearshore", "offshore"))
  
  ## Sonar counts to total salmon by offshore and nearshore.
  self$salmon_counts <- salmon_counts
}

#' Process Qualark Catch
#'
#' @param self R6 speciesCompModel Object.
#' @param test_fishery_counts Data frame in standard format for Whonnock catch.
#'
#' @details input data frame should have columns start_out, full_out, start_in, full_in along with TRIP_DTT, and net_length. The other 
#' required columns are generally `Sockeye Adult Gilled` and `Sockeye Adult Tangled`.
#'
#' @return No object returned, appends test fishery catch to \code{self$test_fishery_catch[[name]]}
#'
#' @export
process_qualark_catch <- function(self, test_fishery_counts){

  test_fishery_counts <- test_fishery_counts |> subset(!is.na(Date))

  if(inherits(test_fishery_counts$`Start Time`, "POSIXct") & inherits(test_fishery_counts$`End Time`, "POSIXct")){
    test_fishery_counts <- test_fishery_counts |> within(soak_time <- as.numeric(difftime(`End Time`, `Start Time`, "min")))
  }else{
  test_fishery_counts <- test_fishery_counts |> 
                         within(start_out <- process_minutes(`Start Time`)) |>
                         within(full_out <- process_minutes(`End Time`)) |>
                         within(soak_time <- full_out - start_out)
  }

  ## Should be close to 5 minutes:
  if(any(test_fishery_counts$soak_time > 10 | test_fishery_counts$soak_time == 0)) cat("[Warning]   Test fishery soak times are not as expected. Double check via $data_list$test_fishery_counts.")

  ## Soak Time and Effort:
  ## *** pvdb Don't like correction factor. Discounting the impact of in-time makes the CPUE relationship linear...
  test_fishery_counts <- test_fishery_counts |> 
                         within(effort <- soak_time*100/1000)

  test_fishery_counts <- test_fishery_counts |> 
                            within(adultchinook <- `Adult Chinook Retained Ad P` + `Adult Chinook Retained Ad A` + `Adult Chinook Retained Ad Unk` + `Adult Chinook Released Ad P` + `Adult Chinook Released Ad A` + `Adult Chinook Released Ad Unk`) |>
                            within(sockeye <- `Adult Sockeye Retained` + `Jack Sockeye Retained` + `Adult Sockeye Released` + `Jack Sockeye Released`) |>
                            within(coho <- `Coho Retained Ad P` + `Coho Retained Ad A` + `Coho Retained Ad Unk` + `Coho Released Ad P` + `Coho Released Ad A` + `Coho Released Ad Unk`) |>
                            within(pink <- `Pink Retained` + `Pink Released`) |>
                            within(chum <- `Chum Retained` + `Chum Released`) |>
                            within(jackchinook <- `Jack Chinook Retained Ad P` + `Jack Chinook Retained Ad A` + `Jack Chinook Retained Ad Unk` + `Jack Chinook Released Ad P` + `Jack Chinook Released Ad A` + `Jack Chinook Released Ad Unk`) |>                            
                            within(chinook <- jackchinook + adultchinook)
  
  test_fishery_counts <- test_fishery_counts |> within(session <- ifelse(`Drift Number` %in% 1:3, "AM", "PM"))
  
  out <- test_fishery_counts |> aggregate(cbind(adultchinook, sockeye, coho, pink, chum, jackchinook, chinook, effort) ~ Date + session, sum)

  out$net_type <- "vmn"
  
  self$test_fishery_catch[["qualark"]] <- out
}
