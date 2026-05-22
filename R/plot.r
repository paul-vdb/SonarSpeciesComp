#' Weighted histogram for plotting length data
#' 
#' @param x Vector of values to pass to histogram.
#' @param wgts Weights of each value.
#' @param range Range that the histogram plot will cover.
#' @param delta The step size to bin by.
#' @param ylim ylimit to be used.
#' @param xlim xlimit to be used.
#' @param xlab xlabel to plot with.
#' @param plot logical to either plot the histogram or not.
#' @param ... Additional arguments for `barplot` function call.
#' @details Taken from `plotrix` package for the weighted histogram, reworked for specific purposes.
#'
#' @return list of range to plot under (`$range`), the binned values (`$lengths`), and histplot `$density`.
#'
#' @export
histplot <- function(x, wgts, range = c(35, 120), delta = 2, ylim = NA, xlim = NA, xlab = NA, plot = TRUE, ...){
    if (missing(x)) 
        stop("Must pass a value to x.")
    if (missing(wgts)) 
        wgts <- rep(1, length(x))
    extra <- delta - diff(range) %% delta
    if(extra != 0) range <- range + c(-extra/2, extra/2)
    n <- diff(range)/delta
    breaks <- seq(range[1], range[2], by = delta)
    
    counts <- NULL
    for ( bin in 1:(length(breaks)-1) ){
      counts <- c(counts, sum(wgts[x >= breaks[bin] & x < breaks[bin + 1]]))
    }
    density <- counts/sum(counts)
    width <- rep(delta, length(density))
    density <- density/width
    
    if (is.na(ylim)) ylim <- c(0, 1.1 * max(density, na.rm = TRUE))
    if (is.na(xlab)) xlab = "Fish Length (cm)"
    mids <- barplot(density, width =  width, col = "grey", space = 0, 
        ylim = ylim, ylab = "Density", xlab = xlab, ...)
    tickpos <- c(mids - width/2, mids[length(mids)] + width[length(width)]/2)
    axis(1, at = tickpos, labels = signif(breaks, 3))

    return(list(range = range, lengths = breaks[-length(breaks)] + delta/2, density = density))
}

#' Plot weighted histogram along with estimated proportions.
#' 
#' @param self R6 speciesCompModel object.
#' @param day Day to plot (default = 1). Day 1 is most recent day.
#' @param include_proportion_labels Include labels in plot (default = FALSE).
#' @param ... Additional options to pass to plot function.
#'
#' @details Taken from `plotrix` package for the weighted histogram, reworked for specific purposes. This function
#' is called directly from `$plot` within the R6 object.
#'
#' @return list of range to plot under (`$range`), the binned values (`$lengths`), and histplot `$density`.
#'
#' @export
plot_mix <- function(self, day = 1, ...){
  if(!any(self$est_date %in% self$params_estimated$N_daily$Date)) 
    stop("Must run 'fitModel' before trying to predict species composition for this date.")

  mu <- self$params_estimated$mu_adjusted
  sigma <- self$params_estimated$sigma
  sigma0 <- self$params_estimated$sigma0
  beta <- self$params_estimated$beta

  dd <- self$est_date - day + 1
  stratum <- self$data_list$pred_df |> subset(Date == dd, select = stratum)
  stratum <- unique(stratum$stratum)
  if( !self$fit_info$adjust_lengths ) 
    self$data_list$length_data$L.cm.modadj <- self$data_list$length_data$L.cm.adj
  else 
    self$data_list$length_data$L.cm.modadj <- self$data_list$length_data$L.cm.adj - self$data_list$X_length %*% beta

  length_df <- self$data_list$length_data |> subset(Date == dd)
  xx <- seq(10, 110, 1)
  species <- self$species_info$species
  p <- self$params_estimated$p_daily |> subset(Date == dd) |> subset(select=-Date) |> as.numeric()
  fx <- do.call("rbind", lapply(1:length(species), FUN = function(i){data.frame(x = xx, f = p[i]*dnorm(xx, mu[i], sqrt(sigma[i]^2 + sigma0^2)), species = species[i])}))
  fx <-  fx |> within(species <- factor(species, levels = self$species_info$species))
  fx_all <- fx |> aggregate(f~x, sum)

  if (require("ggplot2", quietly = TRUE)) {
    plot_h <- ggplot(length_df) + 
      geom_histogram(aes(x = L.cm.modadj, y = ..density.., weight = weights), binwidth = 2, alpha = 0.5, colour = "black") +
      geom_line(data = fx, aes(x = x, y = f, colour = species), linewidth = 1) +
      theme_bw() +
      scale_colour_manual("Species", labels = speciesLabels(species), values = speciesColours(species)) +
      xlab("Fish Length (cm)") + ylab("Density") + 
      geom_line(data = fx_all, aes(x = x, y = f), colour = "black", linetype = 2, linewidth = 1) +
      ggtitle(paste0("Date: ", dd))
    suppressWarnings(print(plot_h))
    return(invisible(NULL))
  } else {

    p <- self$params_estimated$p_daily |> subset(Date == dd)
    N <- self$data_list$total_salmon |> subset(Date == dd)
    hplot <- histplot(length_df$L.cm.modadj, length_df$weights, xlab = "Fish Length (cm)", range = c(10, 110), main = paste(dd, "N =", floor(N$count)), ...)
  
    range <- hplot$range
    for( i in seq_along(species)){
      sppi <- species[i]
      fxi <-  fx |> subset(species == sppi)
      lines(fxi$x - range[1], fxi$f, col = speciesColours(sppi), lty = 1, lwd = 2) ## Have to adjust to scale of bar plot
    }
    lines(fx_all$x-range[1], fx_all$f, col = "black", lwd = 2, lty = 2)
    legend("topright", legend = speciesLabels(self$species_info$species), col = speciesColours(self$species_info$species), 
            lty = rep(1, length(self$species_info$species)), ncol = 2)
  }
  return(invisible(NULL)) 
}

#' Plot Pearson residuals of test fishery model (CPUE - Nq)
#' 
#' @param self R6 speciesCompModel object.
#'
#' @details Plot the Pearson residuals for the test fishery catch, displayed as CPUE - Nq. x-axis labels are species - test fishery - net type.
#'
#' @export
plot_test_fishery <- function(self){
  if(!any(self$est_date %in% self$params_estimated$N_daily$Date)) 
    stop("Must run 'fitModel' before trying to predict species composition for this date.")

  test_catch <- self$data_list$test_fishery_catch
  N <- as.numeric(self$params_estimated$N_daily[cbind(test_catch$day,test_catch$N_index+1)])
  qinv <- self$params_estimated$qinv[test_catch$q_index]
  CPUE <- test_catch$catch/test_catch$effort
  E_CPUE <- N/qinv
  
  test_catch <- test_catch |> within(par <- factor(paste(test_catch$species, test_catch$fishery, test_catch$net_type, sep = "_")))
  test_catch$diff <- CPUE - E_CPUE

  if (require("ggplot2", quietly = TRUE)) {
    plot_h <- ggplot(data = test_catch, aes(x = factor(par), y = diff)) + 
      geom_point() + 
      theme_bw() + 
      ggtitle( paste0("Test Fishery: ", self$est_date) ) +
      xlab("") + ylab("CPUE - Nq") + 
      geom_hline(yintercept = 0, col = 'red', linetype = 2)
      suppressWarnings(print(plot_h))
  }else {
    plot(as.numeric(test_catch$par), CPUE - E_CPUE, xlab = "", ylab = "CPUE - Nq", pch = 16, xaxt = "n", main = paste0("Test Fishery: ", self$est_date))
    axis(1, at = sort(unique(as.numeric(test_catch$par))), labels = gsub("_", "\n ", levels(test_catch$par)), tick = TRUE, padj = 0.5)
    abline(h = 0, col = 'red', lty = 2)
    for( i in 1:length(self$params_estimated$qinv)){
      mtext(text=paste0('Expansion - ', names(self$params_estimated$qinv[i]), ": ",  round(self$params_estimated$qinv[i], 2)), 
        side = 3, adj = 1-0.08, line = i-3, cex = 1)  
    }
  }
  return(invisible(NULL))
}

#' Plot Beam Spreading Trend
#' 
#' @param self R6 speciesCompModel object.
#'
#' @details Plot the change in expected species length against sonar range (m).
#'
#' @export
plot_beam_spreading <- function(self){
  if(!any(self$est_date %in% self$params_estimated$N_daily$Date)) 
    stop("Must run 'fitModel' before trying to predict species composition for this date.")
  
  beta <- self$params_estimated$beta
  L.cm.modadj <- self$data_list$length_data$L.cm.adj - self$data_list$X_length %*% beta
  Xnew <- data.frame(R.m = 1:30, SonarBin = rep(c("Bin1", "Bin2", "Bin3"), each = 10))
  Xnew <- Xnew |> within(beamWidth.cm <- R.m*0.3*pi/180*100)
  Xnew <- Xnew |> subset(SonarBin %in% unique(self$data_list$length_data$SonarBin))
  predmat <- model.matrix(self$data_list$length_adjust_formula, data = Xnew)
  Xnew$adjust <- (predmat %*% beta)[,1]
  
  ## Now adjust species lengths:
  mu <- self$params_estimated$mu_adjusted
  spp <- self$species_info$species

  if (require("ggplot2", quietly = TRUE)) {
    X2 <- NULL
    for( i in seq_along(spp) ) 
      X2 <- rbind(X2, Xnew |> within(species <- spp[i]) |> within(mu <- mu[i] + adjust)) 
      X2 <- X2 |> within(species <- factor(species, levels = spp))

    plot_h <- ggplot(data = X2, aes(x = R.m, y = mu, colour = species)) + 
      geom_line(linewidth = 1) +
      ylab('Expected Length (cm)') + 
      xlab('Range (m)') + 
      ggtitle(paste0("Beam Spreading: ", self$est_date)) + 
      theme_bw() +
      geom_hline(data = data.frame(mu = as.numeric(mu), species = spp), aes(colour = species, yintercept = mu), linetype = 2) +
      scale_colour_manual("Species", labels = speciesLabels(spp), values = speciesColours(spp))
      suppressWarnings(print(plot_h))
  }else{
    plot(0, ylim = c(0,100), xlim = c(0, max(Xnew$R.m)), pch = '', ylab = 'Expected Length (cm)', xlab = 'Range (m)', main = paste0("Beam Spreading: ", self$est_date))
    for( i in seq_along(spp) ){
      lines(Xnew$R.m, mu[i] + Xnew$adjust, col = speciesColours(spp[i]))
      abline(h = mu[i], col = speciesColours(spp[i]), lty = 2)
    }
    legend("topleft", legend = speciesLabels(spp), col = speciesColours(spp), 
            lty = rep(1, length(spp)), ncol = 2)
  }
  return(invisible(NULL))
}

#' Species labels
#'
#' Labels to plot species names.
#'
#' @param species Vector of species names
#'
#' @return Vector of names for display.
#'
#' @export
speciesLabels <- function(species){
  match <- c("smallresident" = "Small Resident", "largeresident" = "Large Resident", "jackchinook" = "Chinook Jack", 
             "pink" = "Pink", "sockeye" = "Sockeye", "coho" = "Coho", "chum" = "Chum", "adultchinook" = "Chinook Adult", "smalladultchinook" = "Chinook Small Adult", 
             "largeadultchinook" = "Chinook Large Adult")
  as.character(match[species])
}

#' Species colours
#'
#' @param species Vector of species names
#'
#' @return Vector of colours matching the species names.
#'
#' @export
speciesColours <- function(species){
  cols <- c("smallresident" = "#FFB100", "largeresident" = "#656837",  "jackchinook" = "#A33CC7", "pink" = "#FF8DA1", 
    "sockeye" = "#CD0000", "coho" = "#228833", "chum" = "#B8604A", "adultchinook" = "#27408B", "smalladultchinook" = "#27408B", "largeadultchinook" = "#27408B")
  cols[species]
}
