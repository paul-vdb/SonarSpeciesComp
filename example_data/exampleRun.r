## EM Algorithm Solver:
library(RTMB)
library(mixtools)
library(R6)
source("../R/utils.r")
source("../R/distributions.r")
source("../R/processInputs.r")
source("../R/EMalgorithm.r")
source("../R/simulate.r")

## Delete me later:
##------------------------------------------------------------------------------------
library(tidyverse)
library(readxl)
setwd("C:/Users/vandambatesp/Documents/GitHub/SonarSpeciesComp/R")
sonarLengths <- read.csv("../example_data/2023MergedLengths_new35cmonly_Oct2024.csv")
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

speciesComp <- speciesCompSummary$new(species = c("largeresident", "jackchinook", "sockeye", "adultchinook"), site = "Mission", estDate = "2023-08-05")
speciesComp$processData(sonarCounts, sonarLengths, testFisheryCounts) ## Process all data:
speciesComp$setDate("2023-08-10", 3)
speciesComp$controlOptimization(verbose = TRUE)

speciesComp$setDate("2023-09-05")
speciesComp$setSpecies(species = c("largeresident", "jackchinook", "pink", "sockeye", "adultchinook"))
speciesComp$setSpeciesLengths(testFisheryLengths = testFisheryLengths, ndays = 10)
speciesComp$setModelProportions(formula = list("pink" = ~ -1 + poly(HourOrder,2) + day:SonarBank:SonarAim + SonarBin, "sockeye" = ~ -1 + day:SonarBank:SonarAim + SonarBin))
speciesComp$setLengthAdjustment(formula = ~ -1 + SonarBin, adjustLengths = TRUE)
speciesComp$setModelParameters(fixedParameters = c("mu", "sigma", "muChinook", "sigmaChinook"),
  parameterValues = list(beta = c(0, 4, 7)), testFisheryWeights = 200)
speciesComp$parameters$mu["pink"] <- 51
speciesComp$fitModel()

tru.params <- speciesComp$estimatedParameters
tru.p <- speciesComp$estimatedDailyProportions

for( i in 1:100 ){
  speciesComp$simulate()
  speciesComp$fitModel(simulatedData = TRUE)
  sim.params <- speciesComp$estimatedParameters
  sim.p <- speciesComp$estimatedDailyProportions
}

plot(sim.p$pink - tru.p$pink)
abline(h = 0, col = 'red')

propdf <- speciesComp$estimatedHourlyProportions %>% pivot_longer(c(largeresident, pink, sockeye, jackchinook, adultchinook), names_to = 'species', values_to = 'estimate')

ggplot(data = propdf %>% filter(day == 2), aes(x = HourOrder, y = estimate*SalmonCount, colour = species, shape = SonarAim)) + 
  geom_point() + geom_line() + 
  facet_wrap(~ SonarBin + SonarBank) + theme_bw()

ggplot(data = propdf %>% filter(day == 2), aes(x = HourOrder, y = estimate, colour = species, shape = SonarAim)) + 
  geom_point() + geom_line() + 
  facet_wrap(~ SonarBin + SonarBank) + theme_bw()

speciesComp$estimatedParameters$beta
exp(speciesComp$estimatedParameters$logsigma0)
exp(speciesComp$estimatedParameters$logq)
speciesComp$parsFixed
speciesComp$adjustLengths

speciesComp$analysisData$testFisheryCounts

speciesComp$parsFixed
speciesComp$parsInit
speciesComp$includeTestFishery
speciesComp$adjustLengths


L <- speciesComp$analysisData$lengthData$L.cm.adj - as.matrix(speciesComp$analysisData$Xlength) %*% speciesComp$parsFixed$beta
plot(density(L))
lines(density(speciesComp$analysisData$lengthData$L.cm.adj), col = 'red')
speciesComp$parameters$sigma["pink"]


speciesComp$analysisData$lengthData %>% filter(weights > 200) %>% ggplot(aes(x = L.cm)) + geom_density() + facet_wrap(~SonarBin)

table(speciesComp$analysisData$lengthData$Hour)
sonarLengths %>% 
  mutate(MissionDate = as.Date(MissionDate, format = "%m/%d/%Y")) %>% 
  filter(MissionDate == "2023-09-03", SonarCode == "A1", SonarBin == "Bin1", Hour == 11)



## Model from scratch:
speciesComp <- speciesCompSummary$new(
            species = c("largeresident", "jackchinook", "sockeye", "adultchinook"), 
            site = "Mission", estDate = "2023-08-05")
speciesComp$processData(sonarCounts, sonarLengths, testFisheryCounts) ## Process all data:
speciesComp$setDate(estDate = "2023-08-10", ndays = 3)
speciesComp$controlOptimization(verbose = TRUE)

speciesComp$setDate("2023-08-01")
speciesComp$setSpecies(species = c("largeresident", "jackchinook", "sockeye", "adultchinook"))
speciesComp$setSpeciesLengths(testFisheryLengths = testFisheryLengths, ndays = 10)
speciesComp$setModelProportions(formula = list("sockeye" = ~ -1 + day:SonarBank:SonarAim + SonarBank:SonarBin, 
                                        "largeresident" = ~ -1 + day:SonarBank:SonarAim + SonarBank:SonarBin)
                                )
# speciesComp$setLengthAdjustment(formula = ~ -1 + SonarBin, adjustLengths = TRUE)
speciesComp$setLengthAdjustment(formula = ~ R.m, adjustLengths = TRUE)

speciesComp$setModelParameters(fixedParameters = c("mu", "sigma", "muChinook", "sigmaChinook", "sigma0"),
  parameterValues = list(beta = c(0, 0), sigma0 = 5), testFisheryWeights = 200)

speciesComp$setPriors(priors = list(beta = function(x){sum(dnorm(x, 0, 0.1, log = TRUE))}, sigma0 =  function(x){sum(dgamma(x, 1, 0.2, log = TRUE))}))
speciesComp$fitModel()

tru.params <- speciesComp$estimatedParameters
tru.p <- speciesComp$estimatedDailyProportions

speciesComp$simulate()
Lnew <- speciesComp$simData$lengthData$L.cm
L <- speciesComp$analysisData$lengthData$L.cm

# speciesComp$plotProps()

ggplot(data = speciesComp$analysisData$lengthData, aes(x = L.cm)) + geom_density(aes(colour = SonarBin, linetype = SonarBank)) + 
  facet_wrap(~day) + 
  theme_bw()
dev.new()
ggplot(data = speciesComp$simData$lengthData, aes(x = L.cm)) + geom_density(aes(colour = SonarBin, linetype = SonarBank)) + 
  facet_wrap(~day) + 
  theme_bw()

speciesComp$analysisData$lengthData %>% group_by(SonarBank, SonarBin, SonarAim, day) %>% reframe(modes = findModes( L.cm ) ) %>%
  ggplot(aes(x = SonarBin, y = modes, colour = SonarBin, shape = SonarBank)) + 
    facet_wrap(~day) + 
    geom_point() + 
    theme_bw()

# speciesComp$setModelParameters(fixedParameters = c("mu", "sigma", "muChinook", "sigmaChinook", "beta"),
  # parameterValues = list(beta = tru.params$beta, sigma0 = 3.5), testFisheryWeights = 200)
speciesComp$fitModel()
ests <- NULL
for( i in 1:100 ){
  speciesComp$simulate()
  speciesComp$fitModel(simulatedData = TRUE)
  # sim.params <- speciesComp$estimatedParameters
  sim.p <- speciesComp$estimatedDailyProportions
  sim.p$iter <- i
  ests <- rbind(ests, sim.p)
}

plot(speciesComp$simData$lengthData$L.cm.adj)
plot(speciesComp$analysisData$lengthData$L.cm.adj)

L.adj <- speciesComp$analysisData$lengthData$L.cm.adj - speciesComp$analysisData$Xlength %*% speciesComp$estimatedParameters$beta

plot(density(L.adj), col = 'red')
lines(density(speciesComp$analysisData$lengthData$L.cm.adj), col = 'blue')

ggplot(data = ests %>% filter(SonarBin == "Bin1", day == 3), aes(x = SonarBank, y = sockeye)) + 
  geom_boxplot() + 
  geom_hline(data = tru.p %>% filter(SonarBin == "Bin1", day == 3), aes(yintercept = sockeye), col = 'red') +
  facet_wrap(~SonarBank + SonarAim) +
  theme_bw()


L.adj <- speciesComp$analysisData$lengthData$L.cm.adj - speciesComp$analysisData$Xlength %*% speciesComp$estimatedParameters$beta
mu <- c(speciesComp$estimatedParameters$mu, speciesComp$estimatedParameters$muChin)
sigma <-  sqrt(c(exp(2*speciesComp$estimatedParameters$logsigma) + exp(2*speciesComp$estimatedParameters$logsigma0), exp(2*speciesComp$estimatedParameters$logsigmaChin) + exp(2*speciesComp$estimatedParameters$logsigma0)))
pn <- t(sapply(1:length(L.adj), FUN = function(x){dnorm(L.adj[x], mu, sigma)}))
mn <- t(apply(pn, 1, FUN = function(x){x/max(x)}))
plot(mn[,1])

plot(L.adj)



fn <- MakeTape(\(x){speciesComp$priorDists$dsigma(x[1]) + speciesComp$priorDists$dsigmaChin(x[2])}, c(0,0))
fn(c(0,0))