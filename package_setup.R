
library(devtools)
library(usethis)
#
# options(
#   usethis.description = list(
#     "Authors@R" = utils::person(
#       "Paul", "van Dam-Bates",
#       email = "paul.vandambates@gmail.com",
#       role = c("aut", "cre")
#     ),
#     License = "MIT + file LICENSE"
#   )
# )
# use_description(fields = list(), check_name = TRUE, roxygen = TRUE)

use_package("RTMB")
use_package("R6")
# use_package("dplyr")
# use_package("ggplot2")

use_build_ignore(c("package_setup.R", "example_data"))
use_vignette("SpeciesCompMethods.qmd", "Technical Model Details")

load("data/Mission2025.Rdata")
use_data(mission_2022)
use_data(mission_2025)

usethis::use_tidy_description()
devtools::document()
devtools::check()
