
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
# use_package("dplyr")
# use_package("ggplot2")

use_build_ignore("package_setup.R")

usethis::use_tidy_description()

document()
