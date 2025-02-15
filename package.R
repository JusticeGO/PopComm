library(usethis)
library(devtools)
library(roxygen2)


has_devel()
available.packages()

dir.create("f:/R_package/PopComm")
setwd("f:/R_package/PopComm")


usethis::create_package("PopComm")
