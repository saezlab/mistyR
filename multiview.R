#Explicit attachment of libraries will be removed in future versions
library(MASS)
library(dplyr)
library(purrr)
library(furrr)
library(readr)
library(stringr)
library(tibble)
library(caret)
library(randomForest)
require(ranger)
library(deldir)
library(distances)
library(digest)
library(rlist)
library(assertthat)


# setup default parallel processing
# the future package is required by furrr
future::plan(future::multiprocess)

# import operators
# magrittr and rlang are imported by dplyr
`%>%` <- magrittr::`%>%`
`!!` <- rlang::`!!`
`:=` <- rlang::`:=`


source("views.R")
source("utils.R")
source("models.R")
source("estimators.R")
