
library(targets)
library(tarchetypes)
library(future)

# targets options ----
tar_option_set(
  packages = c("tibble",
               "readr",
               "tidyr",
               "dplyr",
               "stringi",
               "lubridate",
               "survival",
               "glmnet",
               "purrr",
               "Boruta")
)

# load functions ----
tar_source("R")

# Data folders ----
rawdatadir <- Sys.getenv("RAWDATA_DIR")
rawdatadir2 <- Sys.getenv("RAWDATA_DIR2")

# Pipeline ----
tar_source("targets")
list(
  targets_prep,                            # import and prepare data
  targets_individual_proteins_analysis,    # individual protein associations with MI
  targets_predictors,                      # define predictors
  targets_predictions,                     # predictions and validation
  targets_whole_sample,                    # prediction using whole sample
  targets_ukb                              # external validation
)
