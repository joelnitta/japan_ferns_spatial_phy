# make.R

# Master script for running Pleurorosiopsis analyses.
# This project uses drake to manage workflows.
# For more information about drake, see
# https://ropensci.github.io/drake/

# Setup ----

# For reproducibility
set.seed(9130)

# Set working directory
setwd(here::here())

# Load packages
source("code/packages.R")

# Update drake settings
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# Load functions and plans  ----
source("code/functions.R")
source("code/plan.R")

# Run analyses ----

# If not running in parallel, just call make(plan)
make(plan) 

# OR, uncomment out the below lines to run in parallel.
# 
# future::plan(future::multiprocess)
# 
# tictoc::tic()
# make(
#   plan,
#   parallelism = "future",
#   jobs = 4,
#   prework = list(
#     quote(conflict_prefer("map", "purrr")),
#     quote(conflict_prefer("select", "dplyr")),
#     quote(conflict_prefer("filter", "dplyr")),
#     quote(conflict_prefer("gather", "tidyr"))
#   ))
# tictoc::toc()
