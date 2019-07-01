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
source("code/test_swap_plan.R")

# Load cache
test_swap_cache <- new_cache("test_swap")

# Run analyses ----

# If not running in parallel, just call make(plan)
make(test_swap_plan, cache = test_swap_cache) 

