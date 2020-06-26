# Load packages, functions, and plan
source("R/packages.R")
source("R/functions.R")
source("R/phylo_plan.R")

# Setup cache
phylo_cache <- new_cache("phylo_cache")
options(rstudio_drake_cache = phylo_cache)

# Specify settings for making plan
drake_config(plan, verbose = 1, cache = phylo_cache, seed = 0)
