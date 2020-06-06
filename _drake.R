# Load packages, functions, and plan
source("R/packages.R")
source("R/functions.R")
source("R/plan.R")

# Setup cache
ja_fern_cache <- new_cache("ja_fern_cache")
options(rstudio_drake_cache = ja_fern_cache)

# Specify multicore backend
options(clustermq.scheduler = "multicore")

# Set up drake configuration
# (choose one to use, comment-out the other)

# - multicore
drake_config(
  plan,
  parallelism = "clustermq",
  jobs = 6, # Change to match number of cores available!
  cache = ja_fern_cache,
  seed = 0,
  prework = list(
    quote(conflict_prefer("map", "purrr")),
    quote(conflict_prefer("select", "dplyr")),
    quote(conflict_prefer("filter", "dplyr")),
    quote(conflict_prefer("gather", "tidyr"))
  )
)

# - serial
# drake_config(plan, verbose = 1, cache = ja_fern_cache, seed = 0)
