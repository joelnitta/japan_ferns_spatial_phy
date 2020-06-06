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
  # Resolve conflicts within each job
  prework = list(
    quote(map <- purrr::map),
    quote(select <- dplyr::select),
    quote(filter <- dplyr::filter),
    quote(gather <- tidyr::gather)
  )
)

# - serial
# drake_config(plan, verbose = 1, cache = ja_fern_cache, seed = 0)
