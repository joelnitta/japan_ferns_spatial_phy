# Load packages, functions, and plan
source("R/packages.R")
source("R/functions.R")
source("R/plan.R")

# Setup cache
ja_fern_cache <- new_cache("ja_fern_cache")
options(rstudio_drake_cache = ja_fern_cache)

# Specify parallel back-end
options(clustermq.scheduler = "multicore")

# Specify non-global environment
# to get around captioner modifying global env
# (cf https://github.com/ropensci/drake/issues/749)
envir <- new.env(parent = globalenv())

# Configure settings for making plan
# (choose either parallel or serial, comment-out the other)

# - parallel
drake_config(
  plan,
  parallelism = "clustermq",
  jobs = 6, # Change to match number of cores available!
  cache = ja_fern_cache,
  seed = 0,
  envir = envir
)

# - serial
# drake_config(plan, verbose = 1, cache = ja_fern_cache, seed = 0, envir = envir)
