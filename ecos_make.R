#####
#
# ecos_make.R 
#
# Script for submitting jobs to run ecostructure analyses. This is done
# separately from the rest of make.R so that these resource-heavy jobs can be 
# run in parallel on the cluster.
#
# Make sure make.R has been run first so that all
# input data are up-to-date.
#
# Then, this script should be run in the background
# after loading the R module and other dependencies:
#
# module load /home/nittaj/modulefiles/miniconda
# source activate pacferns
# nohup Rscript ecos_make.R > make.log 2>&1 &
# 
# Note that running ecos_make.R this way will may result in it getting
# killed for using too much memory. The jobs will still get
# submitted and eventually finish though. The one downside is that drake doesn't
# know about this, and still think the targets haven't been built yet.
#
#####

# Set working directory.
setwd(here::here())

# Load packages and resolve namespace conflicts
library(future.batchtools)
source("code/packages.R")
source("code/resolve_conflicts.R")

# Update drake settings.
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# Load functions and plans.
source("code/functions.R")
source("code/ecos_plan.R")

# Setup cache.
ecos_cache <- new_cache("ecos_cache")

# For reproducibility
set.seed(9130)

# Set up parallel backend
future::plan(batchtools_sge, template = "sge_batchtools_hydra.tmpl")

# Make plan
# Submit one job per ecostructure analysis.
# In prework, resolve conflicts and load data objects needed for each job.

make(
  ecos_plan,
  parallelism = "future", 
  jobs = nrow(ecos_plan), 
  cache = ecos_cache,
  prework = list(
    quote(source("code/resolve_conflicts.R")),
    quote(set.seed(9130))
  )
)
