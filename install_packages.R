# This script is used to intall packages for tracking with packrat.
# For more info on installing R packages to docker images with
# packrat, see https://www.joelnitta.com/post/docker-and-packrat/

# This should be run in a conda environment with R 3.5.1

### Initialize packrat ###

# Don't let packrat try to find
# packages to install itself.

install.packages("packrat", repos = "https://cran.rstudio.com/")
packrat::init(
  infer.dependencies = FALSE,
  enter = TRUE,
  restart = FALSE)

### Setup repositories ###

# Install packages that install packages.
install.packages("BiocManager", repos = "https://cran.rstudio.com/")
install.packages("remotes", repos = "https://cran.rstudio.com/")

# Specify repositories so they get included in
# packrat.lock file.
my_repos <- BiocManager::repositories()
my_repos["CRAN"] <- "https://cran.rstudio.com/"
options(repos = my_repos)

### Install CRAN packages ###
cran_packages <- c(
  "ape",
  "assertr",
  "assertthat",
  "conflicted",
  "clustermq", # For running in parallel on the cluster
  "drake",
  "gghighlight",
  "janitor",
  "maps",
  "picante",
  "readxl",
  "scales",
  "scico",
  "tidyverse",
  "viridis",
  "boot", # ecostructure dependency that doesn't get detected automatically
  "slam", # ecostructure dep
  "SQUAREM" # ecostructure dep
  )

install.packages(cran_packages)

### Install bioconductor packages ###
bioc_packages <- c(
  "Biobase" # ecostructure dep
)

BiocManager::install(bioc_packages)

### Install github packages ###
github_packages <- c(
  "joelnitta/jntools",
  "kkdey/methClust", # ecostructure dep
  "kkdey/CountClust", # ecostructure dep
  "kkdey/ecostructure",
  "TaddyLab/maptpx" # ecostructure dep
)

remotes::install_github(github_packages)

### Take snapshot ###

packrat::snapshot(
  snapshot.sources = FALSE,
  ignore.stale = TRUE,
  infer.dependencies = FALSE)
