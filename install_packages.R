# Install packages to a docker image with renv

# First do an initial install to write the Renv lock file
# (change wd on left side of colon as needed):
# docker run --rm -v /Users/joelnitta/repos/japan_ferns_biogeo:/home/japan_ferns_biogeo -w /home/japan_ferns_biogeo rocker/r-ver:3.6.1 bash /home/japan_ferns_biogeo/install_packages.sh

### Initialize renv ###

# Initialize renv, but don't let it try to find packages to install itself.
install.packages("remotes", repos = "https://cran.rstudio.com/")
# Use dev version with most recent bug fixes
remotes::install_github("rstudio/renv")

renv::consent(provided = TRUE)

renv::init(
  bare = TRUE,
  force = TRUE,
  restart = FALSE)

renv::activate()

### Setup repositories ###

# Install packages that install packages.
install.packages("BiocManager", repos = "https://cran.rstudio.com/")
# (Need to do this again because now we're in a fresh renv project)
install.packages("remotes", repos = "https://cran.rstudio.com/")

# Set repos.
my_repos <- BiocManager::repositories()
my_repos["CRAN"] <- "https://cran.rstudio.com/"
options(repos = my_repos)

### Install CRAN packages ###
cran_packages <- c(
  "ape",
  "assertr",
  "assertthat",
  "checkr",
  "conflicted",
  "clustermq", # For running in parallel on the cluster
  "drake",
  "FD",
  "future.batchtools", # For running in parallel on the cluster
  "here",
  "janitor",
  "lwgeom",
  "maps",
  "picante",
  "readxl",
  "tidyverse",
  "usedist"
  )

install.packages(cran_packages)

### Install bioconductor packages ###
bioc_packages <- c("Biobase", "BiocGenerics")

BiocManager::install(bioc_packages, update=FALSE, ask=FALSE)

### Install github packages ###
github_packages <- c(
  "joelnitta/jntools",
  "joelnitta/taxastand",
  "kkdey/CountClust", # ecostructure dep
  "joelnitta/methClust", # ecostructure dep
  "TaddyLab/maptpx", # ecostructure dep
  "joelnitta/ecostructure"
)

remotes::install_github(github_packages)

### Take snapshot ###

renv::snapshot(type = "simple")
