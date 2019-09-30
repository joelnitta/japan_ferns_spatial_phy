# This script writes renv.lock for installing R packages to a docker image.
# It should be run from within the rocker/geospatial:3.5.1 container:
#
# docker run --rm -e DISABLE_AUTH=true -v /Users/joelnitta/Documents/japan_ferns_biogeo:/home/rstudio/project rocker/geospatial:3.5.1 bash /home/rstudio/project/install_packages.sh
#
# Then build the image with
# docker build . -t joelnitta/japan_ferns_biogeo:3.5.1

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
  "gghighlight",
  "ggridges",
  "here",
  "janitor",
  "lwgeom",
  "maps",
  "picante",
  "readxl",
  "scales",
  "scico",
  "sf",
  "tidyverse",
  "usedist",
  "viridis",
  "boot", # ecostructure dependency that doesn't get detected automatically
  "slam", # ecostructure dep
  "SQUAREM" # ecostructure dep
  )

install.packages(cran_packages)

### Install bioconductor packages ###
bioc_packages <- c("Biobase", "BiocGenerics")

BiocManager::install(bioc_packages)

### Install github packages ###
github_packages <- c(
  "joelnitta/jntools",
  "joelnitta/taxastand",
  "thomasp85/patchwork",
  "r-lib/scales", # For degree_format(), which isn't in CRAN scales yet
  "kkdey/methClust", # ecostructure dep
  "kkdey/CountClust", # ecostructure dep
  "TaddyLab/maptpx", # ecostructure dep
  "kkdey/ecostructure"
)

remotes::install_github(github_packages)

### Take snapshot ###

packrat::snapshot(
  snapshot.sources = FALSE,
  ignore.stale = TRUE,
  infer.dependencies = FALSE)
