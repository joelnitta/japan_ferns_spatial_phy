install.packages("BiocManager", repos = "https://cran.rstudio.com/")
my_repos <- BiocManager::repositories()
my_repos["CRAN"] <- "https://cran.rstudio.com/"
options(repos = my_repos)

install.packages("remotes", repos = "https://cran.rstudio.com/")
remotes::install_github("rstudio/renv")
renv::consent(provided = TRUE)
renv::restore(library = "renv")