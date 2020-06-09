FROM rocker/rstudio:4.0.0

ARG DEBIAN_FRONTEND=noninteractive

############################
### Install APT packages ###
############################

# libmagick for animation->magick->phytools
# the rest are same dependencies as rocker/geospatial
# https://hub.docker.com/r/rocker/geospatial/dockerfile
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    lbzip2 \
    libfftw3-dev \
    libgeos-dev \
    libgdal-dev \
    libgsl-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libhdf4-alt-dev \
    libhdf5-dev \
    libjq-dev \
    libpq-dev \
    libproj-dev \
    libprotobuf-dev \
    libnetcdf-dev \
    libsqlite3-dev \
    libssl-dev \
    libudunits2-dev \
    netcdf-bin \
    postgis \
    protobuf-compiler \
    sqlite3 \
    tk-dev \
    unixodbc-dev \
    libgdal-dev \
    libmagick++-dev \
    libzmq3-dev \
    libgmp3-dev \
    libpng-dev

####################################
### Install R packages with renv ###
####################################

# Create directory for renv project library
RUN mkdir renv

# Modify Rprofile.site so renv uses /renv for project library
RUN echo 'Sys.setenv(RENV_PATHS_LIBRARY = "/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Initialize a 'dummy' project and restore the renv library.
# Since the library path is specified as above, the library will be restored to /renv
RUN mkdir tmp/project

COPY ./renv.lock tmp/project

WORKDIR tmp/project

# Don't use cache (the symlinks won't work from Rstudio server)
RUN Rscript -e 'install.packages("renv"); renv::consent(provided = TRUE); renv::settings$use.cache(FALSE); renv::init(bare = TRUE); renv::restore()'

WORKDIR /home/rstudio/
