FROM rocker/r-ver:3.6.1

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
    libgsl0-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libhdf4-alt-dev \
    libhdf5-dev \
    libjq-dev \
    liblwgeom-dev \
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
    libzmq3-dev

####################################
### Install R packages with renv ###
####################################

COPY ./renv.lock ./

COPY ./renv_restore.R ./

RUN mkdir renv

RUN Rscript renv_restore.R

# Modify Rprofile.site so R loads renv library by default
RUN echo '.libPaths("/renv")' >> /usr/local/lib/R/etc/Rprofile.site

WORKDIR /home/
