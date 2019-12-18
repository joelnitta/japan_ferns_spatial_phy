#! /bin/bash
apt-get update
# Install same dependencies as rocker/geospatial
# https://hub.docker.com/r/rocker/geospatial/dockerfile
apt-get install -y --no-install-recommends lbzip2 \
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
    libudunits2-dev \
    libgdal-dev \
    libmagick++-dev \
    libzmq3-dev
Rscript install_packages.R
