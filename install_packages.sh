#! /bin/bash
cd /home/rstudio/project
apt-get update
apt-get install -y --no-install-recommends libudunits2-dev libgdal-dev libmagick++-dev
Rscript install_packages.R
