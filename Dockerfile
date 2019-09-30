FROM rocker/geospatial:3.5.1

ARG DEBIAN_FRONTEND=noninteractive

####################
### APT packages ###
####################

# libudunits2-dev, libgdal-dev for CoordinateCleaner
# libmagick for animation->magick->phytools
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  libmagick++-dev \
  libudunits2-dev \
  libgdal-dev

#######################################
### Install R packages with packrat ###
#######################################

# First install dependencies of methClust and CountClust that for
# some reason don't get installed automatically.
RUN install2.r -e slam \
  SQUAREM \
  boot \
&& Rscript -e 'devtools::install_github("kkdey/CountClust")'

COPY ./packrat/packrat.lock packrat/

RUN Rscript -e 'install.packages("packrat", repos = "https://cran.rstudio.com/"); packrat::restore()'

# Modify Rprofile.site so R loads packrat library by default
RUN echo '.libPaths("/packrat/lib/x86_64-pc-linux-gnu/3.5.1")' >> /usr/local/lib/R/etc/Rprofile.site

#############################
### Other custom software ###
#############################

ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

### gnparser ###
ENV APP_NAME=gnparser
ENV VERSION=0.7.5
ENV DEST=$APPS_HOME/$APP_NAME/$VERSION
RUN wget https://www.dropbox.com/s/7jcrjj0o39vuh3x/$APP_NAME-v$VERSION-linux.tar.gz?dl=1 \
  && tar xf $APP_NAME-v$VERSION-linux.tar.gz?dl=1 \
  && rm $APP_NAME-v$VERSION-linux.tar.gz?dl=1 \
  && mv "$APP_NAME" /usr/local/bin/

WORKDIR /home/rstudio/
