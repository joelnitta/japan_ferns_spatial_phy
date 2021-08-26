FROM rocker/verse:4.0.0

ARG DEBIAN_FRONTEND=noninteractive

############################
### Install APT packages ###
############################

# gcc through libtool for treePL
# libmagick for animation->magick->phytools
# unar for unzipping files with unicode filenames
# the rest are same dependencies as rocker/geospatial
# https://hub.docker.com/r/rocker/geospatial/dockerfile
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    mafft \
    fasttree \
    gcc \
    g++ \
    libnlopt-dev \
    libnlopt0 \
    libcolpack-dev \
    make \
    libomp-dev \
    build-essential \
    autoconf \
    autotools-dev \
    automake \
    libtool \
    libmagick++-dev \
    unar \
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
    libzmq3-dev \
    libgmp3-dev \
    libpng-dev

#############################
### Other custom software ###
#############################

ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

### treePL ###
# patches have been applied since v1.0 but no further version tagged.
# so use the most recent commit hash
ENV VERSION=6e23cdb08c7ec9283ebbee909d783cb0d44e64d5
RUN git clone https://github.com/blackrim/treePL.git \
  && cd $APPS_HOME/treePL \
  && git checkout $VERSION \
  && cd $APPS_HOME/treePL/deps/ \
  && tar xvzf adol-c_git_saved.tar.gz \
  && cd $APPS_HOME/treePL/deps/adol-c/ \
  && ./update_versions.sh \
  && ./configure --with-openmp-flag=-fopenmp --prefix=/usr \
  && make \
  && make install \
  && cd $APPS_HOME/treePL/src \
  && ./configure \
  && make \
  && echo '/usr/lib64' > /etc/ld.so.conf.d/lib64.conf \
  && ldconfig \
  && cp treePL /usr/local/bin

### gnparser ###
WORKDIR $APPS_HOME
ENV APP_NAME=gnparser
ENV VERSION=1.3.3
ENV DEST=$APPS_HOME/$APP_NAME/$VERSION
RUN wget https://github.com/gnames/$APP_NAME/releases/download/v$VERSION/$APP_NAME-v$VERSION-linux.tar.gz \
  && tar xf $APP_NAME-v$VERSION-linux.tar.gz \
  && rm $APP_NAME-v$VERSION-linux.tar.gz \
  && mv "$APP_NAME" /usr/local/bin/
  
### IQ Tree ###
WORKDIR $APPS_HOME
ENV APP_NAME=iqtree
ENV VERSION=1.6.12
ENV DEST=$APPS_HOME/$APP_NAME/$VERSION
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v$VERSION/$APP_NAME-$VERSION-Linux.tar.gz \
  && tar xf $APP_NAME-$VERSION-Linux.tar.gz \
  && rm $APP_NAME-$VERSION-Linux.tar.gz \
  && mv "$APP_NAME-$VERSION-Linux/bin/iqtree" /usr/local/bin/

###########################################
### Install latex package with tiny tex ###
###########################################

# \n\ puts each package name on its own line

RUN printf 'amsmath\n\
atbegshi\n\
atveryend\n\
auxhook\n\
bigintcalc\n\
bitset\n\
booktabs\n\
colortbl\n\
environ\n\
etexcmds\n\
etoolbox\n\
euenc\n\
float\n\
fontspec\n\
geometry\n\
gettitlestring\n\
grffile\n\
hycolor\n\
hyperref\n\
iftex\n\
infwarerr\n\
intcalc\n\
kvdefinekeys\n\
kvoptions\n\
kvsetkeys\n\
latex-amsmath-dev\n\
letltxmacro\n\
ltxcmds\n\
makecell\n\
mdwtools\n\
multirow\n\
pdfescape\n\
pdflscape\n\
pdftexcmds\n\
refcount\n\
rerunfilecheck\n\
setspace\n\
siunitx\n\
stringenc\n\
tabu\n\
threeparttable\n\
threeparttablex\n\
tipa\n\
trimspaces\n\
ulem\n\
unicode-math\n\
uniquecounter\n\
varwidth\n\
wrapfig\n\
xcolor\n\
xunicode\n\
zapfding' >> latex_packages.txt \
  && Rscript -e 'tinytex::tlmgr_update(); tinytex::tlmgr_install(readLines("latex_packages.txt"))'	\
  && rm latex_packages.txt	
	
####################################
### Install R packages with renv ###
####################################

# Create directory for renv project library
RUN mkdir /renv

# Modify Rprofile.site so renv uses /renv for project library, and doesn't use the cache
RUN echo 'Sys.setenv(RENV_PATHS_LIBRARY = "/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Initialize a 'dummy' project and restore the renv library.
# Since the library path is specified as above, the library will be restored to /renv

RUN mkdir -p /tmp/project/renv/local

# Copy needed files: renv.lock and a local package
COPY ./renv.lock /tmp/project

# Install packages to the dummy project using renv::restore()
WORKDIR /tmp/project

# Don't use cache (the symlinks won't work from Rstudio server)
RUN Rscript -e 'devtools::install_github("rstudio/renv@0.13.2-87"); renv::consent(provided = TRUE); renv::settings$use.cache(FALSE); renv::init(bare = TRUE); renv::restore()'

WORKDIR /home/rstudio/

RUN rm -rf /tmp/project
