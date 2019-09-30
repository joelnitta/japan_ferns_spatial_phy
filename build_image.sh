#! /bin/bash
# Take a snapshot then build the image
docker run --rm -e DISABLE_AUTH=true -v /Users/joelnitta/Documents/japan_ferns_biogeo:/home/rstudio/project rocker/geospatial:3.5.1 bash /home/rstudio/project/install_packages.sh
docker build . -t joelnitta/japan_ferns_biogeo:3.5.1
rm .Rprofile
