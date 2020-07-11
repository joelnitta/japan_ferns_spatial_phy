# dev packages
# Not used for running code,
# but included here so renv knows to keep them around
library(remotes)
library(visNetwork)
library(rnaturalearthdata)
# don't track rnaturalearthhires with renv, it gives this error upon attempting snapshot:
# "The following package(s) were installed from an unknown source...
# renv may be unable to restore these packages in the future."
# library(rnaturalearthhires)

# considering these...
library(ggtext)
library(gghighlight)