# dev packages
# Not used for running code,
# but included here so renv knows to keep them around
library(remotes)
library(visNetwork)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(miniUI)
library(lwgeom)
library(tflow)
library(IRanges) # for fuzzyjoin::interval_inner_join()
library(canaper) # FIXME: remove this once canaper is used in _targets.R

# considering these...
library(ggtext)
library(gghighlight)

# needed for patchwork to plot ggplot and base together
library(gridGraphics)