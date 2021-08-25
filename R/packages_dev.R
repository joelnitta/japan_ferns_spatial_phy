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
library(ROI.plugin.glpk)

# considering these...
library(ggtext)
library(gghighlight)

# needed for patchwork to plot ggplot and base together
library(gridGraphics)

# temporarily include here while updating functions.R
# can remove once rgnparser::gn_parse_tidy() is included in a commit to functions.R
library(rgnparser)