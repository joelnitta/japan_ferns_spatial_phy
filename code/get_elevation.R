# Get elevation for grid cells from latitude and longitude.

# Fetching elevation data requires an API key, so just do it once
# and save the results to "data/"

library(tidyverse)
library(rgbif)
library(raster)

# Read in geonames username. See
# rgbif::elevation() for how to set this up.
user <- Sys.getenv("GEONAMES_USER")

# Read in data for 10km grid cells across Japan.
# Original data lack elevation.
all_cells = read_csv("data_raw/2_grid_cells_all.csv") %>%
  rename(longitude = x, latitude = y, secondary_grid_code = id) %>%
  assert(is_uniq, secondary_grid_code)

# Add an elevation column by searching from longitude and latitude.
all_cells_el <- all_cells %>%
  select(decimalLatitude = latitude, decimalLongitude = longitude, everything()) %>%
  as.data.frame() %>%
  # This is where the API key is needed
  elevation(input = ., username = user, elevation_model = "astergdem") %>%
  as_tibble() %>%
  rename(elevation = elevation_geonames) %>%
  mutate(elevation = case_when (
    # There are about 10 grid cells that don't match to any
    # elevation; these are coded as -9999, but convert to NA
    elevation < 0 ~ NaN, 
    TRUE ~ elevation
  ))

write_csv(all_cells_el, "data/all_cells_el.csv")

# An alternate method using the raster package
jpn_el <- raster::getData('alt', country = "JPN")

el_values <-
raster::extract(
  jpn_el,
  all_cells %>% select(longitude, latitude) %>% as.matrix,
  method = "bilinear")

all_cells_with_el <-
  all_cells %>%
  mutate(el = el_values)

# There are over 1000 cells missing elevation data though, so
# don't use this.