# Script to download and prepare environment data for Japan at the
# 2nd-degree mesh (ca. 10km x 10km) scale
#
# This will download the data in chunks as zip files to the
# 'cached_zip' folder temporarily. The final data file will be written
# to 'data/ja_env_data.csv'.
#
# All raw data are at the 3rd-degree mesh scale (1 km x 1 km).
# We are interested in resolution at the 2nd degree mesh level (10 km x 10 km),
# but the data are only available at 3rd degree. So we will sum these
# by 2nd degree mesh ID to convert them to the proper grain.
#
# For more about Japan's mesh grid system for mapping, see
# http://www.stat.go.jp/english/data/mesh/05.html
#
# Also see for ref:
# https://suryu.me/post/mesh_code_mapping/

# Load packages
# kokudosuuchi is a package to interact with 
# Japan's National Land Information Division data API
library(kokudosuuchi) 
library(assertr)
library(tidyverse)

#' Get a list of URLs to download data for a kokudosuuchi dataset
#'
#' @param dataset Name of the dataset
#'
#' @return Character vector; list of URLs for that dataset
#'
#' @examples
#' get_ksj_urls("標高・傾斜度3次メッシュ")
get_ksj_urls <- function (dataset) {
  # Get all KSJ metadata
  ksj_metadata <- kokudosuuchi::getKSJSummary()
  # Get identifier for the dataset we're interested in
  dataset_id <- dplyr::filter(ksj_metadata, stringr::str_detect(title, dataset)) %>%
    dplyr::pull(identifier)
  # Get a tibble of URLs to zip files with this data
  urls <- kokudosuuchi::getKSJURL(dataset_id)
  urls$zipFileUrl
}

# Area data ----

# The raw data contain 16 columns. We are interested in L03a_014 and L03a_015,
# which are square meters of surface for coast and open water.
#
# Description of columns for this dataset:
# http://nlftp.mlit.go.jp/ksj/gml/codelist/LandUseCd-09.html
# Description of the dataset:
# http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-L03-a.html

area_data_raw_all <- 
  get_ksj_urls("土地利用3次メッシュ") %>%
  map(~getKSJData(.,  cache_dir = "cached_zip")) %>%
  purrr::flatten(.)

area_data_raw <-
  area_data_raw_all %>%
  map_df( 
    ~as_tibble(.) %>% select(
      code = 1, 
      contains("L03a_014"), # Coast surface area (sq m)
      contains("L03a_015") # Open water surface area (sq m)
      )
  ) %>%
  rename(coast = L03a_014, water = L03a_015)

area_data <-
  area_data_raw %>%
  # Missing data were apparently coded as 9999999. Consider this 0 water surface area.
  mutate_at(vars(coast, water), ~str_replace_all(., "9999999", "0")) %>%
  mutate_at(vars(coast, water), as.numeric) %>%
  # Add secondary grid code: the first 6 digits of the primary grid code
  mutate(secondary_grid_code = substr(code, 1, 6)) %>%
  select(-code) %>%
  # Group by secondary ID code and take the sum
  group_by(secondary_grid_code) %>%
  summarize_at(vars(coast, water), ~sum(., na.rm = TRUE))

# Elevation data ----
# Includes elevation and slope

# Download all data, flatten into single list
elevation_data_raw_all <- 
  get_ksj_urls("標高・傾斜度3次メッシュ") %>%
  map(~getKSJData(.,  cache_dir = "cached_zip")) %>%
  purrr::flatten(.)

# Extract selected variables
elevation_data_raw <-
  elevation_data_raw_all %>%
  map_df( 
    ~as_tibble(.) %>% select(
      code = 1, 
      contains("G04a_002"), # Mean elevation
      contains("G04a_010") # Mean slope
    )
  ) %>%
  rename(elevation = G04a_002, slope = G04a_010)
  
# Convert to mean values by secondary grid
elevation_data <-
  elevation_data_raw %>%
  mutate_at(vars(elevation, slope), as.numeric) %>%
  # Add secondary grid code: the first 6 digits of the primary grid code
  mutate(secondary_grid_code = substr(code, 1, 6)) %>%
  select(-code) %>%
  # Group by secondary ID code and take the mean
  group_by(secondary_grid_code) %>%
  summarize_at(vars(elevation, slope), ~mean(., na.rm = TRUE))

# Climate data ----
# Monthly and yearly averages of temperature, rainfall, snowfall, sunlight
# For detailed description of columns:
# http://nlftp.mlit.go.jp/ksj/gml/codelist/ClimateCd.html

# Download all data, flatten into single list
climate_data_raw_all <-
  get_ksj_urls("平年値メッシュ") %>%
  map(~getKSJData(.,  cache_dir = "cached_zip")) %>%
  purrr::flatten(.)
  
# Extract selected variables. We will use the monthly means to calculate
# standard deviation in temp (sd in monthtly mean values)
climate_data_raw <-
  climate_data_raw_all %>%
  map_df( 
    ~as_tibble(.) %>% select(
      code = 1, 
      contains("G02_014"), # total yearly rainfall (0.1 mm)
      contains("G02_058"), # total yearly snowfall (cm)
      contains("G02_071"), # total yearly hours of sunlight (0.1 hr)
      contains("G02_084"), # total yearly radiation (0.1 MJm-2)
      contains("G02_017"), # mean temp Jan (0.1 degree centigrade)
      contains("G02_020"), # mean temp Feb
      contains("G02_023"), # mean temp Mar
      contains("G02_026"), # mean temp Apr
      contains("G02_029"), # mean temp May
      contains("G02_032"), # mean temp Jun
      contains("G02_035"), # mean temp Jul
      contains("G02_038"), # mean temp Aug
      contains("G02_041"), # mean temp Sep
      contains("G02_044"), # mean temp Oct
      contains("G02_047"), # mean temp Nov
      contains("G02_050"), # mean temp Dec
      contains("G02_053") # yearly mean temp
    )
  ) %>%
  rename(
    rain = G02_014,
    snow = G02_058,
    sun = G02_071,
    radiation = G02_084,
    temp_01 = G02_017,
    temp_02 = G02_020,
    temp_03 = G02_023,
    temp_04 = G02_026,
    temp_05 = G02_029,
    temp_06 = G02_032,
    temp_07 = G02_035,
    temp_08 = G02_038,
    temp_09 = G02_041,
    temp_10 = G02_044,
    temp_11 = G02_047,
    temp_12 = G02_050,
    temp_year = G02_053
  )

climate_data_select <-
  climate_data_raw %>%
  mutate_at(vars(-code), as.numeric) %>%
  # Check for missing data
  assert(not_na, everything()) %>%
  # Sanity check: rainfall between 0 and 9,000 mm (max record rainfall for Japan 8670 mm)
  assert(function (x) x < 90000 & x > 0, rain) %>% 
  # Sanity check: mean temp between -60 and 30 degrees (min record temp Japan -50)
  assert(function (x) x < 300 & x > -600, contains("temp")) %>%
  # Add secondary grid code: the first 6 digits of the primary grid code
  mutate(secondary_grid_code = substr(code, 1, 6))

# Caculate SD from monthly mean temperature at 3rd degree mesh scale
climate_data_temp_sd <-
  climate_data_select %>%
  select(code, contains("temp")) %>%
  select(-temp_year) %>%
  gather("month", "temp", -code) %>%
  # Group by tertiary ID code and take the SD of all months
  group_by(code) %>%
  summarize(
    temp_sd = sd(temp, na.rm = TRUE)
  ) %>%
  mutate(secondary_grid_code = substr(code, 1, 6)) %>%
  group_by(secondary_grid_code) %>%
  # Group by secondary ID code and take the mean
  summarize(
    temp_sd = mean(temp_sd, na.rm = TRUE)
  )
  
# Join temperature SD with other 3rd degree mesh scale variables
# and summarize as mean values at 2nd degree mesh scale.
climate_data <-
  climate_data_select %>%
  select(secondary_grid_code, rain, snow, sun, radiation, temp_year, -code) %>%
  # Group by secondary ID code and take the mean
  group_by(secondary_grid_code) %>%
  summarize_all(~mean(., na.rm = TRUE)) %>%
  ungroup %>%
  left_join(climate_data_temp_sd) %>%
  # Convert rain to mm and temp to degrees celsius
  # Other units as-is:
  # snow (cm)
  # radiation (0.1 MJm-2)
  # sunlight time (0.1 hr)
  mutate_at(vars(rain, temp_year, temp_sd), ~magrittr::multiply_by(., 0.1))

# Combine datasets ----
ja_env_data <- purrr::reduce(
  list(area_data, 
  elevation_data, 
  climate_data), 
  full_join, by = "secondary_grid_code")

write_csv(ja_env_data, "data/ja_env_data.csv")
