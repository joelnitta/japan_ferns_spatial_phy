# water_area_data
#
# Download area of water surfaces in Japan at 10 x 10 km scale

# Load packages
# kokudosuuchi is a package to interact with 
# Japan's National Land Information Division data API
library(kokudosuuchi) 
library(tidyverse)

# Set working directory
here::here()

# First get list of all available datasets on Kokudo suuchi 
# (National Land Information Division of Japan)
ksj_metadata <- getKSJSummary()

# Identify the dataset we're interested in, "Land usage third-degree mesh data"
dataset_id <- filter(ksj_metadata, stringr::str_detect(title, "土地利用3次メッシュ")) %>%
  pull(identifier)
# Get a tibble of URLs to zip files with this data
urls <- getKSJURL(dataset_id)

# Download all sets of data, read them in, 
# and extract the relevant data into a dataframe.
#
# This will download the data in chunks as zip files to the
# 'cached_zip' folder.
#
# These are 4th-degree mesh (1 km x 1 km), so there's a lot 
# (almost 2000 zip files, each with hundreds of mesh grids).
# We are interested in resolution at the 3rd degree mesh level (10 km x 10 km),
# but the data are only available at 4th degree. So we will sum these
# by 3rd degree mesh ID to convert them to the proper grain.
#
# Description of columns for this dataset:
# http://nlftp.mlit.go.jp/ksj/gml/codelist/LandUseCd-09.html
# Description of the dataset:
# http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-L03-a.html
# Also see for ref:
# https://suryu.me/post/mesh_code_mapping/
# 
# The raw data contain 16 columns. We are interested in L03a_014 and L03a_015,
# which are square meters of surface for coast and open water.

area_data <- map(urls$zipFileUrl, ~getKSJData(.,  cache_dir = "cached_zip")) %>%
  purrr::flatten(.) %>%
  map_df( 
    ~as_tibble(.) %>% select(code = 1, contains("L03a_014"),contains("L03a_015"))
  ) %>%
  rename(coast = L03a_014, water = L03a_015) %>%
  # Missing data were apparently coded as 9999999. Consider this 0 water surface area.
  mutate_at(vars(coast, water), ~str_replace_all(., "9999999", "0")) %>%
  mutate_at(vars(coast, water), as.numeric) %>%
  # Add secondary grid code: the first 6 digits of the primary grid code
  mutate(secondary_grid_code = substr(code, 1, 6)) %>%
  select(-code) %>%
  # Group by secondary ID code and take the sum
  group_by(secondary_grid_code) %>%
  summarize_at(vars(coast, water), ~sum(., na.rm = TRUE))

# Write out the data
write_csv(area_data, "data/2_grid_cells_water_area.csv")
