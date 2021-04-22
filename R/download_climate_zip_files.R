# Download climate data for Japan at the 3rd-degree mesh scale (1 km x 1 km).
#
# The 3rd-degree mesh scale data is provided in a series of zip files (151), 
# instead of all together for the whole country. These are a pain to download 
# manually, so grab URLs to the zip files by parsing Japan's National Land 
# Information Division website with the download links. 
#
# For more about Japan's mesh grid system for mapping, see
# http://www.stat.go.jp/english/data/mesh/05.html
#
# Also see for ref:
# https://suryu.me/post/mesh_code_mapping/

# Load packages
library(tidyverse)

# Read in the website with the download links for climate data
# (last run 2021-04-22)
site_text <- read_lines("https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-G02.html")

# Parse out URLs with the path to the zip files
urls <-
  site_text %>%
  magrittr::extract(str_detect(., "DownLd")) %>%
  str_match("(/ksj/gml/data/[^ ]+zip)") %>%
  magrittr::extract(,2) %>%
  paste0("https://nlftp.mlit.go.jp", .)

files <- str_split(urls, "/") %>% 
  map_chr(last) %>%
  paste0("data_raw/climate_zip_files/", .)

# Download the files
download.file(urls, files)
