# Download zipped data files of conservation areas ----
# 2021-08-26
# Website with protected areas is available here https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-A45.html
# It includes Forest ecosystem protection areas (森林生態系保護地域), which are not included in data
# downloaded from https://www.biodic.go.jp/biodiversity/activity/policy/map/map17/index.html
# The 3rd-degree mesh scale data is provided in a series of zip files, instead of all together
# for the whole country. These are a pain to download manually, so grab URLs to the zip files by
# parsing Japan's National Land Information Division website with the download links. 
# 
# Usage agreement: https://nlftp.mlit.go.jp/ksj/other/agreement.html

source("R/packages.R")
source("R/functions.R")

# Read in the website with the download links for natural protected areas data
site_text <- read_lines("https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-A45.html")

# Parse out URLs with the path to the zip files
# NOTE: this works as of 2021-08-26, but it may not if the website changes
urls <-
  site_text %>%
  magrittr::extract(str_detect(., "DownLd")) %>%
  str_match("(/data/A45[^ ]+zip)") %>%
  magrittr::extract(,2) %>%
  paste0("https://nlftp.mlit.go.jp/ksj/gml", .)

files <- str_split(urls, "/") %>% 
  map_chr(last) %>%
  paste0("data_raw/forest_area_zip_files/", .)

# Download the files
# Watch for warnings about partially downloaded files!
# These may need to be downloaded manually
download.file(urls, files)
