# Filter a list of references in YAML format to those cited in an RMD file
#
# The YAML list should be exported from Zotero
# like this: file -> "export library" -> "Better CSL YAML"

library(tidyverse)
library(yaml)
library(assertr)

# Set path to RMD file
rmd_file <- "ms/manuscript.Rmd"

# Parse RMD file and extract citation keys
citations <- 
  readr::read_lines(rmd_file) %>% 
  stringr::str_split(" |;") %>% 
  unlist %>% 
  magrittr::extract(., stringr::str_detect(., "@")) %>% 
  stringr::str_remove_all("\\[|\\]|\\)|\\(|\\.$|,|\\{|\\}") %>% 
  magrittr::extract(., stringr::str_detect(., "^@|^-@")) %>% 
  stringr::str_remove_all("^@|^-@") %>% 
  unique %>% 
  sort %>%
  tibble(key = .)

# Read in YAML including all references exported from Zotero
ref_yaml <- yaml::read_yaml("data_raw/main_library.yaml")

# Extract citation keys from YAML, filter to only those in the RMD
cite_keys <- map_chr(ref_yaml$references, "id") %>%
  tibble(
    key = .,
    order = 1:length(.)
  ) %>%
  inner_join(citations, by = "key") %>%
  # Stop if any keys are missing
  assert(not_na, everything())

# Filter YAML to only those citation keys in the RMD
ref_yaml_filtered <- list(references = ref_yaml$references[cite_keys$order])

# Write out the YAML file
yaml::write_yaml(ref_yaml_filtered, file = "ms/references.yaml")
