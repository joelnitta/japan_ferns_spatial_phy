# Filter a list of references in YAML format to those cited in an RMD file
#
# The YAML list should be exported from Zotero
# like this: file -> "export library" -> "Better CSL YAML"

library(tidyverse)
library(yaml)
library(assertr)

#' Extract citations (formatted like `[@key]`) from an Rmd file
#'
#' @param rmd_file Character vector; path to Rmd file
#'
#' @return Data frame with one column 'key'
#' 
extract_citations <- function(rmd_file) {
  read_lines(rmd_file) %>%
    stringr::str_split(" |;") %>% 
    unlist %>% 
    magrittr::extract(., stringr::str_detect(., "@")) %>% 
    stringr::str_remove_all("\\[|\\]|\\)|\\(|\\.$|,|\\{|\\}") %>% 
    magrittr::extract(., stringr::str_detect(., "^@|^-@")) %>% 
    stringr::str_remove_all("^@|^-@") %>% 
    unique %>% 
    sort %>%
    tibble(key = .)
}

#' Filter a list of references in YAML format to those occurring in an Rmd file
#'
#' @param rmd_file Character vector; Path to Rmd file(s)
#' @param yaml_in String or list; if string, the path to the YAML file to filter.
#' If list, should be result of reading in a YAML file with yaml::read_yaml()
#' @param yaml_out Path to write filtered YAML reference file
#' @param strict Logical; should input YAML file be required to include
#' all references in rmd file?
#'
#' @return NULL; externally, the filtered YAML will be written to `yaml_out`
#' 
filter_refs_yaml <- function(
  rmd_file, yaml_in = "ms/main_library.yaml", yaml_out = "ms/references.yaml",
  strict = FALSE) {
  
  # Parse RMD file and extract citation keys
  citations <- map_df(rmd_file, extract_citations)
    
  # Read in YAML including all references exported from Zotero
  if(inherits(yaml_in, "character")) {
    ref_yaml <- yaml::read_yaml(yaml_in)
  } else if(inherits(yaml_in, "list")) {
    ref_yaml <- yaml_in
  } else {
    stop("`yaml_in` must be a path to a YAML file (string) or a list read in with yaml::read_yaml()")
  }
  
  # Extract all citation keys from full YAML
  cite_keys_all <- map_chr(ref_yaml$references, "id") %>%
    tibble(
      key = .,
      order = 1:length(.)
    )
  
  # Check that all keys in the yaml are present in the input YAML
  missing <- citations %>% anti_join(cite_keys_all, by = "key")
  
  if (isTRUE(strict)) {
  assertthat::assert_that(
    nrow(missing) == 0,
    msg = glue::glue(
      "The following ref keys are present in the Rmd but missing from the input YAML: {paste(missing$key, collapse = ', ')}"
      ))
  }
  
  cite_keys_filtered <- citations %>% inner_join(cite_keys_all, by = "key")
  
  # Filter YAML to only those citation keys in the RMD
  ref_yaml_filtered <- list(references = ref_yaml$references[cite_keys_filtered$order])
  
  # Write out the YAML file
  yaml::write_yaml(ref_yaml_filtered, file = yaml_out)
}

# Prepare YAML references for filtering ----
# Combine references_other.yaml (manually entered yaml) with 
# main_library.yaml (yaml automatically exported from zotero)
# Prefer manually entered values if duplicated.
refs_other <- read_yaml("ms/references_other.yaml")

cite_keys_other <- map_chr(refs_other$references, "id") %>%
  tibble(key = ., order = 1:length(.))

refs_all <- read_yaml("ms/main_library.yaml")

# remove URLs from main library references
for(i in 1:length(refs_all$references)) {
  refs_all$references[[i]]$URL <- NULL
}

cite_keys_all <- map_chr(refs_all$references, "id") %>%
  tibble(key = ., order = 1:length(.))

cite_keys_keep <- cite_keys_all %>%
  anti_join(cite_keys_other, by = "key")

refs_combined <- list(
  references = c(
    refs_other$references,
    refs_all$references[cite_keys_keep$order])
)

# Fix "Jr." in author names: should be in 'suffix' field, not part of
# first or last name

# to make code easier to read,
# extract first nested item of refs_combined
refs <- refs_combined[[1]]

for(i in seq_along(refs)) {
  for(j in seq_along(refs[[i]]$author)) {
    family_has_jr <- FALSE
    given_has_jr <- FALSE
    if(
      !is.null(refs[[i]]$author[[j]]$family) &
      !is.null(refs[[i]]$author[[j]]$given)) {
      family_has_jr <- str_detect(refs[[i]]$author[[j]]$family, "Jr\\.")
      given_has_jr <- str_detect(refs[[i]]$author[[j]]$given, "Jr\\.")
      refs[[i]]$author[[j]]$family <- str_remove_all(
        refs[[i]]$author[[j]]$family, "Jr\\.") %>%
        str_squish()
      refs[[i]]$author[[j]]$given <- str_remove_all(
        refs[[i]]$author[[j]]$given, "Jr\\.") %>%
        str_squish()
      if(family_has_jr | given_has_jr) {
        refs[[i]]$author[[j]]$suffix <- "Jr."
        # uncomment to show changed name
        # message(glue::glue("Fixed {refs[[i]]$author}"))
      }
    }
  }
}

refs_combined[[1]] <- refs

# Filter YAML references ----
# - main MS
filter_refs_yaml(c("ms/manuscript.Rmd", "ms/SI.Rmd"), refs_combined, "ms/references.yaml")
# - data README (remember to remove `<span class="nocase">` manually)
# also for some reason Kusumoto2017 isn't included, so copy over manually
filter_refs_yaml("ms/data_readme.Rmd", refs_combined, "ms/data_readme_refs.yaml")

# Write out 'aux' file with citation keys to make collection in Zotero
# so reference data are easier to work with.
# The collection can be made with the betterbibtex for Zotero plugin:
# https://retorque.re/zotero-better-bibtex/citing/aux-scanner/
# For our purposes, the aux file just consists of lines formatted like
# \citation{KEY}
# where KEY is the citation key
extract_citations("ms/manuscript.Rmd") %>%
  mutate(latex = glue::glue("\\citation{[key]}", .open = "[", .close = "]")) %>%
  pull(latex) %>%
  write_lines("working/ms_keys.aux")
