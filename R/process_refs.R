### Set up bibliography and citation style files ###

library(tidyverse)
library(assertr)
library(RefManageR)

# Filter bibtex bibliography ----

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
  unique %>% sort

# Read in entire reference library exported from Zotero with BetterBibTex
# (note: this can be a soft-link to ~/Dropbox/bibliography/main_library.bib)
bib <- RefManageR::ReadBib("ms/main_library.bib", check = FALSE)

# Convert reference library to dataframe
bib_df <- as.data.frame(bib) %>% rownames_to_column("cite_key")

# Make sure no citations are missing
tibble(cite_key = citations) %>%
  anti_join(bib_df, by = "cite_key") %>%
  verify(nrow(.) == 0, success_fun = success_logical, error_fun = error_stop)

# Set list of citations to exclude (websites etc since RefManageR doesn't like that format)
citations_exclude <- c("RCoreTeam2020", "Orme2018")

# Filter bibliography to only cited references, write it out
bib_df %>%
  filter(cite_key %in% citations) %>%
  filter(!cite_key %in% citations_exclude) %>% 
  # Filter out website references 
  # (treated separately in references_other.yaml)
  filter(bibtype != "Online") %>%
  # Set URL to NA for articles
  mutate(url = case_when(
    bibtype == "Article" ~ NA_character_,
    TRUE ~ url
  )) %>%
  # Select only needed columns
  select(
    cite_key, bibtype, type, doi,
    title, author, date, 
    journaltitle, volume, issue, number, pages, # article info
    booktitle, editor, publisher, edition, location, address, # book info
    pagetotal, institution, # thesis info
    url # website info
    ) %>%
  # RefManageR::as.BibEntry() needs rownames to be citation keys
  column_to_rownames("cite_key") %>%
  RefManageR::as.BibEntry() %>%
  RefManageR::WriteBib(file = "ms/references.bib")

read_lines("ms/references.bib") %>%
  # Fix "Jr" in Eddie's name (eg, "Stoner, Jr, Winifred Sackville" is correct format in bibtex)
  str_replace_all("James E \\{Watkins Jr\\.\\}", "Watkins, Jr, James E") %>%
  # Fix accents over vowels, which somehow got mangled by RefManageR
  str_replace_all("\\{\\\\a'a\\}", "\\{\\\\'{a}}") %>%
  str_replace_all("\\{\\\\a'e\\}", "\\{\\\\'{e}}") %>%
  str_replace_all("\\{\\\\a'i\\}", "\\{\\\\'{i}}") %>%
  str_replace_all("\\{\\\\a'i\\}", "\\{\\\\'{i}}") %>%
  str_replace_all("\\{\\\\a'\\\\i\\}", "\\{\\\\'{i}}") %>%
  str_replace_all("\\{\\\\a'o\\}", "\\{\\\\'{o}}") %>%
  str_replace_all("\\{\\\\a'u\\}", "\\{\\\\'{u}}") %>%
  write_lines("ms/references.bib")
