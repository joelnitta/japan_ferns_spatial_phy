# Make clean bib files for each Rmd file that only includes cited references

library(jntools)
library(tidyverse)

make_ref_list(
  rmd_file = "ms/manuscript.Rmd", 
  raw_bib = "ms/references_raw.bib",
  final_bib = "ms/references.bib")

# Make some manual fixes to authors in SI bibliography
# (these are institutions, so need double brackets to
# avoid latex thinking they have first and last names)
read_lines("ms/references.bib") %>%
  str_replace(
    "R Core Team",
    "\\{R Core Team\\}") %>%
  write_lines("ms/references.bib")
