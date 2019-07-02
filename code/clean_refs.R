# Clean references exported from Mendeley

library(magrittr)

jntools::clean_bib(here::here("manuscript/references_raw.bib")) %>% 
  write_lines(here::here("manuscript/references.bib"))