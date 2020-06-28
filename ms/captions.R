# Set up captions

# - Figures
fig_nums <- captioner::captioner(prefix = "Fig.")
fig_nums(name = "diversity", caption = "Taxonomic, phylogenetic, and morphological diversity of the ferns of Japan.")
fig_nums(name = "endemism", caption = "Phylogenetic endemism of the ferns of Japan.")
fig_nums(name = "bioregions", caption = "Bioregions of the ferns of Japan.")
fig_nums(name = "ecostructure", caption = "Ecological structure of the ferns of Japan.")

# - Tables
table_nums <- captioner::captioner(prefix = "Table")
table_nums(name = "traits-used", caption = "Fern traits used in this study")

# - SI figures
s_fig_nums <- captioner::captioner(prefix = "Fig. S", auto_space = FALSE, suffix = ": ")
s_fig_nums(name = "k-plot", "Selection of *k* for bioregions.")
s_fig_nums(name = "endemism-restricted", "Phylogenetic endemism of the ferns of Japan, restricted dataset including only taxa endemic to Japan.")

# - SI tables
# s_table_nums <- captioner::captioner(prefix = "Table S", auto_space = FALSE, suffix = ": ")
# s_table_nums(name = "si_table", "")

# Make short versions of citation functions

# - Just the number
figure <- pryr::partial(fig_nums, display = "cite")
table <- pryr::partial(table_nums, display = "cite")
s_figure <- pryr::partial(s_fig_nums, display = "cite")
# s_table <- pryr::partial(s_table_nums, display = "cite")

# - Just the caption
figure_cap <- function(x) {fig_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
table_cap <- function(x) {table_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
s_figure_cap <- function(x) {s_fig_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
# s_table_cap <- function(x) {s_table_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
