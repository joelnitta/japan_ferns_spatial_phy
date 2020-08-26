# Set up captions

# - Figures
fig_nums <- captioner::captioner(prefix = "Fig.")
fig_nums(name = "japan-map", caption = "Map of Japan showing names of places mentioned in the text. Names of the four main islands in bold.")
fig_nums(name = "raw-diversity", caption = "Raw taxonomic (**A**), phylogenetic (**B**), and functional (**C**) diversity of the ferns of Japan.")
fig_nums(name = "rand-diversity", caption = "Results of randomization test for phylogenetic and functional diversity of the ferns of Japan.")
fig_nums(name = "endemism", caption = "Phylogenetic endemism of the ferns of Japan measured using CANAPE (categorical analysis of neo- and paleo-endemism).")
fig_nums(name = "bioregions", caption = "Bioregions of the ferns of Japan.")
fig_nums(name = "div-by-tax-region", caption = "Phylogenetic and morphological diversity of the ferns of Japan by bioregion.")
fig_nums(name = "conserv-status", caption = "Percent of land area with protected status for grid cells with significantly high biodiversity.")

# - Tables
table_nums <- captioner::captioner(prefix = "Table")
table_nums(name = "traits-used", caption = "Fern traits used in this study.")

# - SI figures
s_fig_nums <- captioner::captioner(prefix = "Fig. S", auto_space = FALSE, suffix = ": ")
s_fig_nums(name = "sampling-curve", "Species collection curve for the ferns of Japan fit with the iNEXT package.")
s_fig_nums(name = "grain-size", "Effect of grain size on sampling redundancy.")
s_fig_nums(name = "abundance", caption = "Observed number of specimens (**A**) and sampling redundancy (**B**) per 0.2° grid cell in the ferns of Japan.")
s_fig_nums(name = "k-plot", "Selection of *k* for bioregions.")
s_fig_nums(name = "raw-biplot", "Relationships between observed functional and phylogenetic diversity and taxonomic richness in the ferns of Japan.")
s_fig_nums(name = "endemism-restricted", "Phylogenetic endemism of the ferns of Japan measured using CANAPE (categorical analysis of neo- and paleo-endemism), restricted dataset including only taxa endemic to Japan.")

# - SI tables
s_table_nums <- captioner::captioner(prefix = "Table S", auto_space = FALSE, suffix = ": ")
s_table_nums(name = "phy-sig", caption = "Phylogenetic signal in continuous functional traits of the ferns of Japan.")
s_table_nums(name = "phy-sig-binary", caption = "Phylogenetic signal in quantitative (binary) functional traits of the ferns of Japan.")

# Make short versions of citation functions

# - Just the number
figure <- pryr::partial(fig_nums, display = "cite")
table <- pryr::partial(table_nums, display = "cite")
s_figure <- pryr::partial(s_fig_nums, display = "cite")
s_table <- pryr::partial(s_table_nums, display = "cite")

# - Just the caption
figure_cap <- function(x) {fig_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
table_cap <- function(x) {table_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
s_figure_cap <- function(x) {s_fig_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
s_table_cap <- function(x) {s_table_nums(x) %>% stringr::str_match(": (.*)$") %>% magrittr::extract(,2)}
