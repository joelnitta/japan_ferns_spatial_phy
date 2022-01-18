# Manuscript rendering functions and variables

# Themes ----

map_theme <- function() {
  theme(
    panel.grid.major = element_line(color = "grey60", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "grey60"),
    axis.title = element_blank(),
    legend.position = "bottom"
  )
}

# Color palettes ----

# - Coastline width across all plots
coast_lwd <- 0.05
# - Coast color across all plots
coast_col <- "grey50"

# - Randomization test significance (colorbrewer CVD safe)
colorb_paired <- brewer.pal(n = 8, name = "Paired")
colorb_ylgn <- brewer.pal(n = 5, name = "YlGn")
signif_cols <-
  c(
    "> 0.99" = colorb_paired[[2]], # Dark blue 
    "> 0.975" = colorb_paired[[1]], # Light blue 
    "< 0.025" = colorb_paired[[5]], # Light red
    "< 0.01" = colorb_paired[[6]], #  Dark red
    "not significant" = colorb_ylgn[[1]] # Beige
  )

# - CANAPE (Okabe-Ito CVD safe)
canape_cols <-
  c(
    "neo" = "#D55E00", # red
    "paleo" = "#0072B2", # dark blue
    "mixed" = "#009E73", # green
    "super" = "#F0E442", # yellow
    "not significant" = "grey97" # grey, almost white
  )

# - Bioregions (Okabe-Ito CVD safe)
bioregion_cols <- c( 
  "1" = "#009E73", # green
  "2" = "#F0E442", # yellow
  "3" = "#0072B2", # dark blue
  "4" = "#D55E00", # red
  "5" = "#56B4E9", # light blue
  "6" = "#E69F00", # goldenrod
  "7" = "#CC79A7", # magenta
  "8" = "#000000", # black
  "Other" = "grey80")

# - Protected areas (Okabe-Ito CVD safe)
protection_cols <- c(
  "Medium" = "#F0E442", # yellow
  "High" = "#009E73", # green
  "Total" = "#56B4E9" # light blue
)

protection_lines <- c(
  "Medium" = "dashed", 
  "High" = "dotted",
  "Total" = "solid"
)

# - Deer range (Okabe-Ito CVD safe)
deer_cols <- c( 
  "1978" = "#0072B2", # dark blue
  "2003" = "#F0E442", # yellow
  "Estimated" = "#D55E00" # red
) 

deer_lines <- c( 
  "1978" = "dotted", # dark blue
  "2003" = "dashed", # yellow
  "Estimated" = "solid" # red
) 

# Words ----
# Define formatting for some common custom words that may vary in style between journals
ie <- "*ie*"
eg <- "*eg*"

#' Format R packages names
#' 
#' J. Biogeography: Packages in R should in roman and quotations (e.g. 'vegan') and the relevant reference provided)
#'
#' @param x Name of R package
#'
#' @return Name of R package formatted according to journal requirements
#' 
pack <- function(x) {glue::glue("'{x}' v.{packageVersion(x)}")}

#' Format other software names
#' 
#' J. Biogeography: All software programs should be written in small caps, first written in roman (e.g. MrBayes or BEAST)
#'
#' @param x Name of software
#'
#' @return Name of software formatted according to journal requirements
#' 
# software <- function(x) {glue::glue('<span style="font-variant:small-caps;">{x}</span>')}
software <- function(x) {x}

#' Format R function names
#'
#' @param x Name of R function
#'
#' @return Name of R function formatted according to journal requirements
#' 
func <- function(x) {glue::glue("'{x}'")}

# Pagebreaks ----

#' Pagebreak for MS Word (doc) only
#' 
#' Meant for use within Rmd when rendering.
#' Breaks a page depending if the output format is MS Word
#' 
#' @param rmd_params Parameters set in YAML header or rmarkdown::render(). 
#' Must include `doc_type` (either 'doc' for MS Word or 'pdf' for PDF output)
#' 
pagebreak_doc <- function(rmd_params = params) {ifelse(rmd_params$doc_type == "doc", return("\\clearpage"), return(""))}

#' Pagebreak for PDF only
#' 
#' Meant for use within Rmd when rendering. Requires parameter variable
#' `doc_type` to be defined (either 'doc' for MS Word or 'pdf' for PDF output)
#'
#' Breaks a page depending if the output format is PDF
pagebreak_pdf <- function(rmd_params = params) {ifelse(rmd_params$doc_type == "pdf", return("\\clearpage"), return(""))}

#' Insert pagebreak
#' 
#' Meant for use within Rmd when rendering. Requires parameter variable
#' `doc_type` to be defined (either 'doc' for MS Word or 'pdf' for PDF output)
#'
#' Breaks a page depending if the output format is PDF or MS Word
pagebreak <- function(rmd_params = params) {
  if (rmd_params$doc_type == "doc") {
    run_pagebreak()
  } else if (rmd_params$doc_type == "pdf") {
    return("\\clearpage")
  } else {
    stop("doc_type parameter must be 'doc' or 'pdf'")
  }
}

# Figure output ----

#' Generate a path to save a results file
#' 
#' Only works for figures or tables cited in text, and outputs
#' files to the "results" folder
#'
#' @param result_num Number of result, e.g. "Fig. 2"
#' @param extension Extension to use for file
#'
#' @return String.
#' 
result_file <- function (result_num, extension) {
  
  fs::path(
    here::here("results"),
    result_num %>%
      str_remove_all("\\.") %>% 
      str_replace_all(" ", "_")
  ) %>%
    fs::path_ext_set(extension)
}

# Captions ----

# Follow Journal of Plant Research (JPR) style:
# - All figures are to be numbered using Arabic numerals.
# - Figures should always be cited in text in consecutive numerical order.
# - Figure parts should be denoted by lowercase letters (a, b, c, etc.).
# - If an appendix appears in your article and it contains one or more figures, 
#   continue the consecutive numbering of the main text. 
# - Do not number the appendix figures, "A1, A2, A3, etc." 
# - Figures in online appendices (Electronic Supplementary Material) should, 
#   however, be numbered separately

# - First define the "full" version, which would include a caption
# (except I never use the caption in the function, and instead replace with 'blank')
figure_full <- captioner::captioner(prefix = "Fig. ", auto_space = FALSE)
table_full <- captioner::captioner(prefix = "Table ", auto_space = FALSE)
s_figure_full <- captioner::captioner(prefix = "Fig. S1.", auto_space = FALSE)
s_table_full <- captioner::captioner(prefix = "Table S1.", auto_space = FALSE)

# - Make a short function that prints only the object type and number, e.g., "Fig. 1"
figure <- pryr::partial(figure_full, display = "cite", caption = "blank")
table <- pryr::partial(table_full, display = "cite", caption = "blank")
s_figure <- pryr::partial(s_figure_full, display = "cite", caption = "blank")
s_table <- pryr::partial(s_table_full, display = "cite", caption = "blank")

# - Make a short function that prints only the number (e.g., "1")
figure_num <- function (name) {figure(name) %>% str_remove("Fig. ")}
table_num <- function (name) {table(name) %>% str_remove("Table ")}
s_figure_num <- function (name) {s_figure(name) %>% str_remove("Fig. ")}
s_table_num <- function (name) {s_table(name) %>% str_remove("Table ")}

# Etc ----

#' Rename response variables in data frame
#' 
#' For formatting results tables
#'
#' @param df 
#'
#' @return df with renamed variables
#' 
rename_resp_vars <- function(df) {
  
  lookup <-
    tibble(
      resp_var = c("richness", "fd_obs_z", "pd_obs_z", "rfd_obs_z", "rpd_obs_z"),
      resp_var_print = c("Richness", "SES of FD", "SES of PD", "SES of RFD", "SES of RPD")
    )
  
  df %>%
    left_join(lookup, by = "resp_var") %>%
    select(-resp_var) %>%
    rename(resp_var = resp_var_print)
  
}

#' Rename model type in data frame
#' 
#' For formatting results tables
#'
#' @param df 
#'
#' @return df with renamed model types
#' 
rename_model_type <- function(df) {
  
  lookup <-
    tibble(
      model_type = c("percent_apo", "temp"),
      model_type_print = c("Reproductive", "Environmental")
    )
  
  df %>%
    left_join(lookup, by = "model_type") %>%
    select(-model_type) %>%
    rename(model_type = model_type_print)
  
}

#' Clean a vector of variables for printing
#'
#' @param x Character vector (variables)
#'
#' @return Character vector
#' 
clean_vars <- function(x) {
  x %>%
    gsub("^([^_]+)_obs_z", "SES of \\U\\1", ., perl=TRUE) %>%
    str_replace_all("^percent_apo$", "% apomictic taxa") %>%
    str_replace_all("^temp$", "Temperature") %>%
    str_replace_all("^temp_season$", "Temperature seasonality") %>%
    str_replace_all("^precip_season$", "Precipitation seasonality") %>%
    str_replace_all("^precip$", "Precipitation") %>%
    str_replace_all("^lat_area$", "Area") %>%
    str_replace_all("^richness$", "Richness")
}

#' Make a tibble matching the nodes of a dendrogram to the cluster they belong to
#'
#' @param dendro Dendrogram; list of class "phylo". Dendrogram clustering sites. Part of output of phyloregion::phyloregion()
#' @param tax_clusters Tibble with two columns, "grids" (site names) and "cluster" (ID of cluster for each site)
#' @param cluster_select Value indicating the cluster to select for making the tibble
#'
#' @return Tibble
#' 
map_cluster_to_nodes <- function(dendro, tax_clusters, cluster_select) {
  
  # Get vector of tips matching the selected cluster
  tips <- tax_clusters %>% filter(cluster == cluster_select) %>% pull(grids)
  
  if(length(tips) == 1) {
    # If only one tip, just return that node (tip) number
    nodes <- which(dendro$tip.label == tips)
  } else if (length(tips) > 1) {
    # If multiple tips, get MRCA, then all descendents
    node_mrca <- ape::getMRCA(dendro, tips)
    nodes <- phytools::getDescendants(dendro, node_mrca)
  } else {
    stop(glue::glue("Can't find cluster '{cluster_select}'"))
  }
  
  # Return a tibble mapping nodes to the cluster
  tibble(node = nodes, cluster = cluster_select)
}

# Define function for plotting model fit and actual data
plot_fits <- function(resp_var, resp_var_print, indep_var, signif_var, fit, model_data, signif_cols, ...) {
  
  # Basic plot setup
  p <- ggplot(model_data, aes(x = .data[[indep_var]], y = .data[[resp_var]])) +
    labs(y = resp_var_print) +
    theme_gray(base_size = 9) +
    theme(panel.grid.minor = element_blank())
  
  # Add points, color by significance if not plotting richness
  if (resp_var != "richness") {
    p <- p + 
      geom_point(aes(color = .data[[signif_var]]), size = 0.5, shape = 16) +
      scale_color_manual(values = signif_cols)
  } else {
    p <- p + 
      geom_point(size = 0.5, color = "grey20", alpha = 0.5, shape = 16)
  }
  
  # Add line and ribbon of model fit
  p <- p +
    geom_line(data = fit, alpha = 0.35, color = "blue") +
    geom_ribbon(
      data = fit,
      aes(ymin = low, ymax = up), 
      alpha = 0.15
    ) +
    theme(legend.position = "none")
  
  # Format x-axis ticks and labels
  if (indep_var == "temp") {
    p <- p +
      scale_x_continuous(labels = function(x) x*0.1) +
      labs(x = "Temperature (°C)")
  } else if (indep_var == "percent_apo") {
    p <- p +
      scale_x_continuous(labels = function(x) scales::percent(x, accuracy = 1)) +
      labs(x = "% apomictic taxa")
  }
  
  p
  
}

#' Format text for README description of rectangular data
#'
#' @param data Dataframe
#' @param metadata Tibble with two columns, "col_names" (column names in `data`)
#' and "desc" (description of each column). "col_names" must match column names
#' in `data` exactly (including order)
#'
#' @return Character string formatted for including in data README
#' 
make_data_desc <- function(data, metadata) {

  col_names <- colnames(data)
  
  assertthat::assert_that(
    isTRUE(all.equal(col_names, metadata$col))
  )
  
  metadata %>%
    mutate(text = glue::glue("- {col}: {desc}\n\n")) %>%
    pull(text)

}