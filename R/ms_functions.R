# Manuscript rendering functions and variables ----

# Define formatting for some common custom words that may vary in style between journals
ie <- "*ie*"
eg <- "*eg*"

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

#' Format R packages names
#'
#' @param x Name of R package
#'
#' @return Name of R package formatted according to journal requirements
#' 
pack <- function(x) {glue::glue("'{x}'")}

#' Format R function names
#'
#' @param x Name of R function
#'
#' @return Name of R function formatted according to journal requirements
#' 
func <- function(x) {glue::glue("'{x}'")}

#' Pagebreak for MS Word (doc) only
#' 
#' Meant for use within Rmd when rendering.
#' Breaks a page depending if the output format is MS Word
#' 
#' @param rmd_params Parameters set in YAML header or rmarkdown::render(). 
#' Must include `doc_type` (either 'doc' for MS Word or 'pdf' for PDF output)
#' 
pagebreak_doc <- function(rmd_params = params) {ifelse(rmd_params$doc_type == "doc", run_pagebreak(), return(""))}

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
figure_full <- captioner::captioner(prefix = "Fig.", suffix = "")
table_full <- captioner::captioner(prefix = "Table")
s_figure_full <- captioner::captioner(prefix = "Fig. S", auto_space = FALSE, suffix = "")
s_table_full <- captioner::captioner(prefix = "Table S", auto_space = FALSE, suffix = "")

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
