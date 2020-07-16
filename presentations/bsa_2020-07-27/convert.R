# Convert html slides output by xaringan to pdf, then powerpoint

# requires decktape to be installed

# see https://github.com/gadenbuie/xaringan2powerpoint

setwd(here::here("presentations/bsa_2020-07-27/"))

# Specify slides html file
slides_html <- "index.html"

# "print" HTML to PDF
xaringan::decktape(slides_html, "slides.pdf", docker = FALSE)

# how many pages?
pages <- pdftools::pdf_info("slides.pdf")$pages

# set filenames
filenames <- sprintf("slides/slides_%02d.png", 1:pages)

# create slides/ and convert PDF to PNG files
# overwrite any existing folder: use carefully!
fs::dir_delete("slides")
fs::dir_create("slides")

pdftools::pdf_convert("slides.pdf", filenames = filenames, dpi = 150)

# Template for markdown containing slide images
slide_images <- glue::glue("
---

![]({filenames}){{width=100%, height=100%}}

")
slide_images <- paste(slide_images, collapse = "\n")

# R Markdown -> powerpoint presentation source
md <- glue::glue("
---
output: powerpoint_presentation
---

{slide_images}
")

cat(md, file = "slides_powerpoint.Rmd")

# Render Rmd to powerpoint
rmarkdown::render("slides_powerpoint.Rmd")  ## powerpoint!

# Cleanup
fs::file_delete("slides_powerpoint.Rmd")
fs::dir_delete("slides")

