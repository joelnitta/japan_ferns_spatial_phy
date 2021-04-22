# Install latex packages using tinytex

# This would happen automatically anyways when rendering the Rmd to pdf, 
# but requires downloading packages and may not work during updates of Tex Live. 
# Better to install to the docker image once and keep them there.

tinytex::tlmgr_update()

latex_packages <- c(
  "amsmath",
  "atbegshi",
  "atveryend",
  "auxhook",
  "bigintcalc",
  "bitset",
  "booktabs",
  "etexcmds",
  "etoolbox",
  "euenc",
  "float",
  "fontspec",
  "geometry",
  "gettitlestring",
  "grffile",
  "hycolor",
  "hyperref",
  "iftex",
  "infwarerr",
  "intcalc",
  "kvdefinekeys",
  "kvoptions",
  "kvsetkeys",
  "latex-amsmath-dev",
  "letltxmacro",
  "ltxcmds",
  "mdwtools",
  "pdfescape",
  "pdftexcmds",
  "refcount",
  "rerunfilecheck",
  "setspace",
  "siunitx",
  "stringenc",
  "tipa",
  "unicode-math",
  "uniquecounter",
  "xcolor",
  "xunicode",
  "zapfding"
)

tinytex::tlmgr_install(latex_packages)