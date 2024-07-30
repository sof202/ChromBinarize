required_packages <- c(
  "dplyr",
  "data.table",
  "fitdistrplus",
  "cowplot",
  "grid",
  "gridExtra"
)
installed_packages <- rownames(installed.packages())

install.packages(setdiff(required_packages, installed_packages))
