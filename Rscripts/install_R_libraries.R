required_packages <- c("dplyr", "data.table", "cowplot", "grid", "gridExtra")
installed_packages <- rownames(installed.packages())

install.packages(setdiff(required_packages, installed_packages))
