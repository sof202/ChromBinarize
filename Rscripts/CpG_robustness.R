args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
bed_file_location <- args[2]
min_distance <- as.numeric(args[3])
max_distance <- as.numeric(args[4])
output_file <- args[5]

renv::load(renv_environment)
library(ggplot2)

## ======== ##
##   MAIN   ##
## ======== ##

methylation_data <- chrombinarize::read_bedmethyl(bed_file_location)

cpg_robustness_plot <- chrombinarize::create_cpg_robustness_plot(
  methylation_data,
  min_distance,
  max_distance
)

# needed on some servers to actually create png files
options(bitmapType = "cairo")

ggsave(
  output_file,
  cpg_robustness_plot
)
