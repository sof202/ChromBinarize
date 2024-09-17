args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
input_bed_file <- args[2]
output_directory <- args[3]

renv::load(renv_environment)

## ======== ##
##   MAIN   ##
## ======== ##

methylation_data <- chrombinarize::read_comparison_bedmethyl(input_bed_file)


read_depth_plot <-
  chrombinarize::create_read_depth_plot(methylation_data)
percent_methylation_plot <-
  chrombinarize::create_read_depth_plot(methylation_data)
methylation_correlation_plot <-
  chrombinarize::create_correlation_plot(methylation_data)

options(bitmapType = "cairo")
ggplot2::ggsave(
  file.path(output_directory, "read_depth_change.png"),
  read_depth_plot
)
ggplot2::ggsave(
  file.path(output_directory, "percent_methylation_change.png"),
  percent_methylation_plot
)
ggplot2::ggsave(
  file.path(output_directory, "methylation_correlation.png"),
  methylation_correlation_plot
)
