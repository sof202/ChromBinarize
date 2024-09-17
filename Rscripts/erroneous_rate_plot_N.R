args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
bed_file_location <- args[2]
mark <- args[3]
max_read_depth <- as.numeric(args[4])
output_file <- args[5]

renv::load(renv_environment)
library(ggplot2)

## ============ ##
##   PLOTTING   ##
## ============ ##
concatenate_plots <- function(methylation_data,
                              max_read_depth,
                              percent_thresholds) {
  plot_list <- list()
  index <- 1
  for (percent in percent_thresholds) {
    p_plot <-
      chrombinarize::create_error_rate_plot(
        methylation_data,
        max_read_depth,
        percent
      )
    plot_list[[index]] <- p_plot
    index <- index + 1
  }
  p_plot_grid <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 3,
    labels = percent_thresholds,
    label_x = 0.65,
    label_y = 0.7
  )
  return(p_plot_grid)
}

## ======== ##
##   MAIN   ##
## ======== ##

methylation_data <- chrombinarize::read_bedmethyl(bed_file_location)

methylation_data <- dplyr::filter(methylation_data, mark_name == mark)

percent_thresholds <- c(0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)

p_plot_grid <-
  concatenate_plots(
    methylation_data,
    max_read_depth,
    percent_thresholds
  )

y_axis <- grid::textGrob("Probability of erroneous reads",
  gp = grid::gpar(fontface = "bold", fontsize = 15), rot = 90
)

x_axis <- grid::textGrob("Minimum read depth considered",
  gp = grid::gpar(fontface = "bold", fontsize = 15)
)

p_plot_grid <- gridExtra::grid.arrange(
  gridExtra::arrangeGrob(p_plot_grid, left = y_axis, bottom = x_axis)
)

# needed on some servers to actually create png files
options(bitmapType = "cairo")

ggsave(
  output_file,
  plot = p_plot_grid
)
