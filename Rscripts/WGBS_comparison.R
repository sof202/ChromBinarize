args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
input_bed_file <- args[2] # Expects 7 columns
output_directory <- args[3]

renv::load(renv_environment)
library(ggplot2)

## ===================== ##
##   DATA MANIPULATION   ##
## ===================== ##
methylation_data <- chrombinarize::read_comparison_bedmethyl(input_bed_file)

# Here we combine 5mC and 5hmC as BS captures both signals without separating
# them
methylation_5mc <- dplyr::filter(
  methylation_data,
  mark_name == "m"
)
methylation_5hmc <- dplyr::filter(
  methylation_data,
  mark_name == "h"
) |>
  dplyr::mutate(
    mark_name = "m"
  )

combine_percentages <- function(reads1, reads2, percent1, percent2) {
  return((reads1 * percent1 + reads2 * percent2) / (reads1 + reads2))
}

methylation_data <- dplyr::full_join(
  methylation_5mc,
  methylation_5hmc,
  by = c(
    "chr",
    "start",
    "end",
    "mark_name",
    "BS_read_depth",
    "BS_percent_methylation"
  )
) |>
  dplyr::mutate(
    ONT_percent_methylation = combine_percentages(
      ONT_read_depth.x,
      ONT_read_depth.y,
      ONT_percent_methylation.x,
      ONT_percent_methylation.y
    ),
    ONT_read_depth = ONT_read_depth.x + ONT_read_depth.y
  ) |>
  dplyr::select(
    c(
      "chr",
      "start",
      "end",
      "mark_name",
      "ONT_read_depth",
      "ONT_percent_methylation",
      "BS_read_depth",
      "BS_percent_methylation"
    )
  )

methylation_data <- methylation_data |>
  dplyr::mutate(
    "absolute_change_percent_methylation" = abs(
      ONT_percent_methylation - BS_percent_methylation
    ),
    "absolute_change_read_depth" = abs(ONT_N - BS_N)
  ) |>
  dplyr::filter(
    ONT_N >= 30, # Filtering by read depth should give better results
    BS_N >= 30,
    mark == "m" # still have hydroxy rows at this point.
  )

# Purely for read depth histogram (makes plot more interpretable)
filtered_methylation_data <- methylation_data |>
  dplyr::filter(
    absolute_change_read_depth < 1000
  )

## ============ ##
##   PLOTTING   ##
## ============ ##

read_depth_plot <-
  ggplot(filtered_methylation_data, aes(x = absolute_change_read_depth)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(x = "absolute change in read depth")

percent_methylation_plot <-
  ggplot(methylation_data, aes(x = absolute_change_percent_methylation)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(x = "absolute change in percent methylation")

breaks <- c(1, 10, 100, 1000, 10000)

methylation_correlation_plot <-
  ggplot(
    methylation_data,
    aes(x = ONT_percent_methylation, y = BS_percent_methylation)
  ) +
  geom_bin2d(bins = 100) +
  scale_fill_gradient(
    name = "count", trans = "log",
    breaks = breaks, labels = breaks
  ) +
  theme_bw() +
  labs(x = "percent methylation in ONT", y = "percent methylation in WGBS")


options(bitmapType = "cairo")
ggsave(
  file.path(output_directory, "read_depth_change.png"),
  read_depth_plot
)
ggsave(
  file.path(output_directory, "percent_methylation_change.png"),
  percent_methylation_plot
)
ggsave(
  file.path(output_directory, "methylation_correlation.png"),
  methylation_correlation_plot
)
