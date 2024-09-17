args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
input_bed_file <- args[2] # Expects 7 columns
output_directory <- args[3]

renv::load(renv_environment)

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

## ============ ##
##   PLOTTING   ##
## ============ ##

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
