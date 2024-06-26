args <- commandArgs(trailingOnly = TRUE)
input_bed_file <- args[1] # Expects 7 columns
output_directory <- args[2]

library(ggplot2)

## ===================== ##
##   DATA MANIPULATION   ##
## ===================== ##
methylation_data <- data.table::fread(input_bed_file)
colnames(methylation_data) <- c(
  "chr",
  "start",
  "end",
  "mark",
  "ONT_N",
  "ONT_percent_methylation",
  "oxBS_N",
  "oxBS_percent_methylation"
)

# Here we extract the hydroxymethylation
methylation_data <- methylation_data |>
  dplyr::filter(
    mark == "h"
  )

methylation_data <- methylation_data |>
  dplyr::mutate(
    "absolute_change_percent_methylation" = abs(
      ONT_percent_methylation - oxBS_percent_methylation
    ),
    "absolute_change_read_depth" = abs(ONT_N - oxBS_N)
  ) |>
  dplyr::filter(
    ONT_N >= 30, # Filtering by read depth should give better results
    oxBS_N >= 30,
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
    aes(x = ONT_percent_methylation, y = oxBS_percent_methylation)
  ) +
  geom_bin2d(bins = 100) +
  scale_fill_gradient(name = "count", trans = "log",
                        breaks = breaks, labels = breaks) +
  theme_bw() +
  labs(x = "percent methylation in ONT", y = "percent methylation in oxBS")


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
