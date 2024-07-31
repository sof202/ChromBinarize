args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
bed_file_location <- args[2]
min_distance <- as.numeric(args[3])
max_distance <- as.numeric(args[4])
output_file <- args[5]

renv::load(renv_environment)
library(ggplot2)

## ===================== ##
##   DATA MANIPULATION   ##
## ===================== ##

read_methylation_data <- function(bed_file_location) {
  methylation_data <- data.table::fread(bed_file_location)
  names(methylation_data) <-
    c("Chr", "start", "end", "name", "size", "strand", "percent_methylation")
  return(methylation_data)
}

add_lead_and_lag <- function(methylation_data) {
  methylation_data <- dplyr::mutate(
    methylation_data,
    "leading_percent_methylation" = dplyr::lead(percent_methylation),
    "leading_start_position" = dplyr::lead(start),
    "lagging_percent_methylation" = dplyr::lag(percent_methylation),
    "lagging_start_position" = dplyr::lag(start)
  )
  return(methylation_data)
}

exclude_distant_neighbours <- function(start,
                                       neighbour_start,
                                       max_distance,
                                       percent_methylation) {
  distance <- abs(start - neighbour_start)
  return(ifelse(
    distance < max_distance && distance > min_distance,
    percent_methylation,
    NA
  ))
}

calculate_mean_methylation <- function(lead_methylation, lag_methylation) {
  return(mean(c(lead_methylation, lag_methylation), na.rm = TRUE))
}

## ============ ##
##   PLOTTING   ##
## ============ ##

create_cpg_robustness_plot <- function(methylation_data) {
  methylation_correlation <- cor(
    methylation_data[["percent_methylation"]],
    methylation_data[["average_surrounding_methylation"]]
  )
  cpg_robustness_plot <-
    ggplot(
      methylation_data,
      aes(x = percent_methylation, y = average_surrounding_methylation)
    ) +
    geom_bin2d(bins = 100) +
    scale_fill_continuous(type = "viridis") +
    geom_abline(intercept = 0, slope = 1, color = "black") +
    labs(
      x = "Percent methylation",
      y = "Average surrounding percent methylation"
    ) +
    scale_fill_gradientn(colors = c("grey", "red", "yellow")) +
    annotate(
      "text",
      x = 95,
      y = 5,
      label = paste0("R^2: ", methylation_correlation)
    ) +
    theme_bw()

  return(cpg_robustness_plot)
}

## ======== ##
##   MAIN   ##
## ======== ##

methylation_data <- read_methylation_data(bed_file_location)

methylation_data <-
  add_lead_and_lag(methylation_data) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    "true_lead_percent_methylation" = exclude_distant_neighbours(
      start,
      leading_start_position,
      max_distance,
      leading_percent_methylation
    )
  ) |>
  dplyr::mutate(
    "true_lag_percent_methylation" = exclude_distant_neighbours(
      start,
      lagging_start_position,
      max_distance,
      lagging_percent_methylation
    )
  ) |>
  dplyr::mutate(
    "average_surrounding_methylation" = calculate_mean_methylation(
      true_lead_percent_methylation,
      true_lag_percent_methylation
    )
  ) |>
  dplyr::filter(
    !is.na(average_surrounding_methylation)
  )

cpg_robustness_plot <- create_cpg_robustness_plot(methylation_data)
write.table(
  methylation_data,
  quote = FALSE,
  row.names = FALSE
)

# needed on some servers to actually create png files
options(bitmapType = "cairo")

ggsave(
  output_file,
  cpg_robustness_plot
)
