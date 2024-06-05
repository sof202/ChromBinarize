args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
mark <- args[2] # m for 5mc and h for 5hmc
max_read_depth <- as.numeric(args[3]) # you don't need to consider N = 2000
plot_type <- as.name(args[4])
output_file <- args[5]

library(ggplot2)

## =============== ##
##   SUBSET DATA   ##
## =============== ##

subet_methylation_data <- function(methylation_data, methylation_threshold, plot_type) {
  if (plot_type == "p1") {
    methylation_data <- methylation_data |>
      dplyr::filter(percent_methylation >= methylation_threshold)
    return(methylation_data)
  }
  methylation_data <- methylation_data |>
    dplyr::filter(percent_methylation <= methylation_threshold)
  return(methylation_data)
}

## ============================ ##
##   PROBABILITY CALCULATIONS   ##
## ============================ ##

calculate_error_rate <- function(methylation_data, n, plot_type) {
  if (plot_type == "p1") {
    methylation_data <- methylation_data |>
      dplyr::filter(read_depth >= n) |>
      dplyr::mutate(
      "incorrectly_classified" = read_depth * (100 - percent_methylation) / 100
      )
  } else {
    methylation_data <- methylation_data |>
      dplyr::filter(read_depth >= n) |>
      dplyr::mutate(
      "incorrectly_classified" = read_depth * percent_methylation / 100
      )
  }
  incorrect_total <- sum(methylation_data$incorrectly_classified)
  total_reads <- sum(methylation_data$read_depth)
  error_rate <- incorrect_total / total_reads

  return(error_rate)
}

## ====================== ##
##   PLOTTING FUNCTIONS   ##
## ====================== ##

create_p_plot_data <- function(max_read_depth, methylation_data, plot_type) {
  n_values <- 1:max_read_depth

  if (plot_type == "p1") {
    p1 <- lapply(n_values, function(n) {
      calculate_error_rate(methylation_data, n)
    })
    stats_table <- data.table::data.table(
      "n" = n_values,
      "p1" = unlist(p1)
    )
    return(stats_table)
  }

  p2 <- lapply(n_values, function(n) {
    calculate_error_rate(methylation_data, n)
  })
  stats_table <- data.table::data.table(
    "n" = n_values,
    "p2" = unlist(p2)
  )
  return(stats_table)
}

create_p_plot <- function(max_read_depth, methylation_data, plot_type) {
  stats_table <- create_plot_data(max_read_depth, methylation_data, plot_type)
  p_plot <- ggplot(stats_table, aes(x = n, y = !!plot_type)) +
    geom_point(color = "black") +
    labs(
      x = "min read depth considered",
      y = "probability of erroneous read",
      title = "probability of erroneous read for binomial distribution"
    ) +
    theme_bw()
  return(p_plot)
}

concatenate_plots <- function(max_read_depth, methylation_data, plot_type, percent_thresholds) {
  plot_list <- NULL
  for (percent in percent_thresholds) {
    subsetted_methylation_data <- subset_methylation_data(methylation_data, percent, plot_type)
    plot_list <- c(plot_list, create_p_plot(max_read_depth, subsetted_methylation_data, plot_type))
  }
  p_plot_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = 3, labels = percent_thresholds)
  return(p_plot_grid)
}

## ======== ##
##   MAIN   ##
## ======== ##

methylation_data <- data.table::fread(bed_file)

names(methylation_data) <- c(
  "chr",
  "start",
  "end",
  "mark_name",
  "read_depth",
  "strand",
  "percent_methylation"
)

methylation_data <- dplyr::filter(methylation_data, mark_name == mark)

if (plot_type == "p1") {
  percent_thresholds <- c(80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5, 100)
} else if (plot_type == "p1") {
  percent_thresholds <- c(0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)
} else {
  stop("please provide p1 or p2 as the plot type")
}

p_plot_grid <- concatenate_plots(max_read_depth, methylation_data, plot_type, percent_thresholds)

# needed on some servers to actually create png files
options(bitmapType = "cairo")

ggsave(
  output_file,
  plot = p_plot_grid
)
