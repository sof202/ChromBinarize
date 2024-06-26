# Required packages:
#  ggplot2
#  data.table
#  dplyr
#  cowplot
#  grid
#  gridExtra

args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
mark <- args[2] # m for 5mc and h for 5hmc
plot_type <- as.name(args[3])
output_file <- args[4]

library(ggplot2)
`%nin%` <- Negate(`%in%`)

## =============== ##
##   SUBSET DATA   ##
## =============== ##

subset_methylation_data <- function(methylation_data, read_threshold) {
  methylation_data <- methylation_data |>
    dplyr::filter(read_depth >= read_threshold)
  return(methylation_data)
}

## ============================ ##
##   PROBABILITY CALCULATIONS   ##
## ============================ ##

calculate_error_rate <- function(methylation_data, percent, plot_type) {
  if (plot_type == "p1") {
    methylation_data <- methylation_data |>
      dplyr::filter(percent_methylation >= percent) |>
      dplyr::mutate(
        "incorrectly_classified" = read_depth * (100 - percent_methylation) / 100
      )
  } else {
    methylation_data <- methylation_data |>
      dplyr::filter(percent_methylation <= percent) |>
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

create_p_plot_data <- function(methylation_data, plot_type) {
  if (plot_type == "p1") {
    percent_values <- 60:100
    p1 <- lapply(percent_values, function(percent) {
      calculate_error_rate(methylation_data, percent, plot_type)
    })
    stats_table <- data.table::data.table(
      "percent" = percent_values,
      "p1" = unlist(p1)
    )
    return(stats_table)
  }

  percent_values <- 40:0
  p2 <- lapply(percent_values, function(percent) {
    calculate_error_rate(methylation_data, percent, plot_type)
  })
  stats_table <- data.table::data.table(
    "percent" = percent_values,
    "p2" = unlist(p2)
  )
  return(stats_table)
}

create_p_plot <- function(methylation_data, plot_type) {
  stats_table <- create_p_plot_data(methylation_data, plot_type)
  p_plot <- ggplot(stats_table, aes(x = percent, y = !!plot_type)) +
    geom_point(color = "black") +
    labs(
      x = "",
      y = ""
    ) +
    theme_bw()
  return(p_plot)
}

concatenate_plots <- function(n_thresholds, methylation_data, plot_type) {
  plot_list <- list()
  index <- 1
  for (n_threshold in n_thresholds) {
    subsetted_methylation_data <-
      subset_methylation_data(methylation_data, n_threshold)
    p_plot <- create_p_plot(subsetted_methylation_data, plot_type)
    plot_list[[index]] <- p_plot
    index <- index + 1
  }
  p_plot_grid <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 3,
    labels = n_thresholds,
    label_x = 0.65,
    label_y = 0.9
  )
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

if (as.character(plot_type) %nin% c("p1", "p2")) {
  stop("Please specify a plot type from p1 or p2")
}

if (mark == "m") {
  n_thresholds <- c(100, 150, 200, 250, 300, 350, 400, 450, 500)
} else {
  n_thresholds <- c(10, 20, 30, 40, 50, 60) # hydroxy is usually with less reads
}

p_plot_grid <- concatenate_plots(n_thresholds, methylation_data, plot_type)

y_axis <- grid::textGrob("Probability of erroneous Read",
  gp = grid::gpar(fontface = "bold", fontsize = 15), rot = 90
)

x_axis <- grid::textGrob("Minimum percent methylated considered",
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
