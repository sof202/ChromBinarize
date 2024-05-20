args <- commandArgs(trailingOnly = TRUE)
bin_count_data_location <- args[1]
output_folder <- args[2]

library(ggplot2)
library(magrittr)

smooth_signal <- function(bin_counts, iteration) {
  reference_column <- paste0("iteration_", iteration - 1)
  new_column_name <- paste0("iteration_", iteration)
  bin_counts <- bin_counts %>%
    dplyr::mutate(
      leading = dplyr::lead(!!as.symbol(reference_column)),
      lagging = dplyr::lag(!!as.symbol(reference_column))
    ) %>%
    dplyr::mutate(
      new = rowMeans(
        .[, c(!!reference_column, "leading", "lagging")],
        na.rm = TRUE
      )
    ) %>%
    dplyr::select(-c("leading", "lagging")) %>%
    dplyr::rename(!!new_column_name := new)
  return(bin_counts)
}

rle_post_processing <- function(row_counts) {
  names(row_counts) <- c("type", "count")

  # The following makes histrogram plots easier to understand
  row_counts$type <- as.factor(row_counts$type)
  row_counts <- dplyr::filter(row_counts, count < 50)
  return(row_counts)
}


run_length_encoding <- function(data) {
  run_counts <- data.table::data.table()
  row <- 1
  while (row < length(data)) {
    count <- 0
    if (data[[row]] > 0) {
      while (row < length(data) && data[[row]] > 0) {
        row <- row + 1
        count <- count + 1
      }
      run_counts <- rbind(run_counts, list("present", count), use.names = FALSE)
    } else {
      while (row < length(data) && data[[row]] == 0) {
        row <- row + 1
        count <- count + 1
      }
      run_counts <- rbind(run_counts, list("absent", count), use.names = FALSE)
    }
  }
  run_counts <- rle_post_processing(run_counts)
  return(run_counts)
}

## ======== ##
##   MAIN   ##
## ======== ##

bin_counts <- data.table::fread(bin_count_data_location)
bin_counts <- bin_counts |>
  dplyr::rename(
    "chr" = "V1",
    "start" = "V2",
    "end" = "V3",
    "iteration_0" = "V4"
  )

# RLE requires a numerical vector as input
# We smooth the signal to account for *small* gaps between methylated CpGs
for (i in 1:3) {
  bin_counts <- smooth_signal(bin_counts, i)
}

# TODO: Make this look nicer
unsmoothed_signal <- unlist(bin_counts[, "iteration_0"])
once_smoothed_signal <- unlist(bin_counts[, "iteration_1"])
twice_smoothed_signal <- unlist(bin_counts[, "iteration_2"])
thrice_smoothed_signal <- unlist(bin_counts[, "iteration_3"])

unsmoothed_runs <- run_length_encoding(unsmoothed_signal)
once_smoothed_runs <- run_length_encoding(once_smoothed_signal)
twice_smoothed_runs <- run_length_encoding(twice_smoothed_signal)
thrice_smoothed_runs <- run_length_encoding(thrice_smoothed_signal)


## ============ ##
##   PLOTTING   ##
## ============ ##

create_stacked_histogram <- function(signal_data, plot_title) {
  histogram_plot <-
    ggplot(signal_data, aes(x = count, fill = type)) +
    geom_histogram(binwidth = 1) +
    scale_fill_manual(values = c(
      "absent" = "lightgreen",
      "present" = "darkgreen"
    )) +
    labs(
      title = plot_title,
      x = "Run length",
      y = "Frequency"
    ) +
    scale_y_continuous(trans = "log10") +
    theme_light()

  return(histogram_plot)
}

ggsave(
  paste0(output_folder, "/unsmoothed_run_counts.png"),
  create_stacked_histogram(unsmoothed_runs, "Unsmoothed Runs")
)
ggsave(
  paste0(output_folder, "/once_smoothed_run_counts.png"),
  create_stacked_histogram(once_smoothed_runs, "Once Smoothed Runs")
)
ggsave(
  paste0(output_folder, "/twice_smoothed_run_counts.png"),
  create_stacked_histogram(twice_smoothed_runs, "Twice Smoothed Runs")
)
ggsave(
  paste0(output_folder, "/thrice_smoothed_run_counts.png"),
  create_stacked_histogram(thrice_smoothed_runs, "Thrice Smoothed Runs")
)
