#' @title Generate Data for `create_error_rate_plot()`
#'
#' @description Create a data table with the data required for the
#'  `create_error_rate_plot()` function.
#'
#' @inheritParams create_error_rate_plot
#'
#' @return A data table with columns
#'  - n: The minimum read depth considered (integer)
#'  - error_rate: The error rate estimated (numeric)
#'
#' @examples
#' create_error_rate_data(methylation_data, 1000, 5)
#'         n error_rate
#'     <int>      <num>
#'  1:     1          0.05
#'  2:     2          0.049
#'  3:     3          0.049
#'  4:     4          0.048
#'  5:     5          0.046
#'  6:     6          0.046
#'  7:     7          0.04
#'  8:     8          0.04
#'  9:     9          0.039
#' 10:    10          0.039
create_error_rate_data <- function(methylation_data,
                                   max_read_depth = 1000,
                                   percent_threshold = 5) {
  n_values <- 1:max_read_depth
  error_rates <- lapply(n_values, function(n) {
    methylation_data <- dplyr::filter(methylation_data, read_depth >= n)
    estimate_error_rate(methylation_data, n, percent_threshold)
  })
  error_rate_data <- data.table::data.table(
    "n" = n_values,
    "error_rate" = unlist(error_rates)
  )
  return(error_rate_data)
}

#' @title Generate a Plot of Error Rates
#'
#' @description Find how the predicted error rate for your methylation data
#'  changes as the read depth and percent methylation thresholds change
#'
#' @inheritParams bedmethyl_format
#' @param max_read_depth The maximum read depth threshold to view (integer).
#'  Some data sets have a single site with really high read depth, plotting
#'  all the way up to this point is a waste of computation and won't give you
#'  meaningful data (estimating an error rate using a single site). Defaults to
#'  1000.
#'
#' @return A plot (ggplot) of type `geom_point` showcasing how the estimated
#'  error rate changes with the minimum read depth considered
#'
#' @examples
#' create_error_rate_plot(methylation_data, 1000, 5)
#'
#' @export
create_error_rate_plot <- function(methylation_data,
                                   max_read_depth = 1000,
                                   percent_threshold = 5) {
  max_read_depth <- min(max(methylation_data[["read_depth"]]), max_read_depth)

  error_rate_data <- create_error_rate_data(
    methylation_data,
    max_read_depth,
    percent_threshold
  )
  error_rate_plot <-
    ggplot2::ggplot(
      error_rate_data,
      ggplot2::aes(x = n, y = error_rate)
    ) +
    ggplot2::geom_point(color = "black") +
    ggplot2::labs(
      x = "Smallest read depth considered",
      y = "Estimated error rate"
    ) +
    ggplot2::theme_bw()
  return(error_rate_plot)
}
