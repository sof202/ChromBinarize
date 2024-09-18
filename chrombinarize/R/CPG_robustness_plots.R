calculate_nearby_methylation <- function(methylation_data,
                                         current_start,
                                         min_distance,
                                         max_distance) {
  max_lower_bound <- current_start - max_distance
  max_upper_bound <- current_start + max_distance
  min_lower_bound <- current_start - min_distance
  min_upper_bound <- current_start + min_distance

  average_count <-
    dplyr::filter(methylation_data, start != current_start) |>
    dplyr::filter(start <= max_upper_bound & start >= min_upper_bound) |>
    dplyr::filter(start <= min_lower_bound & start >= max_lower_bound) |>
    dplyr::summarise("average" = mean(percent_methylation, na.rm = TRUE)) |>
    dplyr::pull()

  return(average_count)
}

create_cpg_robustness_data <- function(methylation_data,
                                       min_distance,
                                       max_distance) {
  cpg_robustness_data <- methylation_data |>
    dplyr::rowwise() |>
    dplyr::mutate("surrounding_methylation" = calculate_nearby_methylation(
      methylation_data,
      start,
      min_distance,
      max_distance
    )) |>
    dplyr::filter(!is.na(surrounding_methylation)) |>
    data.table::as.data.table()

  return(cpg_robustness_data)
}

create_cpg_robustness_plot <- function(methylation_data,
                                       min_distance,
                                       max_distance) {
  methylation_data <- create_cpg_robustness_data(
    methylation_data,
    min_distance,
    max_distance
  )
  methylation_correlation <- cor(
    methylation_data[["percent_methylation"]],
    methylation_data[["surrounding_methylation"]]
  )

  cpg_robustness_plot <-
    ggplot2::ggplot(
      methylation_data,
      ggplot2::aes(x = percent_methylation, y = surrounding_methylation)
    ) +
    ggplot2::geom_bin2d(bins = 100) +
    ggplot2::scale_fill_continuous(type = "viridis") +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "black") +
    ggplot2::labs(
      x = "Percent methylation",
      y = "Average surrounding percent methylation"
    ) +
    ggplot2::scale_fill_gradientn(colors = c("grey", "red", "yellow")) +
    ggplot2::annotate(
      "text",
      x = 95,
      y = 5,
      label = paste0("R^2: ", methylation_correlation)
    ) +
    ggplot2::theme_bw()

  return(cpg_robustness_plot)
}
