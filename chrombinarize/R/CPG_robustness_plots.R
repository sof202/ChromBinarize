#' @title Calculate the Average Methylation for Nearby CpGs
#'
#' @description Given distances of what defines two CpGs being 'nearby' this
#'  function calculates the average percent methylation for such CpGs.
#'
#' @inheritParams create_cpg_robustness_plot
#' @param chromosome The chromosome that the CpG lies on (string)
#' @param current_start The start position of the CpG in question (integer)
#'
#' @return The average percent methylation for 'nearby' CpGs to the input CpG
#'  (numeric)
#'
#' @examples
#' # Given the following input data:
#'       chr start   end   name read_depth strand percent_methylation
#'    <char> <num> <num> <char>      <num> <char>               <num>
#' 1:   chr1    10    11      m         21      .                   4
#' 2:   chr1    15    16      m        123      .                   33
#' 3:   chr1    20    21      m         12      .                   21
#' 4:   chr1    25    26      m         54      .                   50
#'
#' calculate_nearby_methylation(methylation_data, "chr1", 10, 0, 5)
#'  33
#' calculate_nearby_methylation(methylation_data, "chr1", 20, 0, 10)
#'  25
#' calculate_nearby_methylation(methylation_data, "chr1", 20, 6, 10)
#'  27
calculate_nearby_methylation <- function(methylation_data,
                                         chromosome,
                                         current_start,
                                         min_distance = 0,
                                         max_distance = 30) {
  max_lower_bound <- current_start - max_distance
  max_upper_bound <- current_start + max_distance
  min_lower_bound <- current_start - min_distance
  min_upper_bound <- current_start + min_distance

  nearby_counts <- methylation_data |>
    dplyr::filter(chr == chromosome) |>
    dplyr::filter(start != current_start) |>
    dplyr::filter(
      (start <= max_upper_bound & start >= min_upper_bound) |
        (start <= min_lower_bound & start >= max_lower_bound)
    ) |>
    dplyr::pull(percent_methylation)

  return(mean(nearby_counts, na.rm = TRUE))
}

#' @title Create Data for `create_cpg_robustness_plot()`
#'
#' @description Adds the column: average methylation of surrounding CpGs to the
#'  inputted methylation data.
#'
#' @inheritParams create_cpg_robustness_plot
#'
#' @return A modified version of the input data.table with the added column:
#' - surrounding_methylation: Average methylation signal of 'nearby' CpGs
#'  (numeric)
#'
#' @examples
#' create_cpg_robustness_data(bedmethyl_data, 0, 5)
#'
#'       chr start   end   name read_depth strand percent_methylation
#'    <char> <num> <num> <char>      <num> <char>               <num>
#' 1:   chr1    10    11      m         21      .                   2
#' 2:   chr1    15    16      m        123      .                   3
#' 3:   chr1    20    21      m         12      .                   4
#' 4:   chr1    25    26      m         54      .                   5
#'    surrounding_methylation
#'                      <num>
#' 1:                       4
#' 2:                       5
#' 3:                       2
#' 4:                       3
create_cpg_robustness_data <- function(methylation_data,
                                       min_distance = 0,
                                       max_distance = 30) {
  cpg_robustness_data <- methylation_data |>
    dplyr::rowwise() |>
    dplyr::mutate("surrounding_methylation" = calculate_nearby_methylation(
      methylation_data,
      chr,
      start,
      min_distance,
      max_distance
    )) |>
    dplyr::filter(!is.na(surrounding_methylation)) |>
    data.table::as.data.table()

  return(cpg_robustness_data)
}


#' @title Create a Plot Showcasing Methylation Stability
#'
#' @description Generates a density plot that shows how similar the methylation
#'   signal is for surrounding CpG sites for each CpG.
#'
#' @inheritParams bedmethyl_format
#' @param min_distance The smallest distance a CpG can be away from another to
#'  be called nearby. Integer valued, defaults to 0.
#' @param max_distance The largest distance a CpG can be away from another to
#'  be called nearby. Integer valued, defaults to 30.
#'
#' @return A plot (ggplot) of type `geom_bin2d` showing how the methylation
#'  signal of CpGs correlate with 'nearby' CpGs.
#'
#' @details A plot where most points lie on the main diagonal is one where
#'  nearby CpGs are generally very similar (in terms of methylation status)
#'
#'  It might seem odd that min_distance is a parameter here. This is mainly
#'  here to give the user a little more control. If you want to see how
#'  dissimilar far away CpGs are from each other, you would want to increase
#'  this value to remove closer CpGs.
#'
#' @export
create_cpg_robustness_plot <- function(methylation_data,
                                       min_distance = 0,
                                       max_distance = 30) {
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
