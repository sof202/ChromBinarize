#' @title Add Absolute Change to BEDMethyl Data
#'
#' @description Adds two columns:
#' - "absolute change in read depth"
#' - "absolute change in percent methylation"
#'
#' to a `comparison_bedmethyl` data table.
#' This is primarily for the purpose of using the
#' `create_read_depth_plot()` and `create_percent_comparison_plot()` functions.
#'
#' @inheritParams create_read_depth_plot
#'
#' @return A modified version of the `comparison_bedmethyl` data.table with
#'    additional columns:
#' - absolute_change_percent_methylation: absolute change in percent
#'    methylation (numeric)
#' - absolute_change_read_depth: absolute change in read depth (integer)
#'
#' The returned data.table retains the following columns:
#' - chr: chromosome name (string)
#' - start: starting base pair position (integer)
#' - end: ending base pair position (integer)
#' - mark_name: "m" for 5mC and "h" for 5hmC (string)
#' - ONT_read_depth: read depth for ONT (integer)
#' - ONT_percent_methylation: percentage of methylated reads in ONT
#'    (numeric)
#' - BS_read_depth: read depth for BS (integer)
#' - BS_percent_methylation: percentage of methylated reads in BS (numeric)
add_absolute_change_columns <- function(comparison_bedmethyl) {
  comparison_bedmethyl <- dplyr::mutate(
    comparison_bedmethyl,
    "absolute_change_percent_methylation" = abs(
      ONT_percent_methylation - BS_percent_methylation
    ),
    "absolute_change_read_depth" = abs(ONT_N - BS_N)
  )
  return(comparison_bedmethyl)
}



#' @title Generate a Histogram for the Change in Read Depth
#'
#' @description Using comparative methylation data from BS-Seq and ONT, a
#'  histogram is produced to showcase the distribution of the absolute change
#'  in read depth seen between each dataset.
#'
#' @param comparison_bedmethyl A data.table with the following columns:
#' - chr: chromosome name (string)
#' - start: starting base pair position (integer)
#' - end: ending base pair position (integer)
#' - mark_name: "m" for 5mC and "h" for 5hmC (string)
#' - read_depth: read depth (integer)
#' - percent_methylation: percentage of reads observed to be methylated
#'    (numeric)
#' @param mark The name of the mark to inspect ("m" or "h")
#' @param read_depth_filter A filter to use on the read depth in your bedmethyl
#'  file. This is in place to reduce random fluctuations (common in low count
#'  data)
#'
#' @return A plot (ggplot) of type `geom_histogram` showing the distribution
#'  of the absolute change in percent methylation between the datasets
#'
#' @examples
#' # Read in comparison_bedmethyl file
#' bedmethyl <- read_comparison_bedmethyl("path/to/bedmethyl.bed")
#'
#' # Create histogram
#' create_read_depth_plot(bedmethyl)
#'
#' @export
create_read_depth_plot <- function(comparison_bedmethyl,
                                   mark = c("m", "h"),
                                   read_depth_filter = 30) {
  mark <- match.arg(mark)

  comparison_bedmethyl <- add_absolute_change_columns(comparison_bedmethyl)
  comparison_bedmethyl <- dplyr::filter(
    comparison_bedmethyl,
    ONT_N >= !!read_depth_filter,
    BS_N >= !!read_depth_filter,
    mark == !!mark
  )

  # read depth can get very large if lots of samples are merged. If this is the
  # case, then the histogram becomes very difficult to interpret. We remove
  # highly covered sites for this reason.
  comparison_bedmethyl <- dplyr::filter(
    comparison_bedmethyl,
    absolute_change_read_depth < 1000
  )

  read_depth_plot <-
    ggplot2::ggplot(
      comparison_bedmethyl,
      ggplot2::aes(x = absolute_change_read_depth)
    ) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "absolute change in read depth")
  return(read_depth_plot)
}

#' @title Generate a Histogram for the Change in Percent Methylation
#'
#' @description Using comparative methylation data from BS-Seq and ONT, a
#'  histogram is produced to showcase the distribution of the absolute change
#'  in percentage methylation seen between each dataset.
#'
#' @inheritParams create_read_depth_plot
#'
#' @return A plot (ggplot) of type `geom_histogram` showing the distribution
#'  of the absolute change in percent methylation between the datasets
#'
#' @inherit create_read_depth_plot details
#'
#' @examples
#' # Read in comparison_bedmethyl file
#' bedmethyl <- read_comparison_bedmethyl("path/to/bedmethyl.bed")
#'
#' # Create histogram
#' create_read_depth_plot(bedmethyl)
#'
#' @export
create_percent_comparison_plot <- function(comparison_bedmethyl,
                                           mark = c("m", "h"),
                                           read_depth_filter = 30) {
  mark <- match.arg(mark)

  comparison_bedmethyl <- add_absolute_change_columns(comparison_bedmethyl)
  comparison_bedmethyl <- dplyr::filter(
    comparison_bedmethyl,
    ONT_N >= !!read_depth_filter,
    BS_N >= !!read_depth_filter,
    mark == !!mark
  )

  percent_comparison_plot <-
    ggplot2::ggplot(
      comparison_bedmethyl,
      ggplot2::aes(x = absolute_change_percent_methylation)
    ) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "absolute change in percent methylation")
  return(percent_comparison_plot)
}


#' @title Generate a Density Heatmap Correlation Plot
#'
#' @description Using comparitive methylation data from BS-Seq and ONT, a
#'  density heatmap is produced between the percent methylation seen at each
#'  CpG site for both data types.
#'
#' @inheritParams create_read_depth_plot
#'
#' @return A plot (ggplot) of type `geom_bin2d` showing how the percent
#' methylation in BS and ONT compare
#'
#' @details
#'  The most ideal plot would be one where all colouration is displayed on the
#'  main diagonal (datasets agree with one another). This is assuming of course
#'  that both datasets are from the same study/sample etc.
#' @export
create_correlation_plot <- function(comparison_bedmethyl,
                                    mark = c("m", "h"),
                                    read_depth_filter = 30) {
  mark <- match.arg(mark)

  comparison_bedmethyl <- dplyr::filter(
    comparison_bedmethyl,
    ONT_N >= !!read_depth_filter,
    BS_N >= !!read_depth_filter,
    mark == !!mark
  )
  # For better interpretability, a log scale is used for colouration
  breaks <- c(1, 10, 100, 1000, 10000)

  methylation_correlation_plot <-
    ggplot2::ggplot(
      comparison_bedmethyl,
      ggplot2::aes(x = ONT_percent_methylation, y = BS_percent_methylation)
    ) +
    ggplot2::geom_bin2d(bins = 100) +
    ggplot2::scale_fill_gradient(
      name = "count",
      trans = "log",
      breaks = breaks,
      labels = breaks
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "percent methylation in ONT",
      y = "percent methylation in WGBS"
    )
  return(methylation_correlation_plot)
}
