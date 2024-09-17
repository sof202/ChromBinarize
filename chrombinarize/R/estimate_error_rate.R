#' @title Estimate the Error Rate of Methylation Data
#'
#' @description This function estimates the error rate of methylation data
#'  using predicted unmethylated base pair positions. Any methylation signal
#'  in the input file is considered incorrect and contributes to the estimated
#'  error rate.
#'
#' @param methylation_data
#' A data.table with columns:
#' - chr: chromosome name (string)
#' - start: starting base pair position (integer)
#' - end: ending base pair position (integer)
#' - mark_name: "m" for 5mC and "h" for 5hmC (string)
#' - read_depth: read depth for the site (integer)
#' - strand: "+", "-" or "." (character)
#' - percent_methylation: percentage of reads reported as methylated (numeric)
#' @param read_depth_threshold The minimum read depth to consider (integer).
#'  Higher values are recommended. Defaults to 30.
#' @param percent_threshold The maximum percent methylation to consider
#'  (numeric). Lower values are recommended. Defaults to 5.
#'
#' @return An estimated error rate (numeric) expressed as a probability
#'  (0 to 1).
#'
#' @details The input file should only contain positions with a low percentage
#'  of reads called as methylated. By default, ChromBinarize uses a reference
#'  set where the read depth is at least **500** and the percentage of reads
#'  being methylated is less than **5%**. It can also use only CpGs in CGIs,
#'  which are typically unmethylated regions of the genome.
#'
#'  **Warning**: If you are too stringent on unmethylated positions, the
#'  estimated error rate will approach 0.
#'  This may not be helpful. Aim for a larger error rate that is also accurate.
#'
#' @examples
#' estimate_error_rate(bedmethyl_data, 40, 5)
#' estimate_error_rate(bedmethyl_data, 500, 2)
#' estimate_error_rate(bedmethyl_data)
#'
#' # Bad usage
#' estimate_error_rate(bedmethyl_data, 1, 5)
#' estimate_error_rate(bedmethyl_data, 30, 0.1)
#'
#' @export
estimate_error_rate <- function(methylation_data,
                                read_depth_threshold = 30,
                                percent_threshold = 5) {
  unmethylated_positions <-
    dplyr::filter(
      methylation_data,
      read_depth >= !!read_depth_threshold,
      percent_methylation <= !!percent_threshold
    ) |>
    dplyr::mutate(
      "incorrectly_methylated" = read_depth * percent_methylation / 100
    )

  incorrectly_methylated_total <-
    sum(unmethylated_positions[["incorrectly_methylated"]])
  total_reads <- sum(unmethylated_positions[["read_depth"]])
  estimated_error_rate <-
    incorrectly_methylated_total / total_reads

  return(estimated_error_rate)
}
