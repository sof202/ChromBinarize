#' @title Estimate the Error Rate of Methylation Data
#'
#' @description This function estimates the error rate of methylation data
#'  using predicted unmethylated base pair positions. Any methylation signal
#'  in the input file is considered incorrect and contributes to the estimated
#'  error rate.
#'
#' @param reference_set_path A file path (string) to a file detailing
#'  CpG positions in the epigenome that do not exhibit methylation.
#'  The file is expected to have columns: "read depth at base pair position"
#'  and "percentage of such reads called as methylated", e.g.:
#'
#'      34	5.9
#'      30	5
#'      42	4.4
#'      51	5.9
#'      22	0
#'
#' @return An estimated error rate (numeric) expressed as a proportion (0 to 1).
#'
#' @details The input file should only contain positions with a low percentage
#'  of reads called as methylated. By default, ChromBinarize uses a reference
#'  set where the read depth is at least **500** and the percentage of reads
#'  being methylated is less than **5%**. It can also use only CpGs in CGIs,
#'  which are typically unmethylated regions of the genome.
#'
#'  **Warning**: If you are too stringent on unmethylated positions, the
#'  estimated error rate may approach **0**.
#'  This may not be helpful. Aim for a larger error rate that is also accurate.
#'
#' @examples
#' # Estimate error rate using a TSV file
#' estimate_error_rate("path/to/reference_set.tsv")
#'
#' # Estimate error rate using a CSV file
#' estimate_error_rate("path/to/reference_set.csv")
#'
#' @export
estimate_error_rate <- function(reference_set_path) {
  unmethylated_positions <- data.table::fread(
    reference_set_path,
    col.names = c("reads", "percent_methylated")
  )

  unmethylated_positions <-
    dplyr::mutate(
      unmethylated_positions,
      "incorrectly_methylated" = reads * percent_methylated / 100
    )

  incorrectly_methylated_total <-
    sum(unmethylated_positions[["incorrectly_methylated"]])
  total_reads <- sum(unmethylated_positions[["reads"]])
  error_rate <-
    incorrectly_methylated_total / total_reads

  return(error_rate)
}
