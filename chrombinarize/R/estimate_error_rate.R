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
