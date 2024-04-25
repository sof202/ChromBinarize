args <- commandArgs(trailingOnly = TRUE)
bin_counts_file <- args[1]
output_directory <- args[2]

remove_zero_bins <- function(bin_counts) {
  return(dplyr::filter(bin_counts, count > 0))
}

get_average_bin_count <- function(bin_counts) {
  # 0 count bins dominate the methylome, we want to discern between high
  # density and low density only
  bin_counts <- remove_zero_bins(bin_counts)
  return(mean(bin_counts$count))
}

is_significantly_methylated <-
  function(bin_count, average_bin_count, threshold = 0.0001) {
    p_value <- ppois(bin_count, average_bin_count, lower.tail = FALSE)
    if (p_value < threshold) {
      return(1)
    }
    return(0)
  }


bin_counts <- data.table::fread(bin_counts_file)
names(bin_counts) <- c("chr", "start", "end", "count")
average_bin_count <- get_average_bin_count(bin_counts)

bin_counts <- bin_counts |>
  dplyr::mutate(
    "significantly_methylated" = is_significantly_methylated(
      count,
      average_bin_count
    )
  ) |>
  dplyr::select(
    significantly_methylated
  )

data.table::fwrite(
  bin_counts,
  file = output_directory,
  append = TRUE
)
