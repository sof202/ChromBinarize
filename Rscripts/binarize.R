args <- commandArgs(trailingOnly = TRUE)
bin_counts_file <- args[1]
dense_output_file <- args[2]
sparse_output_file <- args[3]

## ==================== ##
##  BIN IDENTIFICATION  ##
## ==================== ##
remove_zero_bins <- function(bin_counts) {
  return(dplyr::filter(bin_counts, count > 0))
}

get_average_bin_count <- function(bin_counts) {
  # 0 count bins dominate the methylome, we want to discern between high
  # density and low density only
  bin_counts <- remove_zero_bins(bin_counts)
  return(mean(bin_counts$count))
}

is_densely_methylated <-
  function(bin_count, average_bin_count, threshold = 0.0001) {
    return(
      as.numeric( # Numeric is used as T/F is not as transferable as 0/1
        ppois(bin_count, average_bin_count, lower.tail = FALSE) < threshold
      )
    )
  }


## ===================== ##
##   DATA MANIPULATION   ##
## ===================== ##
bin_counts <- data.table::fread(bin_counts_file)
names(bin_counts) <- c("chr", "start", "end", "count")
average_bin_count <- get_average_bin_count(bin_counts)

bin_counts <- bin_counts |>
  dplyr::mutate(
    "densely_methylated" = is_densely_methylated(count, average_bin_count)
  ) |>
  dplyr::mutate(
    "methylation_present" = count > 0
  )

densely_methylated_bins <- dplyr::select(bin_counts, densely_methylated)

sparsely_methylated_bins <- bin_counts |>
  dplyr::filter(!densely_methylated) |>
  dplyr::select(methylation_present)

## =========== ##
##   OUTPUTS   ##
## =========== ##
data.table::fwrite(
  densely_methylated_bins,
  file = dense_output_file,
  append = TRUE
)

data.table::fwrite(
  sparsely_methylated_bins,
  file = sparse_output_file,
  append = TRUE
)
