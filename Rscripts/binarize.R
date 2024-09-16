args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
bin_counts_file <- args[2]
dense_output_file <- args[3]
sparse_output_file <- args[4]
bin_size <- as.numeric(args[5])
beta_threshold <- as.numeric(args[6])

renv::load(renv_environment)

bin_counts <- chrombinarize::determine_dense_bins(
  bin_counts_file,
  beta_threshold
)

bin_counts <- dplyr::mutate(
  bin_counts,
  "methylation_present" = as.numeric(count > 0)
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
