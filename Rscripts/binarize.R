args <- commandArgs(trailingOnly = TRUE)
bin_counts_file <- args[1]
dense_output_file <- args[2]
sparse_output_file <- args[3]
bin_size <- as.numeric(args[4])
beta_threshold <- as.numeric(args[5])

## ======================== ##
##   DISTRIBUTION FITTING   ##
## ======================== ##

get_beta_parameters <- function(bin_counts) {
  fit <- fitdistrplus::fitdist(
    bin_counts[["density"]],
    "beta",
    method = "mle"
  )

  shape1 <- coef(fit)[["shape1"]]
  shape2 <- coef(fit)[["shape2"]]

  return(c(shape1, shape2))
}

## ==================== ##
##  BIN IDENTIFICATION  ##
## ==================== ##

remove_zero_bins <- function(bin_counts) {
  return(dplyr::filter(bin_counts, count > 0))
}

is_densely_methylated <-
  function(bin_count, shape1, shape2, threshold) {
    return(
      as.numeric( # Numeric is used as T/F is not as transferable as 0/1
        pbeta(bin_count, shape1, shape2, lower.tail = FALSE) < threshold
      )
    )
  }

## ======== ##
##   MAIN   ##
## ======== ##

bin_counts <- data.table::fread(bin_counts_file)
names(bin_counts) <- c("chr", "start", "end", "count")

bin_counts <- dplyr::mutate(bin_counts,
  "density" = count / bin_size
)

# We want to discern between sparsely and densely methylated bins. As such
# we remove any bins with zero signal as these bins will massively skew our
# beta distribution to the left.
beta_parameters <- get_beta_parameters(
  remove_zero_bins(bin_counts)
)
shape1 <- beta_parameters[[1]]
shape2 <- beta_parameters[[2]]

bin_counts <- bin_counts |>
  dplyr::mutate(
    "densely_methylated" = is_densely_methylated(
      count,
      shape1,
      shape2,
      beta_threshold
    )
  ) |>
  dplyr::mutate(
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
