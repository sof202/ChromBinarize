#' @title Estimates Shape Parameters of Beta Distribution
#'
#' @description Using the fitdistrplus package, the parameters for the beta
#'  distribution (lovingly called shape1 and shape2) are estimated using
#'  the reported density of reads in each window.
#'
#' @param bin_densities
#'  A numeric vector, detailing the collection of read densities across all
#'  genomic regions considered.
#'
#' @return A numeric vector containing estimated parameters shape1 and shape2
#'
#' @examples
#'  bin_densities <- c(0.001, 0.09, 0.1, 0.001, 0.002)
#'  get_beta_parameters(bin_densities)
#'  [1] 0.3834632 9.7126634
get_beta_parameters <- function(bin_densities) {
  fit <- fitdistrplus::fitdist(
    bin_densities,
    "beta",
    method = "mle"
  )

  shape1 <- coef(fit)[["shape1"]]
  shape2 <- coef(fit)[["shape2"]]

  return(c(shape1, shape2))
}

#' @title Determine if a Bin is Densely Methylated
#'
#' @description Using the beta distribution and the given bin density the upper
#'  tail of the beta distibution is calculated and compared against the given
#'  threshold.
#'
#' @param bin_density
#'  The density of the bin being inspected (numeric)
#' @param shape1, shape2
#'  Beta distriubution parameters (numeric)
#' @param threshold
#'  The value required for the p-value from the beta distribution to be
#'  considered significant
#'
#' @return A Boolean value (represented by a 0 or 1) where 1 indicates that
#'  the bin density provided is significantly above what is expected (and 0
#'  is the opposite of course)
#'
#' @examples
#'  # A low density likely returns 0
#'  is_densely_methylated(0.001, 0.3, 9.5, 0.001)
#'  0
#'
#'  # A high density likely returns 1
#'  is_densely_methylated(0.4, 0.3, 9.5, 0.001)
#'  1
is_densely_methylated <-
  function(bin_density, shape1, shape2, threshold = 0.001) {
    return(as.numeric(
      pbeta(bin_density, shape1, shape2, lower.tail = FALSE) < threshold
    ))
  }


#' @title Converts Bin Counts to Data Table
#'
#' @description Reads in a BED file detailing the number of reads in each
#'  region as a data.table with the necessary column names for
#'  `determine_dense_bins()`
#'
#' @inheritParams determine_dense_bins
#'
#' @return A data.table with the columns: "chr" (chromosome, string)
#'  "start" (start position, integer), "end" (end position, integer) and
#'  "count" (number of reads in bin, integer)
#'
#' @examples
#'  read_bin_counts_file("path/to/bin_counts.bed")
#'
#'        chr start   end count
#'     <char> <int> <int> <int>
#'  1:   chr1     0   200     3
#'  2:   chr1   200   400     4
#'  3:   chr1   400   600     0
#'  4:   chr1   600   800     1
#'  5:   chr1   800  1000   199
#'  ...
read_bin_counts_file <- function(bin_counts_file) {
  if (!file.exists(bin_counts_file)) stop("ERROR: File does not exist.")

  bin_counts <- suppressWarnings(data.table::fread(
    bin_counts_file,
    col.names = c("chr", "start", "end", "count"),
    colClasses = c("character", "integer", "integer", "integer")
  ))

  verify_column_class(
    bin_counts[["chr"]],
    is.character,
    "The first column (chromosome name) must be a string"
  )
  verify_column_class(
    bin_counts[["start"]],
    is.integer,
    "The second column (start) must be an integer"
  )
  verify_column_class(
    bin_counts[["end"]],
    is.integer,
    "The third column (end) must be an integer"
  )
  verify_column_class(
    bin_counts[["count"]],
    is.integer,
    "The fourth column (count) must be an integer"
  )
  return(bin_counts)
}


#' @title Determines Densely Methylated Regions
#'
#' @description Using the beta distribution, this finds bins with a higher
#'  than expected density of reads.
#'
#' @param bin_counts_file
#'  A file path (string) to a file detailing the number of reads in each bin
#'  This file must have columns: "chromosome", "start position", "end position"
#'  and "number of reads in region". i.e. a BED3+1 file where the final column
#'  is the number of reads in the region
#' @param beta_threshold
#'  A numeric value that determines the significance threshold required for
#'  a bin to be densely populated with reads. Defaults to 0.001
#'
#' @return A data.table detailing which bins in the input file are densely
#'  populated with methylated sites. The data.table has columns:
#'  "chromosome", "start position", "end position", "number of reads in region"
#'  and "is densely methylated" (Boolean value)
#'
#' @examples
#' determine_dense_bins("path/to/bin_counts.bed", 200, 0.0001)
#'
#'        chr start   end count density densely_methylated
#'     <char> <int> <int> <int>   <num>              <num>
#'  1:   chr1     0   200     3   0.015                  0
#'  2:   chr1   200   400     4   0.020                  0
#'  3:   chr1   400   600     0   0.000                  0
#'  4:   chr1   600   800     1   0.005                  0
#'  5:   chr1   800  1000   199   0.995                  1
#'  ...
#'
#' @export
determine_dense_bins <- function(bin_counts_file, beta_threshold = 0.001) {
  bin_counts <- read_bin_counts_file(bin_counts_file)
  bin_size <- bin_counts[["end"]][1] - bin_counts[["start"]][1]
  bin_counts <- dplyr::mutate(bin_counts,
    "density" = count / bin_size
  )

  # We want to discern between sparsely and densely methylated bins. As such
  # we remove any bins with zero signal as these bins will massively skew our
  # beta distribution to the left. Most of the genome is not methylated
  beta_parameters <- get_beta_parameters(
    dplyr::filter(bin_counts, count > 0) |> dplyr::pull(density)
  )
  shape1 <- beta_parameters[[1]]
  shape2 <- beta_parameters[[2]]

  bin_counts <- bin_counts |>
    dplyr::mutate(
      "densely_methylated" = is_densely_methylated(
        density,
        shape1,
        shape2,
        beta_threshold
      )
    )

  return(bin_counts)
}
