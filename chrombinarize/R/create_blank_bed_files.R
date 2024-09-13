#' @title Create genomic windows for a chromosome
#'
#' @description Given a chromosome and its size this creates a data.table
#'  that details regions of a fixed length spanning the chromosome
#'
#' @param chromsome The name of the chromsome (string)
#' @param chromosome_sizes a data.table detailing the length of each chromosome
#'  with columns "chromosome name" and "chromosome size"
#' @param bin_size The desired size of each region (integer)
#'
#' @return A data.table with columns "chromosome",
#' "starting base pair position" and "end base pair position"
#'
#' @examples
#'  create_bins("22", chromosome_sizes, 200)
#'
#'  Generates a data.table of the form:
#'  chr22	0	200
#'  chr22	200	400
#'  chr22	400	600
#'  ...
create_bins <- function(chromosome, chromosome_sizes, bin_size) {
  chromosome_length <- chromosome_sizes[[paste0("chr", chromosome)]]
  bin_starts <- seq(0, chromosome_length, bin_size)
  bins <- data.table::data.table("start" = bin_starts)
  chromsome_name <- paste0("chr", chromosome)
  bins <- bins |>
    dplyr::mutate(
      "chr" = chromsome_name,
      "end" = start + bin_size
    ) |>
    dplyr::select(chr, start, end)
  return(bins)
}

#' @title A wrapper for `create_bins()`
#'
#' @description Using a desired bin size, creates data.tables in the form of
#'  a bed file for chromosomes 1 to X. Each row is the next subsequent window
#'  equal to the bin size.
#'
#' @param chromosome_sizes a data.table detailing the length of each chromosome
#'  with columns "chromosome name" and "chromosome size"
#' @param bin_size The desired size of each region (integer)
#'
#' @return A data.table with columns "chromosome",
#' "starting base pair position" and "end base pair position"
#'
#' @examples
#'  create_bins("22", chromosome_sizes, 200)
#'
#'  Generates a data.table of the form:
#'  chr22	0	200
#'  chr22	200	400
#'  chr22	400	600
#'  ...
#'
#' @export
create_blank_bed_files <- function(chromosome_sizes, bin_size) {
  chromosomes <- c(seq(1, 22), "X")
  blank_bed_files <- lapply(chromosomes, function(chromosome) {
    do.call(create_bins, c(list(
      chromosome = chromosome,
      chromosome_sizes = chromosome_sizes,
      bin_size = bin_size
    )))
  })
  names(blank_bed_files) <- chromosomes
  return(blank_bed_files)
}
