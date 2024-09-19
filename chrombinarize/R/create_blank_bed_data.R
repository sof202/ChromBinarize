#' @title Read and Process a Chromosome Lengths File
#'
#' @description Convert a file detailing the sizes of each chromosome
#'  into a vector to be used by `create_blank_bed_data()`
#'
#' @inheritParams create_blank_bed_data
#'
#' @details It is assumed that the input file is coming from ChromHMM's
#'  included files. However, if this is not the case, you will need to ensure
#'  that chromosomes are given in the form "chr22", "chrX" etc. (chr prefixed)
#'
#' @return An integer vector with values equal to the given chromosome lengths.
#' Each entry has a name equal to the chromosome name given in the input file
#'
#' @examples
#' process_chromosome_sizes("path/to/chromosome_sizes.txt")
#'
#'      chr1                  chr2                  chr3
#' 249250621             243199373             198022430
#'      chr4                  chr5                  chr6
#' 191154276             180915260             171115067
#'      chr7                  chrX                  chr8
#' 159138663             155270560             146364022
#' ...
process_chromosome_sizes <- function(chromosome_sizes_file) {
  chromosome_sizes_table <- data.table::fread(
    chromosome_sizes_file,
    col.names = c("chromosome", "chromosome_size"),
    colClasses = c("character", "integer")
  )

  verify_column_class(
    chromosome_sizes_table[["chromosome"]],
    is.character,
    "The first column (chromosome name) must be a string"
  )
  verify_column_class(
    chromosome_sizes_table[["size"]],
    is.integer,
    "The second column (chromosome size) must be a string"
  )

  chromosome_sizes <- chromosome_sizes_table[["chromosome_size"]]
  names(chromosome_sizes) <- chromosome_sizes_table[["chromosome"]]
  return(chromosome_sizes)
}

#' @title Create Genomic Windows for a Chromosome
#'
#' @description Given a chromosome and its size this creates a data.table
#'  that details regions of a fixed length spanning the chromosome
#'
#' @param chromosome_name The name of the chromsome (string).
#'  This must be of the form "chrxyz" as ChromHMM expects this form.
#' @param chromosome_length The length of the given chromosome (integer)
#' @inheritParams create_blank_bed_data
#'
#' @return A data.table with columns "chromosome",
#' "starting base pair position" and "end base pair position"
#'
#' @examples
#'  create_bins("chr22", "1000000", 200)
#'
#'  chr22	0	200
#'  chr22	200	400
#'  chr22	400	600
#'  ...
#'  chr22	999800	1000000
create_bins <- function(chromosome_name, chromosome_length, bin_size) {
  if (!startsWith(chromosome_name, "chr")) {
    stop("chromosome name must start with the string 'chr'")
  }
  bin_starts <- seq(0, chromosome_length, bin_size)
  bins <- data.table::data.table("start" = bin_starts)
  bins <- bins |>
    dplyr::mutate(
      "chr" = chromosome_name,
      "end" = start + bin_size
    ) |>
    dplyr::select(chr, start, end)
  return(bins)
}

#' @title A Wrapper for `create_bins()`
#'
#' @description Using the given bin size, this creates data.tables for each
#'  chromosome requested in the form of a bed file. Each row gives the
#'  coordinates of a genomic window each equal to the bin size.
#'
#' @param chromosome_sizes_file
#'  A file path (string) to a text file detailing the length of each chromosome.
#'  The file should contain two columns:
#' - chromosome_name: A string detailing the chromosome name (ex: "chr22")
#' - chromosome_size: An integer detailing the length of the chromosome
#' @param bin_size The desired size of each region (integer). Defaults to 200.
#' @param chromosomes A character vector of all chromosomes to create bins for.
#'  Defaults to 1,2,3,...,22,X
#'
#' @return A list of data.tables with columns "chromosome",
#'  "starting base pair position" and "end base pair position". Each list item
#'  is given the corresponding chromosome index as its name.
#'
#' @examples
#'  create_blank_bed_data("path/to/chromosome_sizes.txt", 100, seq(1,22))
#'
#'  $`1`
#'             chr     start       end
#'          <char>     <num>     <num>
#'       1:   chr1         0       200
#'       3:   chr1       400       600
#'       4:   chr1       600       800
#'       5:   chr1       800      1000
#'      ---
#' 1246250:   chr1 249249800 249250000
#' 1246251:   chr1 249250000 249250200
#' 1246252:   chr1 249250200 249250400
#' 1246253:   chr1 249250400 249250600
#' 1246254:   chr1 249250600 249250800
#' ...
#'  $`22`
#'            chr    start      end
#'         <char>    <num>    <num>
#'      1:  chr22        0      200
#'      2:  chr22      200      400
#'      3:  chr22      400      600
#'      4:  chr22      600      800
#'      5:  chr22      800     1000
#'     ---
#' 256519:  chr22 51303600 51303800
#' 256520:  chr22 51303800 51304000
#' 256521:  chr22 51304000 51304200
#' 256522:  chr22 51304200 51304400
#' 256523:  chr22 51304400 51304600
#'
#' @export
create_blank_bed_data <- function(chromosome_sizes_file,
                                  bin_size = 200,
                                  chromosomes = c(seq(1, 22), "X")) {
  chromosome_sizes <- process_chromosome_sizes(chromosome_sizes_file)

  blank_bed_data <- lapply(chromosomes, function(chromosome) {
    chromosome_name <- paste0("chr", chromosome)
    chromosome_length <- chromosome_sizes[[paste0("chr", chromosome)]]
    do.call(create_bins, c(list(
      chromosome_name = chromosome_name,
      chromosome_length = chromosome_length,
      bin_size = bin_size
    )))
  })
  names(blank_bed_data) <- chromosomes
  return(blank_bed_data)
}
