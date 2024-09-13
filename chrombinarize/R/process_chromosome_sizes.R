#' @title Read and process a chromosome lengths file
#'
#' @description Convert a file detailing the sizes of each chromosome
#'  into a datatable with the correct column headings
#'
#' @param chromosome_sizes_file
#'  A file path (string) to a text file detailing the length of each chromosome.
#'  The file should contain two columns: chromosome name and chromosome length.
#'
#' @return A data.table with columns "chromsome" and "size" for each
#'  chromosome given in the input file
#'
#' @examples
#' process_chromosome_sizes("path/to/chromosome_sizes.txt")
#'
#' @export
process_chromosome_sizes <- function(chromosome_sizes_file) {
  chromosome_sizes_table <- data.table::fread(chromosome_sizes_file)
  names(chromosome_sizes_table) <- c("chromosome", "size")
  chromosome_sizes <- chromosome_sizes_table$size
  names(chromosome_sizes) <- chromosome_sizes_table$chromosome
  return(chromosome_sizes)
}
