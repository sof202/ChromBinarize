args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
chromosome_sizes_file <- args[2]
bin_size <- as.numeric(args[3])
output_directory <- args[4]

renv::load(renv_environment)

process_chromosome_sizes <- function(chromosome_sizes_file) {
  chromosome_sizes_table <- data.table::fread(chromosome_sizes_file)
  names(chromosome_sizes_table) <- c("chromosome", "size")
  chromosome_sizes <- chromosome_sizes_table$size
  names(chromosome_sizes) <- chromosome_sizes_table$chromosome
  return(chromosome_sizes)
}

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

chromosome_sizes <- process_chromosome_sizes(chromosome_sizes_file)

blank_bed_files <-
  create_blank_bed_files(
    chromosome_sizes = chromosome_sizes,
    bin_size = bin_size
  )

invisible(lapply(names(blank_bed_files), function(file) {
  do.call(data.table::fwrite, c(list(
    x = blank_bed_files[[file]],
    file = paste0(output_directory, "/chromosome", file, ".bed"),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )))
}))
