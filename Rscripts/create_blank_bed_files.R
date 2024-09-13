args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
chromosome_sizes_file <- args[2]
bin_size <- as.numeric(args[3])
output_directory <- args[4]

renv::load(renv_environment)

chromosome_sizes <-
  chrombinarize::process_chromosome_sizes(chromosome_sizes_file)

blank_bed_files <-
  chrombinarize::create_blank_bed_files(
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
