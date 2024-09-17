args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
bed_file_location <- args[2]
read_depth_threshold <- as.numeric(args[3])
percent_threshold <- as.numeric(args[4])
output_file_name <- args[5]

renv::load(renv_environment)

methylation_data <- chrombinarize::read_bedmethyl(bed_file_location)
binomial_p <- chrombinarize::estimate_error_rate(
  methylation_data,
  read_depth_threshold,
  percent_threshold
)

methylation_data <- methylation_data |>
  dplyr::mutate(
    read_depth = as.numeric(read_depth),
    percent_methylation = as.numeric(percent_methylation),
    number_of_methylated_reads = read_depth * percent_methylation / 100,
    "likelihood of methylated reads being erroneous" = pbinom(
      number_of_methylated_reads,
      read_depth,
      binomial_p,
      lower.tail = FALSE
    )
  ) |>
  dplyr::select(
    -number_of_methylated_reads
  )

data.table::fwrite(methylation_data,
  file = file.path(folder, output_file_name),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
