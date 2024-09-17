args <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
folder <- args[2]
reference_set <- args[3]
input_bed_file <- args[4]
output_file_name <- args[5]

renv::load(renv_environment)

reference_set_path <- file.path(folder, reference_set)
binomial_p <- chrombinarize::estimate_error_rate(reference_set_path)

bed_file_location <- file.path(folder, input_bed_file)
methylation_data <- chrombinarize::read_bedmethyl(bed_file_location)

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
