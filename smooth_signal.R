# This script will see how clustered the CpGs that are methylated are
# Ideally, methylation should be highly clustered around areas of high
# methylation, and less clustered around areas with singletons.

args <- commandArgs(trailingOnly = TRUE)

binarized_file_location <- args[1]
iterations <- as.numeric(args[2])
output_folder <- args[3]

library(magrittr)

binarized_data <- data.table::fread(
  binarized_file_location
)
binarized_data <- binarized_data |>
  dplyr::rename("chr" = V1, "start" = V2, "end" = V3, "iteration_0" = V4) |>
  dplyr::mutate(iteration_0 = as.numeric(iteration_0))

smooth_binarized_signal <- function(binarized_data, iteration) {
  reference_column <- paste0("iteration_", iteration - 1)
  new_column_name <- paste0("iteration_", iteration)
  binarized_data <- binarized_data %>%
    dplyr::mutate(
      leading = dplyr::lead(!!as.symbol(reference_column)),
      lagging = dplyr::lag(!!as.symbol(reference_column))
    ) %>%
    dplyr::mutate(
      new = rowMeans(
        .[, c(!!reference_column, "leading", "lagging")],
        na.rm = TRUE
      )
    ) %>%
    dplyr::select(-c("leading", "lagging")) %>%
    dplyr::rename(!!new_column_name := new)
  return(binarized_data)
}

for (i in 1:iterations) {
  binarized_data <- smooth_binarized_signal(binarized_data, i)
}

data.table::fwrite(binarized_data,
  file = paste0(output_folder, "/smoothed_data.txt"),
  col.names = FALSE,
  sep = " "
)


## ============ ##
##   PLOTTING   ##
## ============ ##

library(ggplot2)

binarized_data <- binarized_data |>
  dplyr::select(-c("end", "chr"))

binarized_data

line_plot <- binarized_data %>%
  tidyr::gather("key", "value", -start) %>%
  ggplot(aes(x = start, y = value, color = key)) +
  geom_line()

histogram_plot <- binarized_data %>%
  ggplot(aes(x = iteration_5)) +
  geom_histogram()

ggsave(paste0(output_folder, "/collated_lines.png"), line_plot)

ggsave(paste0(output_folder, "/histogram.png"), histogram_plot)
