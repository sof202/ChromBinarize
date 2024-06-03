args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
mark <- args[2] # m for 5mc and h for 5hmc

# percentage, 95 corresponds to <= 5 and >= 95
N_threshold <- as.numeric(args[3])
plot_1_name <- args[4]
plot_2_name <- args[5]

library(ggplot2)

methylation_data <- data.table::fread(bed_file)

names(methylation_data) <- c(
  "chr",
  "start",
  "end",
  "mark_name",
  "read_depth",
  "strand",
  "percent_methylation"
)

methylation_data <- methylation_data |>
  dplyr::filter(
    mark_name == mark,
    read_depth >= N_threshold 
  )


erroneous_unmethylated_p <- function(methylated_positions, percent) {
  methylated_positions <- methylated_positions |>
    dplyr::filter(percent_methylation >= percent) |>
    dplyr::mutate(
      "incorrectly_unmethylated" =
        read_depth * (100 - percent_methylation) / 100
    )
  incorrectly_unmethylated_total <-
    sum(methylated_positions$incorrectly_unmethylated)
  total_reads <- sum(methylated_positions$read_depth)
  erroneous_unmethylated_p <-
    incorrectly_unmethylated_total / total_reads

  return(erroneous_unmethylated_p)
}

erroneous_methylated_p <- function(unmethylated_positions, percent) {
  unmethylated_positions <- unmethylated_positions |>
    dplyr::filter(percent_methylation <= 100 - percent) |>
    dplyr::mutate(
      "incorrectly_methylated" =
        read_depth * percent_methylation / 100
    )
  incorrectly_methylated_total <-
    sum(unmethylated_positions$incorrectly_methylated)
  total_reads <- sum(unmethylated_positions$read_depth)
  erroneous_methylated_p <-
    incorrectly_methylated_total / total_reads

  return(erroneous_methylated_p)
}

percent_values <- 66:100
p1 <- lapply(percent_values, function(percent) {
  erroneous_unmethylated_p(methylation_data, percent)
})
p2 <- lapply(percent_values, function(percent) {
  erroneous_methylated_p(methylation_data, percent)
})

stats_table <- data.table::data.table(
  "percent" = percent_values,
  "p1" = unlist(p1),
  "p2" = unlist(p2)
)

p1_plot <- ggplot(stats_table, aes(x = percent, y = p1)) +
  geom_point(color = "black") +
  labs(
    x = "minimum percent methylation considered",
    y = "probability of erroneous unmethylation",
    title = "probability of erroneous unmethylation for binomial distribution"
  ) +
  theme_bw()

p2_plot <- ggplot(stats_table, aes(x = 100 - percent, y = p2)) +
  geom_point(color = "black") +
  labs(
    x = "maximum percent methylation considered",
    y = "probability of erroneous methylation",
    title = "probability of erroneous methylation for binomial distribution"
  ) +
  theme_bw() +
  scale_x_reverse()

# needed on some servers to actually create png files
options(bitmapType = "cairo")

ggsave(
  plot_1_name,
  plot = p1_plot
)

ggsave(
  plot_2_name,
  plot = p2_plot
)
