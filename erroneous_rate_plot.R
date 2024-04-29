args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
mark <- args[2] # m for 5mc and h for 5hmc
max_read_depth <- as.numeric(args[3]) # you don't need to consider N = 2000

# percentage, 95 corresponds to <= 5 and >= 95
methylation_threshold <- as.numeric(args[4])

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

methylated_positions <- methylation_data |>
  dplyr::filter(
    mark_name == mark,
    percent_methylation >= methylation_threshold
  )

unmethylated_positions <- methylation_data |>
  dplyr::filter(
    mark_name == mark,
    percent_methylation <= 100 - methylation_threshold
  )

erroneous_unmethylated_p <- function(methylated_positions, n) {
  methylated_positions <- methylated_positions |>
    dplyr::filter(read_depth >= n) |>
    dplyr::mutate(
      "incorrectly_unmethylated" =
        read_depth * (100 - percent_methylation) / 100
    )
  incorrectly_unmethylated_total <-
    sum(methylated_positions$incorrectly_unmethylated)
  total_reads <- sum(methylated_positions$reads)
  erroneous_unmethylated_p <-
    incorrectly_unmethylated_total / total_reads

  return(erroneous_unmethylated_p)
}

erroneous_methylated_p <- function(unmethylated_positions, n) {
  unmethylated_positions <- unmethylated_positions |>
    dplyr::filter(read_depth >= n) |>
    dplyr::mutate(
      "incorrectly_methylated" =
        read_depth * percent_methylation / 100
    )
  incorrectly_methylated_total <-
    sum(unmethylated_positions$incorrectly_methylated)
  total_reads <- sum(unmethylated_positions$reads)
  erroneous_methylated_p <-
    incorrectly_methylated_total / total_reads

  return(erroneous_methylated_p)
}


stats_table <- data.table::data.table("n" = 1:max_read_depth)

stats_table <- stats_table |>
  dplyr::mutate(
    "p1" = erroneous_unmethylated_p(unmethylated_positions, n),
    "p2" = erroneous_methylated_p(methylated_positions, n)
  )


p1_plot <- ggplot(stats_table, aes(x = n, y = p1)) +
  geom_smooth(color = "red", span = 1.2) +
  geom_point(color = "black") +
  labs(
    x = "min read depth considered",
    y = "probability of erroneous unmethylation"
  ) +
  theme_bw()

p2_plot <- ggplot(stats_table, aes(x = n, y = p2)) +
  geom_smooth(color = "red", span = 1.2) +
  geom_point(color = "black") +
  labs(
    x = "min read depth considered",
    y = "probability of erroneous methylation"
  ) +
  theme_bw()

options(bitmapType = "cairo")

ggsave(
  "p1.png",
  plot = p1_plot
)

ggsave(
  "p2.png",
  plot = p2_plot
)
