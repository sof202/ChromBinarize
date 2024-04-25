args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]

## =============== ##
##   OBTAINING P   ##
## =============== ##

# Sites we are reasonably sure are in fact methylated
methylated_positions <- data.table::fread(paste0(folder, "/methylated.csv"))

# Sites we are reasonably sure are in fact unmethylated
unmethylated_positions <- data.table::fread(paste0(folder, "/unmethylated.csv"))

names(methylated_positions) <- c("reads", "percent_methylated")
names(unmethylated_positions) <- c("reads", "percent_methylated")

methylated_positions <- dplyr::mutate(methylated_positions,
  "incorrectly_unmethylated" = reads * (100 - percent_methylated) / 100
)
unmethylated_positions <-
  dplyr::mutate(unmethylated_positions,
    "incorrectly_methylated" = reads * percent_methylated / 100
  )

incorrectly_unmethylated_total <-
  sum(methylated_positions$incorrectly_unmethylated)
total_reads <- sum(methylated_positions$reads)
erroneous_unmethylated_p <-
  incorrectly_unmethylated_total / total_reads

incorrectly_methylated_total <-
  sum(unmethylated_positions$incorrectly_methylated)
total_reads <- sum(unmethylated_positions$reads)
erroneous_methylated_p <-
  incorrectly_methylated_total / total_reads

## ======= ##
##   MAIN  ##
## ======= ##

methylation_data <- data.table::fread(paste0(folder, "/filtered_reads.bed"))
names(methylation_data) <-
  c("Chr", "start", "end", "name", "size", "strand", "percent_methylation")

reverse_binomial <- function(n, p, percent_methylation) {
  x <- ceiling(n * percent_methylation / 100)
  return(1 - pbinom(x - 1, n, p))
}

methylation_data <- dplyr::mutate(methylation_data,
  size = as.numeric(size),
  percent_methylation = as.numeric(percent_methylation),
  "likelihood of unmethylated reads being erroneous" =
    reverse_binomial(
      size,
      erroneous_unmethylated_p,
      (100 - percent_methylation)
    ),
  "likelihood of methylated reads being erroneous" =
    reverse_binomial(
      size,
      erroneous_methylated_p,
      percent_methylation
    )
)

data.table::fwrite(methylation_data,
  file = paste0(folder, "/processed_reads.bed"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
)
