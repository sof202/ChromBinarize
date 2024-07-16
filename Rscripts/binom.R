args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]

## =============== ##
##   OBTAINING P   ##
## =============== ##

# Sites we are reasonably sure are in fact unmethylated
unmethylated_positions <-
  data.table::fread(file.path(folder, "unmethylated_reads.bed"))

names(unmethylated_positions) <- c("reads", "percent_methylated")

unmethylated_positions <-
  dplyr::mutate(
    unmethylated_positions,
    "incorrectly_methylated" = reads * percent_methylated / 100
  )

incorrectly_methylated_total <-
  sum(unmethylated_positions$incorrectly_methylated)
total_reads <- sum(unmethylated_positions$reads)
erroneous_methylated_p <-
  incorrectly_methylated_total / total_reads

## ======= ##
##   MAIN  ##
## ======= ##

reverse_binomial <- function(n, p, percent_methylation) {
  x <- ceiling(n * percent_methylation / 100)
  return(1 - pbinom(x - 1, n, p))
}

methylation_data <- data.table::fread(file.path(folder, "filtered_reads.bed"))
names(methylation_data) <- c(
  "Chr",
  "start",
  "end",
  "name",
  "sample_size",
  "strand",
  "percent_methylation"
)

methylation_data <- dplyr::mutate(methylation_data,
  sample_size = as.numeric(sample_size),
  percent_methylation = as.numeric(percent_methylation),
  "likelihood of methylated reads being erroneous" =
    reverse_binomial(
      sample_size,
      erroneous_methylated_p,
      percent_methylation
    )
)

data.table::fwrite(methylation_data,
  file = file.path(folder, "processed_reads.bed"),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
