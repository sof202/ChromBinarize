#' @export
reverse_binomial <- function(n, p, percent_methylation) {
  x <- ceiling(n * percent_methylation / 100)
  return(1 - pbinom(x - 1, n, p))
}
