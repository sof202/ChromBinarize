# Test data is a subset of hg19.txt that contains chromosomes 1-22, X and Y
chromosome_sizes_test_file <-
  test_path("test_data", "example_chromsizes.txt")

# This test bed data only includes chromosome 1 as it allows for a smaller
# rds file
expected_small_bed_data <-
  readRDS(test_path("test_data", "expected_bed_data.Rds"))

expected_chromosome_lengths <-
  readRDS(test_path("test_data", "expected_chromosome_lengths.Rds"))

testthat::test_that("Processing function returns correct vector", {
  chromosome_lengths <-
    chrombinarize:::process_chromosome_sizes(chromosome_sizes_test_file)

  expect_true(is.integer(chromosome_lengths))
  expect_equal(
    names(chromosome_lengths),
    paste0("chr", c(seq(1, 22), "X", "Y"))
  )
  expect_equal(chromosome_lengths, expected_chromosome_lengths)
})

testthat::test_that("Main function returns correct output for valid input", {
  bed_data <- create_blank_bed_data(
    chromosome_sizes_test_file,
    200,
    c(seq(1, 22), "X", "Y")
  )

  small_bed_data <- create_blank_bed_data(
    chromosome_sizes_test_file,
    200,
    "Y"
  )

  expect_true(is.list(bed_data))
  expect_equal(names(bed_data), c(seq(1, 22), "X", "Y"))
  expect_true(all(vapply(bed_data, inherits, logical(1), "data.table")))
  expect_equal(small_bed_data, expected_small_bed_data)
})
