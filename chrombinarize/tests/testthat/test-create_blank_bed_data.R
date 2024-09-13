# Test data is a subset of hg19.txt that contains chromosomes 1-22, X and Y
chromosome_sizes_test_file <-
  testthat::test_path("test_data", "example_chromsizes.txt")

# This test bed data only includes chromosome 1 as it allows for a smaller
# rds file
expected_small_bed_data <-
  readRDS(testthat::test_path("test_data", "expected_bed_data.Rds"))

testthat::test_that("Function returns correct output for valid input", {
  bed_data <- create_blank_bed_data(
    chromosome_sizes_test_file,
    200,
    c(seq(1, 22), "X", "Y")
  )

  small_bed_data <- create_blank_bed_data(
    chromosome_sizes_test_file,
    200,
    1
  )

  expect_true(is.list(bed_data))
  expect_equal(names(bed_data), c(seq(1, 22), "X", "Y"))
  expect_true(all(vapply(bed_data, inherits, logical(1), "data.table")))
  expect_equal(small_bed_data, expected_small_bed_data)
})
