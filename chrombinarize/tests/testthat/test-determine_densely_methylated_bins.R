bin_counts_test_file <-
  testthat::test_path("test_data", "example_bin_counts.bed")

expected_dense_bins <-
  readRDS(testthat::test_path("test_data", "expected_dense_bins.Rds"))

test_that("Basic functionality test", {
  result <- determine_dense_bins(bin_counts_test_file)
  expect_equal(ncol(result), 6)
  expect_equal(nrow(result), 51)
  expect_equal(
    result[["densely_methylated"]],
    expected_dense_bins[["densely_methylated"]]
  )
})

test_that("Throws error if file cannot be read", {
  expect_error(read_bin_counts_file("non-existent-file"))
})

test_that("Throws error if file has the incorrect number of columns", {
  tempfile_path <- withr::local_tempfile(
    pattern = "not_enough_columns.txt",
    lines = c("string,string", "string,string")
  )
  expect_error(read_bin_counts_file(tempfile_path))

  tempfile_path <- withr::local_tempfile(
    pattern = "too_many_columns.txt",
    lines = c(
      "string,string,string,string,string",
      "string,string,string,string,string"
    )
  )
  expect_error(read_bin_counts_file(tempfile_path))
})


test_that("Throws error if file has incorrect column classes", {
  tempfile_path <- withr::local_tempfile(
    pattern = "end_column_wrong.txt",
    lines = c("chr1,0,string,2")
  )
  expect_error(read_bin_counts_file(tempfile_path))

  tempfile_path <- withr::local_tempfile(
    pattern = "start_column_wrong.txt",
    lines = c("chr1,string,0,2")
  )
  expect_error(read_bin_counts_file(tempfile_path))

  tempfile_path <- withr::local_tempfile(
    pattern = "count_column_wrong.txt",
    lines = c("chr1,0,0,string")
  )
  expect_error(read_bin_counts_file(tempfile_path))
})
