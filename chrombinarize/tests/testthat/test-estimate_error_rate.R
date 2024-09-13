example_input_file <-
  testthat::test_path("test_data", "example_unmethylated_reads.bed")
expected_error_rate <-
  readRDS(testthat::test_path("test_data", "expected_error_rate.Rds"))

test_that("error rate is being calculated correctly", {
  calculated_error_rate <- estimate_error_rate(example_input_file)

  expect_true(is.numeric(calculated_error_rate))
  expect_equal(calculated_error_rate, expected_error_rate)
})


test_that("error rate is zero for all unmethylated reads", {
  tempfile_path <- withr::local_tempfile(
    pattern = "allZeros.bed",
    lines = c("10\t0", "20\t0")
  )
  calculated_error_rate <- estimate_error_rate(tempfile_path)
  expect_equal(calculated_error_rate, 0)
})
test_that("error rate is one for all methylated reads", {
  tempfile_path <- withr::local_tempfile(
    pattern = "allZeros.bed",
    lines = c("10\t100", "20\t100")
  )
  calculated_error_rate <- estimate_error_rate(tempfile_path)
  expect_equal(calculated_error_rate, 1)
})

test_that("function throws error for incorrect file format", {
  tempfile_path <- withr::local_tempfile(
    pattern = "allZeros.bed",
    lines = c("string,string")
  )
  expect_error(estimate_error_rate(tempfile_path))
})
