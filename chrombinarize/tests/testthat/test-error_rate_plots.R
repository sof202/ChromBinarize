test_that("Expected value is obtained", {
  methylation_data <- readRDS(test_path(
    "test_data",
    "read_bedmethyl_files",
    "expected_bedmethyl_data.Rds"
  ))
  actual_error_rate_data <-
    chrombinarize:::create_error_rate_data(methylation_data, 60, 5)
  expected_error_rate_data <- readRDS(test_path(
    "test_data",
    "error_rate_plot_data",
    "expected_error_rate_data.Rds"
  ))
  expect_equal(actual_error_rate_data, expected_error_rate_data)
})
