test_that("error rate is being calculated correctly", {
  normal_bedmethyl <- readRDS(test_path(
    "test_data",
    "estimate_error_rate",
    "normal_bedmethyl.Rds"
  ))
  actual_error_rate <- estimate_error_rate(normal_bedmethyl, 30, 5)
  expected_error_rate <- readRDS(test_path(
    "test_data",
    "estimate_error_rate",
    "normal_bedmethyl_error_rate.Rds"
  ))

  expect_true(is.numeric(actual_error_rate))
  expect_equal(actual_error_rate, expected_error_rate)
})

test_that("error rate is zero if all reads are unmethylated", {
  all_zeros <- readRDS(test_path(
    "test_data",
    "estimate_error_rate",
    "all_zeros.Rds"
  ))
  actual_error_rate <- estimate_error_rate(all_zeros, 30, 5)
  expect_equal(actual_error_rate, 0)
})

test_that("error rate is NaN if no reads hit threshold", {
  all_one_hundreds <- readRDS(test_path(
    "test_data",
    "estimate_error_rate",
    "all_one_hundreds.Rds"
  ))
  calculated_error_rate <- estimate_error_rate(all_one_hundreds, 30, 5)
  expect_equal(calculated_error_rate, NaN)

  all_low_read_depth <- readRDS(test_path(
    "test_data",
    "estimate_error_rate",
    "all_low_read_depth.Rds"
  ))
  calculated_error_rate <- estimate_error_rate(all_low_read_depth, 30, 5)
  expect_equal(calculated_error_rate, NaN)
})
