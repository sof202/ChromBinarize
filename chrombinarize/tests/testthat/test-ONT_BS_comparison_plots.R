test_that("New columns added are as expected", {
  methylation_data <- readRDS(test_path(
    "test_data",
    "read_comparative_bedmethyl_files",
    "expected_comparitive_bedmethyl.Rds"
  ))
  actual_updated_columns <- add_absolute_change_columns(methylation_data)
  expected_updated_columns <- readRDS(test_path(
    "test_data",
    "ONT_BS_comparison",
    "expected_added_columns.Rds"
  ))

  expect_true(ncol(actual_updated_columns) == 10L)
  expect_identical(
    colnames(actual_updated_columns),
    colnames(expected_updated_columns)
  )
  expect_true(
    all(actual_updated_columns[["absolute_change_read_depth"]] >= 0)
  )
  expect_true(
    all(actual_updated_columns[["absolute_change_percent_methylation"]] >= 0)
  )
  expect_true(
    all(actual_updated_columns[["absolute_change_percent_methylation"]] <= 100)
  )
  expect_equal(actual_updated_columns, expected_updated_columns)
})
