test_that("bedmethyl with correct format can be read", {
  example_bedmethyl_path <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "example_bedmethyl.bed"
  )
  expected_result <- readRDS(test_path(
    "test_data",
    "read_bedmethyl_files",
    "expected_bedmethyl_data.Rds"
  ))
  actual_result <- read_bedmethyl(example_bedmethyl_path)
  expect_equal(actual_result, expected_result)
})

test_that("bedmethyl with incorrect format returns an error", {
  # There must be exactly 7 columns
  wrong_number_of_columns <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "wrong_column_number.bed"
  )
  expect_error(
    read_bedmethyl(wrong_number_of_columns),
    info = "There must be exactly 7 columns"
  )

  wrong_column_type <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "wrong_column_types.bed"
  )
  expect_error(
    read_bedmethyl(wrong_column_type),
    info = "Columns must be of the specified types listed in the function"
  )

  invalid_mark_names <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "invalid_mark_names.bed"
  )
  expect_error(
    read_bedmethyl(invalid_mark_names),
    info = "Mark names must be 'h' or 'm' only"
  )

  negative_values <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "negative_values.bed"
  )
  expect_error(
    read_bedmethyl(negative_values),
    info = "No numeric column should have a negative value"
  )

  bad_region_size <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "bad_region_size"
  )
  expect_error(
    read_bedmethyl(bad_region_size),
    info = "No region should have negative length"
  )

  bad_percent_methylation <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "bad_percent_methylation"
  )
  expect_error(
    read_bedmethyl(bad_percent_methylation),
    info = "methylation shouldn't be greater than 100"
  )
})

test_that("Comparitive bedmethyl with correct format can be read", {
  example_bedmethyl_path <- test_path(
    "test_data",
    "read_comparison_files",
    "example_comparitive_bedmethyl.bed"
  )
  expected_result <- readRDS(test_path(
    "test_data",
    "read_comparison_files",
    "expected_comparitive_bedmethyl.Rds"
  ))
  actual_result <- read_comparison_bedmethyl(example_bedmethyl_path)
  expect_equal(actual_result, expected_result)
})

test_that("Comparitive bedmethyl with incorrect format returns an error", {
  wrong_number_of_columns <- test_path(
    "test_data",
    "read_comparison_files",
    "wrong_column_number.bed"
  )
  expect_error(
    read_bedmethyl(wrong_number_of_columns),
    info = "There must be exactly 8 columns"
  )

  wrong_column_type <- test_path(
    "test_data",
    "read_comparison_files",
    "wrong_column_types.bed"
  )
  expect_error(
    read_bedmethyl(wrong_column_type),
    info = "Columns must be of the specified types listed in the function"
  )

  invalid_mark_names <- test_path(
    "test_data",
    "read_comparison_files",
    "invalid_mark_names.bed"
  )
  expect_error(
    read_bedmethyl(invalid_mark_names),
    info = "Mark names must be 'h' or 'm' only"
  )

  negative_values <- test_path(
    "test_data",
    "read_comparison_files",
    "negative_values.bed"
  )
  expect_error(
    read_bedmethyl(negative_values),
    info = "No numeric column should have a negative value"
  )

  bad_region_size <- test_path(
    "test_data",
    "read_comparison_files",
    "bad_region_size"
  )
  expect_error(
    read_bedmethyl(bad_region_size),
    info = "No region should have negative length"
  )

  bad_percent_methylation <- test_path(
    "test_data",
    "read_comparison_files",
    "bad_percent_methylation"
  )
  expect_error(
    read_bedmethyl(bad_percent_methylation),
    info = "methylation shouldn't be greater than 100"
  )
})
