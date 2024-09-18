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
  wrong_number_of_columns <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "wrong_column_number.bed"
  )
  expect_error(read_bedmethyl(wrong_number_of_columns))

  wrong_column_type <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "wrong_column_types.bed"
  )
  expect_error(read_bedmethyl(wrong_column_type))

  invalid_mark_names <- test_path(
    "test_data",
    "read_bedmethyl_files",
    "invalid_mark_names.bed"
  )
  expect_error(read_bedmethyl(invalid_mark_names))
})

test_that("Comparitive bedmethyl with correct format can be read", {
  example_bedmethyl_path <- test_path(
    "test_data",
    "read_comparative_bedmethyl_files",
    "example_comparitive_bedmethyl.bed"
  )
  expected_result <- readRDS(test_path(
    "test_data",
    "read_comparative_bedmethyl_files",
    "expected_comparitive_bedmethyl.Rds"
  ))
  actual_result <- read_comparison_bedmethyl(example_bedmethyl_path)
  expect_equal(actual_result, expected_result)
})

test_that("Comparitive bedmethyl with incorrect format returns an error", {
  wrong_number_of_columns <- test_path(
    "test_data",
    "read_comparative_bedmethyl_files",
    "wrong_column_number.bed"
  )
  expect_error(read_bedmethyl(wrong_number_of_columns))

  wrong_column_type <- test_path(
    "test_data",
    "read_comparative_bedmethyl_files",
    "wrong_column_types.bed"
  )
  expect_error(read_bedmethyl(wrong_column_type))

  invalid_mark_names <- test_path(
    "test_data",
    "read_comparative_bedmethyl_files",
    "invalid_mark_names.bed"
  )
  expect_error(read_bedmethyl(invalid_mark_names))
})
