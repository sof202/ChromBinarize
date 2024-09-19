test_that("Normal functionality is obtained for data creation", {
  methylation_data <- readRDS(test_path(
    "test_data",
    "read_bedmethyl_files",
    "expected_bedmethyl_data.Rds"
  ))
  actual_cpg_robustness_data <-
    chrombinarize:::create_cpg_robustness_data(methylation_data, 0, 1000)
  expected_cpg_robustness_data <- readRDS(test_path(
    "test_data",
    "cpg_robustness",
    "expected_cpg_robustness_data.Rds"
  ))
  expect_equal(actual_cpg_robustness_data, expected_cpg_robustness_data)

  # NaN values should be removed to avoid plotting errors
  expect_false(
    any(is.na(actual_cpg_robustness_data[["surrounding_methylation"]]))
  )
})

test_that("Edge cases are handled as expected", {
  methylation_data <- readRDS(test_path(
    "test_data",
    "read_bedmethyl_files",
    "expected_bedmethyl_data.Rds"
  ))

  # Filtering steps will always leave an empty data frame
  expect_no_error(
    chrombinarize:::calculate_nearby_methylation(
      methylation_data,
      "chr1",
      1000,
      0,
      0
    )
  )

  # Chromosome doesn't exist
  expect_no_error(
    chrombinarize:::calculate_nearby_methylation(
      methylation_data,
      "chrZ",
      1000,
      0,
      30
    )
  )

  # Base pair position doesn't exist
  expect_no_error(
    chrombinarize:::calculate_nearby_methylation(
      methylation_data,
      "chr1",
      -1,
      0,
      30
    )
  )
})
