mock_comparison_bedmethyl <- data.table::data.table(
  "chr" = rep("chr1", 8),
  "start" = c(0, 0, 200, 200, 400, 400, 600, 600),
  "end" = c(200, 200, 400, 400, 600, 600, 800, 800),
  "mark_name" = c("h", "m", "h", "m", "h", "m", "h", "m"),
  "ONT_read_depth" = c(6, 6, 99, 99, 12, 12, 32, 32),
  "ONT_percent_methylation" = c(5, 5, 10, 20, 0, 50, 30, 31),
  "BS_read_depth" = rep(0, 8),
  "BS_percent_methylation" = rep(0, 8)
)

expected_output <- data.table::data.table(
  "chr" = rep("chr1", 4),
  "start" = c(0, 200, 400, 600),
  "end" = c(200, 400, 600, 800),
  "mark_name" = c("m", "m", "m", "m"),
  "ONT_read_depth" = c(6, 99, 12, 32),
  "ONT_percent_methylation" = c(10, 30, 50, 61),
  "BS_read_depth" = rep(0, 4),
  "BS_percent_methylation" = rep(0, 4)
)

test_that("Normal functionality", {
  actual_output <- combine_5mc_5hmc(mock_comparison_bedmethyl)
  expect_true(all.equal(actual_output, expected_output))
  expect_true(
    all(actual_output[["mark_name"]] == "m"),
    info = "The mark_name column contains values other than 'm'"
  )
  expect_true(
    all(vapply(
      actual_output[["ONT_percent_methylation"]],
      inherits,
      logical(1),
      "numeric"
    )),
    info = "ONT_percent_methylation is expected to be a numeric column"
  )
})
