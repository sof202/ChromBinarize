#' @export
read_methylation_data <- function(bed_file_location) {
  if (!file.exists(bed_file_location)) stop("ERROR: File does not exist.")

  methylation_data <- data.table::fread(
    bed_file_location,
    col.names = c(
      "chr",
      "start",
      "end",
      "name",
      "read_depth",
      "strand",
      "percent_methylation"
    ),
    colClasses = c(
      "character",
      "integer",
      "integer",
      "character",
      "integer",
      "character",
      "numeric"
    )
  )

  verify_column_class(
    methylation_data[["chr"]],
    is.character,
    "The first column (chromosome name) must be a string"
  )
  verify_column_class(
    methylation_data[["start"]],
    is.integer,
    "The second column (start) must be an integer"
  )
  verify_column_class(
    methylation_data[["end"]],
    is.integer,
    "The third column (end) must be an integer"
  )
  verify_column_class(
    methylation_data[["name"]],
    is.character,
    "The fourth column (name) must be a string"
  )
  verify_column_class(
    methylation_data[["read_depth"]],
    is.integer,
    "The fifth column (read_depth) must be an integer"
  )
  verify_column_class(
    methylation_data[["strand"]],
    is.character,
    "The sixth column (strand) must be a string"
  )
  verify_column_class(
    methylation_data[["percent_methylation"]],
    is.numeric,
    "The seventh column (percent_methylation) must be a numeric (float)"
  )

  return(methylation_data)
}
