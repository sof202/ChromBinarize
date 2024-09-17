#' @title Read in a BEDMethyl File
#'
#' @description Converts a BEDMethyl file into a data table that can be used
#'  by many other functions in the pipeline.
#'
#' @param bed_file_location A file path (string) to the BEDMethyl file to be
#'  processed.
#'
#' @return A data.table with columns: "chr" (chromosome name, string), "start",
#'  (starting base pair position, int), "end" (ending base pair position, int),
#'  "mark_name" (m for 5mC and h for 5hmC), "read_depth" (read depth, int),
#'  "strand" (Either "+", "-" or "." for unknown strand) and
#'  "percent_methylation" (The percentage of such reads that were observed to
#'  be methylated, numeric)
#'
#' @details This function also checks that your BEDmethyl file is of the
#'  correct form. If your file has the incorrect number of columns or incorrect
#'  column structure, a helpful error message will be returned instead.
#'
#' @examples
#' read_bedmethyl("path/to/bedmethyl_file.bed")
#'
#'        chr start   end   name read_depth strand percent_methylation
#'     <char> <int> <int> <char>      <int> <char>               <num>
#'  1:   chr1   100   200      m         55      +                  20
#'  2:   chr1   300   400      h         20      +                  25
#'  3:   chr2   150   250      m         10      -                  20
#'  4:   chr2   500   600      h          1      +                   0
#'  5:   chr3   800   900      m         43      -                  49
#'  6:   chr3   950  1050      h         72      +                 100
#'
#' @export
read_bedmethyl <- function(bed_file_location) {
  if (!file.exists(bed_file_location)) stop("ERROR: File does not exist.")

  methylation_data <- data.table::fread(
    bed_file_location,
    col.names = c(
      "chr",
      "start",
      "end",
      "mark_name",
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

  is_mh <- function(ch) {
    return(ch == "m" || ch == "h")
  }

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
    methylation_data[["mark_name"]],
    is_mh,
    "The fourth column (mark_name) must only contain strings 'm' or 'h'"
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

#' @title Read in a comparison BEDMethyl File
#'
#' @description Converts a comparison BEDMethyl file into a data table that 
#'  can be used by the comparison scripts in the pipeline.
#'
#' @param bed_file_location A file path (string) to the comparison BEDMethyl
#' file to be processed.
#'
#' @return A data.table with columns: "chr" (chromosome name, string), "start",
#'  (starting base pair position, int), "end" (ending base pair position, int),
#'  "mark_name" (m for 5mC and h for 5hmC), "read_depth" (read depth, int),
#'  "strand" (Either "+", "-" or "." for unknown strand) and
#'  "percent_methylation" (The percentage of such reads that were observed to
#'  be methylated, numeric)
#'
#' @details This function also checks that your BEDmethyl file is of the
#'  correct form. If your file has the incorrect number of columns or incorrect
#'  column structure, a helpful error message will be returned instead.
#'
#' @examples
#' read_comparison_bedmethyl("path/to/bedmethyl_file.bed")
#'
#'        chr start   end mark_name ONT_read_depth ONT_percent_methylation
#'     <char> <int> <int>    <char>          <int>                   <num>
#'  1:   chr1   100   200         m             12                     8.0
#'  2:   chr1   300   400         h             22                    13.6
#'  3:   chr2   150   250         m             32                    50.0
#'  4:   chr2   500   600         h             25                     0.0
#'  5:   chr3   800   900         m             11                     9.0
#'  6:   chr3   950  1050         h             12                     8.0
#'     BS_read_depth BS_percent_methylation
#'             <int>                  <num>
#'  1:            14                   20.0
#'  2:             3                   66.7
#'  3:            40                  100.0
#'  4:            35                  100.0
#'  5:            32                   50.0
#'  6:            12                  100.0
#'
#' @export
read_comparison_bedmethyl <- function(bed_file_location) {
  if (!file.exists(bed_file_location)) stop("ERROR: File does not exist.")

  methylation_data <- data.table::fread(
    bed_file_location,
    col.names = c(
      "chr",
      "start",
      "end",
      "mark_name",
      "ONT_read_depth",
      "ONT_percent_methylation",
      "BS_read_depth",
      "BS_percent_methylation"
    ),
    colClasses = c(
      "character",
      "integer",
      "integer",
      "character",
      "integer",
      "numeric",
      "integer",
      "numeric"
    )
  )

  is_mh <- function(ch) {
    return(ch == "m" || ch == "h")
  }

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
    methylation_data[["mark_name"]],
    is_mh,
    "The fourth column (mark_name) must only contain strings 'm' or 'h'"
  )
  verify_column_class(
    methylation_data[["ONT_read_depth"]],
    is.integer,
    "The fifth column (ONT_read_depth) must be an integer"
  )
  verify_column_class(
    methylation_data[["ONT_percent_methylation"]],
    is.numeric,
    "The sixth column (ONT_percent_methylation) must be a numeric (float)"
  )
  verify_column_class(
    methylation_data[["BS_read_depth"]],
    is.integer,
    "The seventh column (BS_read_depth) must be an integer"
  )
  verify_column_class(
    methylation_data[["BS_percent_methylation"]],
    is.numeric,
    "The eighth column (BS_percent_methylation) must be a numeric (float)"
  )

  return(methylation_data)
}
