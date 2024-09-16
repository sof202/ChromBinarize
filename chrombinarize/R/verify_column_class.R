#' @title Check if a Column has the Correct Class
#'
#' @description To ensure that files are read in correctly, columns can be
#'  checked for if they are of the desired class.
#'
#' @param column A vector of any class
#' @param type_checker A base R type checker (i.e. `is.integer()`)
#' @param stop_message A string to be outputted in case of failure
#'
#' @return A custom stop error message if the column doesn't adhere to the
#'  desired class.
#'
#' @examples
#' # Returns nothing
#' verify_column_class(c(1, 2, 3, 4), is.integer, "column is not an integer")
#'
#' # Returns error message
#' verify_column_class(
#'   c("these", "are", "characters"),
#'   is.integer,
#'   "column is not an integer"
#' )
#' "column is not an integer" (and terminates call stack)
verify_column_class <- function(column, type_checker, stop_message) {
  if (!all(vapply(column, type_checker, logical(1)))) {
    stop(stop_message)
  }
}
