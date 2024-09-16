verify_column_class <- function(column, type_checker, stop_message) {
  if (!all(vapply(column, type_checker, logical(1)))) {
    stop(stop_message)
  }
}
