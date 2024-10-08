% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verify_column_class.R
\name{verify_column_class}
\alias{verify_column_class}
\title{Check Column has Correct Class}
\usage{
verify_column_class(column, type_checker, stop_message)
}
\arguments{
\item{column}{The vector (extracted column) to inspect}

\item{type_checker}{A type checker (i.e. \code{is.integer()})}

\item{stop_message}{A string to be outputted in case of failure}
}
\value{
A custom stop error message if the column doesn't adhere to the
desired class.
}
\description{
To ensure that files are read in correctly, columns can be
checked for if they are of the desired class.
}
\examples{
# Returns nothing
verify_column_class(c(1, 2, 3, 4), is.integer, "column is not an integer")

# Returns error message
verify_column_class(
  c("these", "are", "characters"),
  is.integer,
  "column is not an integer"
)
Error: column is not an integer (and terminates call stack)
}
