% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{empty_tibble}
\alias{empty_tibble}
\title{Create empty tibble with column names and column types.}
\usage{
empty_tibble(nrow = 0, names, type = "number")
}
\arguments{
\item{nrow}{\code{Numeric} value defining number of desired rows. \code{Default = 0}.}

\item{names}{A vector of \code{characters} for desired column names.}

\item{type}{A \code{character} value for the desired type of all columns.
Takes \code{number}, 'character\verb{, or }factor\code{. }Default = number`.}
}
\description{
Create empty tibble with column names and column types.
}
\examples{
# empty tibble with zero rows and character columns
empty_tibble(nrow = 0,
             names = c("col1", "col2", "col3"),
             type = "character")

# empty tibble with three rows and numeric columns
empty_tibble(nrow = 3,
             names = c("col1", "col2", "col3"),
             type = "number")

# empty tibble with five rows and factor columns
empty_tibble(nrow = 5,
             names = c("col1", "col2", "col3"),
             type = "factor")

}
