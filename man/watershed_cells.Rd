% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Segmentation.R
\name{watershed_cells}
\alias{watershed_cells}
\title{Watershed function}
\usage{
watershed_cells(imclean, tol = 0.1, verbose = FALSE)
}
\arguments{
\item{imclean}{Image of class 'Image'}

\item{tol}{tolerance used for watershedding \link{default: 0.1}}

\item{verbose}{Print messages}
}
\description{
Take an image of class 'Image' as input and applies a watershed
algorithm of split merged cells.
}
