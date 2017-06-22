\name{colNorm}
\alias{colNorm}
\title{Column norm of a matrix}
\description{
  \code{colNorm} provides a column norm of a matrix to provide an
  example of using \pkg{RcppGSL} with \code{GNU GSL} library.
}
\usage{
colNorm(M)
}
\arguments{
  \item{M}{a numeric matrix.}
}
\details{
  The column norm of \code{M} is returned.

  This example reprises an example in section 8.4.13 of the GSL
  reference manual. 
}

\value{
  \code{colNorm} returns a vector each element of which corresponds to
  the vector norm of the corresponding column of \code{M}.
}
\references{GNU GSL project: \url{http://www.gnu.org/software/gsl}}
\author{
  The GNU GSL library is being written by team of authors with the
  overall development, design and implementation lead by Brian Gough and
  Gerard Jungman.  RcppGSL is written by Romain Francois and Dirk Eddelbuettel.
}
\examples{
  ## see Section 8.4.13 of the GSL manual
  ## create M as a sum of two outer products
  M <- outer(sin(0:9), rep(1,10), "*") + outer(rep(1, 10), cos(0:9), "*")
  print(colNorm(M))
  ## same result using just R
  print(apply(M, 2, function(x) sqrt(sum(x^2))))
}

