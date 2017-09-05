## colNorm.R: R wrapper to Rcpp/GSL colNorm implementation
##
## Copyright (C)  2010 Dirk Eddelbuettel and Romain Francois
##
## This file is part of RcppGSL.
##
## RcppGSL is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppGSL is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppGSL.  If not, see <http://www.gnu.org/licenses/>.

## old call now shadowed by auto-generated colNorm() in RcppExports.R
colNorm_old <- function(M) {
    stopifnot(is.matrix(M))
    res <- .Call("colNorm_old", M, package="RcppGSLExample")
}

