\name{curve3d}
\alias{curve3d}
\title{Plot a 3D surface representing a 2D curve}
\description{
  Two-dimensional analogue of \code{curve}: generates a surface
  and plots it
}
\usage{
curve3d(expr, from = c(0, 0), to = c(1, 1), n = c(41, 41), add = FALSE,
zlab = NULL, log = NULL, sys3d = c("persp", "wireframe", "rgl",
"contour", "image", "none"), ...)
}
\arguments{
  \item{expr}{a mathematical expression using \code{x} and \code{y} as
    the independent variables}
  \item{from}{minimum values for \code{x} and \code{y}}
  \item{to}{maximum values for \code{x} and \code{y}}
  \item{n}{number of grid points in each direction}
  \item{add}{(logical) add to an existing plot?  (only possible for
    contour plots or rgl)}
  \item{zlab}{z label}
  \item{log}{(character): \code{"x"}, \code{"y"}, or \code{"xy"} for
    logarithmic axes}
  \item{sys3d}{3D plotting system to use: one of
    \code{"persp"}, \code{"wireframe"},
    \code{"rgl"}, \code{"contour"}, \code{"image"},
   or \code{"none"}}
  \item{\dots}{additional arguments to the plotting functions}
}
\value{
  invisibly, a list of
  \item{x}{x values}
  \item{y}{y values}
  \item{z}{z matrix}
}
\author{Ben Bolker}
\seealso{\code{\link{outer}}, \code{\link{curve}}}
\examples{

}
\keyword{hplot}
