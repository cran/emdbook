\name{gridsearch2d}
\alias{gridsearch2d}
\title{interactive grid search in 2D}
\description{
  For a specified objective function, computes a
  2D surface and allows the user to zoom in and
  visualize the surface
}
\usage{
gridsearch2d(fun, v1min, v2min, v1max, v2max, n1 = 20, n2 = 20, logz = FALSE, sys3d = c("both", "contour", "image"), ...)
}
\arguments{
  \item{fun}{A function taking two arguments}
  \item{v1min}{Starting minimum value for variable 1}
  \item{v2min}{Starting minimum value for variable 2}
  \item{v1max}{Starting maximum value for variable 1}
  \item{v2max}{Starting maximum value for variable 2}
  \item{n1}{Grid size for variable 1}
  \item{n2}{Grid size for variable 2}
  \item{logz}{Log scale? (see details)}
  \item{sys3d}{"both", "contour", or "image" for image type}
  \item{\dots}{Other arguments to graphing functions}
}
\details{
  The log scale used is log10(m-min(m)+mindm), where mindm
  is the difference between the minimum and the next-smallest
  value or 1e-10, whichever is larger.

  The function initially computes the value over the
  function and displays it as a contour, image, or overlaid
  contour ("both").  It then prompts the user to enter
  a new range by clicking on its corners.  If the user
  accepts the new range, it recomputes and redisplays (otherwise
  it just redisplays to erase the box).  When stopping,
  the user may opt to click on a final estimate.
}
\value{
  If requested, a list of the x and y coordinates of
  the final point
}
\author{Ben Bolker}
\keyword{misc}
