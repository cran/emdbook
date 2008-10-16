\name{traceplot.mcmc}
\alias{traceplot.mcmc}
\title{trace plot for BUGS run, in lattice graphics}
\description{
   Produces a trace plot (values of parameters for different
chains plotting against BUGS iterations) in the lattice package
}
\usage{
traceplot.mcmc(x, data = NULL, outer, aspect = "xy", default.scales = list(y = list(relation = "free")), start = 1, thin = 1, main = attr(x, "title"), xlab = "", ylab = "", plot.points = "rug", ..., subset = coda:::thinned.indices(x, start = start, thin = thin))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{mcmc} object (BUGS run)}
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{outer}{(stub)}
  \item{aspect}{aspect ratio setting (see \code{\link{xyplot}}}
  \item{default.scales}{scale setting (see \code{\link{xyplot}}}
  \item{start}{starting iteration number}
  \item{thin}{thinning factor}
  \item{main}{main title}
  \item{xlab}{x label}
  \item{ylab}{y label}
  \item{plot.points}{??}
  \item{\dots}{additional arguments to \code{\link{xyplot}}}
  \item{subset}{subset of data}
}
\value{
  produces a lattice plot 
}
\author{Ben Bolker}
\seealso{\code{\link{xyplot}}}
\keyword{hplot}
