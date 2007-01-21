\name{lump.mcmc.list}
\alias{lump.mcmc.list}
\title{Combine multiple MCM chains}
\description{
  Converts an \code{mcmc.list} object to an \code{mcmc}
  object by combining chains for all variables
}
\usage{
lump.mcmc.list(x)
}
\arguments{
  \item{x}{an \code{mcmc.list} object}
}
\value{
  an \code{mcmc} object 
}
\author{Ben Bolker}
\seealso{\code{coda} package}
\keyword{misc}
