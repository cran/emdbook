\name{scinot}
\alias{scinot}
\title{Scientific notation as LaTeX/expression()}
\description{
  Takes a number and returns a version formatted in LaTeX
  (suitable for use with \code{Sexpr()} in an Sweave document)
  or in \code{expression()} (suitable for plotting)
}
\usage{
scinot(x, format = c("latex", "expression"), delim="$",
pref="", ...)
}
\arguments{
  \item{x}{a numeric vector (of length 1)}
  \item{format}{produce LaTeX or expression() format?}
  \item{delim}{delimiter to add at beginning and end (latex only)}
  \item{pref}{text to put before expression (expression only)}
  \item{\dots}{additional arguments to \code{\link{formatC}}}
}
\value{
  a character vector (if \code{latex}) or expression (if \code{expression})
}
\author{Ben Bolker}
\seealso{\code{\link{formatC}}, \code{\link{expression}},
  \code{\link{plotmath}},
  \code{latexSN} in \code{Hmisc} package}
\examples{
scinot(1e-5)
scinot(1e-5,digits=0)
scinot(1e-5,"expression")
scinot(1e-5,"expression",pref="p=")
}
\keyword{misc}
