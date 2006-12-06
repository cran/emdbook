\name{dbetabinom}
\alias{dbetabinom}
\alias{rbetabinom}
\title{Beta-binomial distribution}
\description{
   Density function and random variate generator 
for the beta-binomial function, parameterized in
terms of probability and overdispersion
}
\usage{
dbetabinom(x, prob, size, theta, log = FALSE)
rbetabinom(n, prob, size, theta)
}
\arguments{
  \item{x}{a numeric vector of values}
  \item{prob}{numeric vector: mean probability of underlying beta distribution}
  \item{size}{integer: number of samples}
  \item{theta}{overdispersion parameter}
  \item{log}{(logical) return log probability density?} 
  \item{n}{integer number of random variates to return}
}
\details{
   The beta-binomial distribution is the result of
compounding a beta distribution of probabilities with
a binomial sampling process.  The density function is ...
  v <- lchoose(size,x)-lbeta(theta*(1-prob),theta*prob)+lbeta(size-x+theta*(1-prob),x+theta*prob)
(alternative parameterization)
}
\value{
   A vector of probability densities or random deviates.
}
\references{ ~put references to the literature/web site here ~ }
\author{Ben Bolker}
\seealso{\code{\link{dbeta}}, \code{\link{dbinom}}}
\examples{
  set.seed(100)
  n <- 9
  z <- rbetabinom(1000, 0.5, size=n, theta=4)
  plot(table(z)/length(z),ylim=c(0,0.34),col="gray",lwd=4)
  points(0:n,dbinom(0:n,size=n,prob=0.5),col=2,pch=16,type="b")
  points(0:n,dbetabinom(0:n,size=n,theta=4,
           prob=0.5),col=3,pch=17,type="b")
}
\keyword{misc}
