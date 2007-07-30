\name{lambertW}
\alias{lambertW}
\title{Lambert W function}
\description{
  Computes the Lambert W function, giving efficient solutions to the equation x*exp(x)==z
}
\usage{
lambertW(z, b = 0, maxiter = 10, eps = .Machine$double.eps, min.imag = 1e-09)
}
\arguments{
  \item{z}{(complex) vector of values for which to compute the function}
  \item{b}{(integer) vector of branches: b=0 specifies the principal
    branch, 0 and -1 are the ones that can take non-complex values}
  \item{maxiter}{maximum numbers of iterations for convergence}
  \item{eps}{convergence tolerance}
  \item{min.imag}{maximum magnitude of imaginary part to chop when
    returning solutions}
    }
\details{
Compute the Lambert W function of z.  This function satisfies
W(z)*exp(W(z)) = z, and can thus be used to express solutions
of transcendental equations involving exponentials or logarithms.
In ecology, the Lambert W can be used to solve the so-called
"Rogers equation" for predator functional response with depletion.
}
\value{
  Complex or real vector of solutions.
}
\references{Corless, Gonnet, Hare, Jeffrey, and Knuth (1996), "On the Lambert
W Function", Advances in Computational Mathematics 5(4):329-359}
\author{Nici Schraudolph <schraudo@inf.ethz.ch> (original
  version (c) 1998), Ben Bolker (R translation)
}
\seealso{\code{?Lambert} in the \code{gsl} package by Robin Hankin}
\note{
This implementation should return values within 2.5*eps of its
counterpart in Maple V, release 3 or later.  Please report any
discrepancies to the author or translator.
}
\examples{
curve(lambertW(x),from=0,to=10)
pvec = seq(-1,1,length=40)
m = outer(pvec,pvec,function(x,y)Re(lambertW(x+y*1i)))
persp(pvec,pvec,m,
      theta=290,shade=0.5,zlab="lambertW")
num1 = uniroot(function(x) {x*exp(x)-1},lower=0,upper=1,tol=1e-9)
abs(lambertW(1)-num1$root)<1e-9
###
## Rogers random predator equation:
rogers.pred <- function(N0,a,h,T) {
   N0 - lambertW(a*h*N0*exp(-a*(T-h*N0)))/(a*h)
}
holling2.pred <- function(N0,a,h) {
  a*N0/(1+a*h*N0)
}
curve(rogers.pred(x,a=1,h=0.2,T=1),from=0,to=60,
  ylab="Number eaten/unit time",xlab="Initial number",ylim=c(0,5),
  main="Predation: a=1, h=0.2")
curve(rogers.pred(x,a=1,h=0.2,T=5)/5,add=TRUE,lty=2,from=0)
curve(rogers.pred(x,a=1,h=0.2,T=0.2)*5,add=TRUE,lty=3,from=0)
curve(rogers.pred(x,a=1,h=0.2,T=10)/10,add=TRUE,lty=4,from=0)
curve(holling2.pred(x,a=1,h=0.2),add=TRUE,lty=1,lwd=2,from=0)
abline(h=5)
legend(30,2,
   c(paste("Rogers, T=",c(0.2,1,5,10),sep=""),
    "Holling type II"),lwd=c(rep(1,4),2),lty=c(3,1,2,4,1))
## final size of an epidemic
finalsize = function(R0) {
   1+1/R0*lambertW(-R0*exp(-R0))
}
curve(finalsize,from=1,to=10,xlab=expression(R[0]),ylab="Final size")
}
\keyword{math}
