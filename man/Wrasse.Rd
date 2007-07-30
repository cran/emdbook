\name{Wrasse}
\alias{Wrasse}
\alias{WrasseSettlement}
\alias{WrasseRecruitment}
\alias{WrasseRecruitment_sum}
\docType{data}
\title{Reef fish (wrasse) data}
\description{
  Two data sets on \emph{Dascyllus trimaculatus}
(three-spotted wrasse), one on the distribution of settlement densities
  to empty anemones across time and space, the other on survival
(recruitment) of arriving settlers as a function of experimentally manipulated densities
  from Schmitt et al. (1999)
}
\usage{data(WrasseSettlement)
data(WrasseRecruitment)
data(WrasseRecruitment_sum)
}
\format{
  Three data frames:
  \describe{
    \item{\code{site}}{settlement site (location)}
    \item{\code{pulse}}{monthly settlement pulse}
    \item{\code{obs}}{observation within pulse}
    \item{\code{density}}{density of settlers per 0.1 m2 anemone}
    \item{\code{area}}{anemone area in cm2}
    \item{\code{init}}{initial settler density}
    \item{\code{surv}}{surviving density after 6 months}
    \item{\code{settler.den}}{target experimental density of settlers on
      experimental anemones}
    \item{\code{surv.den}}{mean surviving density after 6 months, by target density}
    \item{\code{SE}}{standard error of survivor density, by target density}
}
}
\source{
  Schmitt et al. (1999), "Quantifying the effects of multiple processes
  on local abundance", Ecology Letters 2:294-303.  (Original data kindly provided
  by Schmitt and Holbrook.)
}
\examples{
}
\keyword{datasets}
