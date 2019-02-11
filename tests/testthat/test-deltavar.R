library(testthat)
library(emdbook)
fakegamma <- function(...) gamma(...)

## to test whether
test_msg <- function(msg,silent=FALSE) {
   grepl(msg,try(x ,silent=silent))
}
french_msg <- function() {
    test_msg("introuvable", silent=TRUE)
}
context("deltavar",

        test_that("numeric derivs fallback", {
            Sys.setlocale("LC_MESSAGES", "C")
            m <- c(scale=0.8,shape=12)
            S <- matrix(c(0.015,0.125,0.125,8.97), nrow=2)
            expect_warning(d1 <- deltavar(scale*fakegamma(1+1/shape),
                                  meanval=m,
                                  Sigma=S))
            Sys.setenv(
