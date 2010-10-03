\name{Select with shortcuts}

\alias{pickFisher}
\alias{pickSimes}

\title{Confidence limits for the number of true hypotheses, with shortcuts.}

\description{Calculates confidence limits for the number of true hypotheses on the basis of the closed testing procedure using specific local tests that allow shortcuts.}

\usage{
pickFisher (all, select = seq_along(all), alpha=0.05, silent=FALSE)

pickSimes (all, select = seq_along(all), alpha=0.05, hommel = FALSE, silent=FALSE)
}

\arguments{
    \item{all}{The vector of p-values for all tested hypotheses.}
    \item{select}{The indexing vector of the p-values of the hypotheses to be selected. May be any type of appropriate indexing vector (integers, logical, or character).}
    \item{alpha}{The significance level of the test procedure.}
    \item{silent}{If \code{FALSE}, prints the result to the screen.}
    \item{hommel}{If \code{TRUE}, uses Hommel's (1983) variant of Simes' inequality. This variant is more conservative, but is valid under any form of dependence of the p-values.}
}

\details{The results of these functions are identical to applying \code{\link{closed}} and \code{\link{pick}}, for specific choices of the local test, but are computationally more efficient. \code{pickFisher} uses local tests based on Fisher combinations. This local test is only valid if p-values of true hypotheses are independent. \code{pickSimes} uses a local test based on Simes' inequality. It is valid if p-values of true hypotheses are independent but also under some forms of positive correlations. The Hommel variant of the Simes local test is valid under any dependence structure of the p-values.}

\value{The function returns the upper confidence limit for the number of true hypotheses among the set \code{reject}. The lower confidence limit is always equal to 0.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\examples{
  # Fisher's method
  pvalues <- c(p1 = 0.051, p2 = 0.064, p3 = 0.097, p4 = 0.108)
  pickFisher(pvalues, c("p1", "p3"))
  pickFisher(pvalues, -4)

  # Simes method
  pickSimes(pvalues)
}

\keyword{htest}
