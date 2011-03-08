\name{Select with shortcuts}

\alias{pickFisher}
\alias{pickSimes}
\alias{curveFisher}
\alias{curveSimes}


\title{Confidence limits for the number of true hypotheses, with shortcuts.}

\description{Calculates confidence limits for the number of true hypotheses on the basis of the closed testing procedure using specific local tests that allow shortcuts.}

\usage{
pickFisher (p, select = seq_along(p), alpha=0.05, silent=FALSE)

curveFisher (p, select = seq_along(p), alpha=0.05, plot = TRUE)

pickSimes (p, select = seq_along(p), alpha=0.05, hommel=FALSE, silent=FALSE)

curveSimes (p, select = seq_along(p), alpha=0.05, hommel=FALSE, plot = TRUE)
}

\arguments{
    \item{p}{The vector of p-values for all tested hypotheses.}
    \item{select}{The indexing vector of the p-values of the hypotheses to be selected. May be any type of appropriate indexing vector (integers, logical, or character).}
    \item{alpha}{The significance level of the test procedure.}
    \item{silent}{If \code{FALSE}, prints verbose result to the screen.}
    \item{plot}{If \code{TRUE}, plots the curve of correct rejections versus total rejections.}
    \item{hommel}{If \code{TRUE}, uses Hommel's (1983) variant of Simes' inequality. This variant is more conservative, but is valid under any form of dependence of the p-values.}
}

\details{The results of the \code{pickFisher} and \code{pickSimes} functions are identical to applying \code{\link{closed}} and \code{\link{pick}}, for specific choices of the local test, but are computationally more efficient. \code{pickFisher} uses local tests based on Fisher combinations. This local test is only valid if p-values of true hypotheses are independent. \code{pickSimes} uses a local test based on Simes' inequality. It is valid if p-values of true hypotheses are independent but also under some forms of positive correlations. The Hommel variant of the Simes local test is valid under any dependence structure of the p-values.}

\value{For \code{pickFisher} and \code{pickSimes}, the function returns the lower confidence limit for the number of false hypotheses (correct rejection) among the set \code{reject}. The upper confidence limit is always equal to the number of rejections made. \code{curveFisher} and \code{curveSimes} return the same confidence limit, but for selecting only the top 1,2,3,... hypotheses from the selected set.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\examples{
  # Fisher's method
  data(NAEP)
  pickFisher(NAEP, c("NH", "NC", "IA"))
  pickFisher(NAEP, 1:7)
  curveFisher(NAEP)

  # Simes method
  pickSimes(NAEP, c("NH", "NC", "IA"))
  pickSimes(NAEP, 1:7)
  curveSimes(NAEP)
}

\keyword{htest}
