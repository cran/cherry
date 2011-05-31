\name{Pick}

\alias{pick}

\title{Confidence limits for the number of true hypotheses.}

\description{Calculates confidence limits for the number of true hypotheses on the basis of the closed testing procedure.}

\usage{
pick (closure, reject, silent=FALSE)
}

\arguments{
    \item{closure}{An object of class \code{\link{closure}}, typically created through a call to \code{\link{closed}} .}
    \item{reject}{A character vector identifying the hypotheses to be rejected. Must be a subvector of \code{hypotheses(closure)}.}
    \item{silent}{If \code{FALSE}, prints the result to the screen.}
}

\details{The function \code{pick} calculates a confidence interval for the number of true hypotheses among a selected set of hypotheses.}

\value{The function returns the upper confidence limit for the number of true hypotheses among the set \code{reject}. The lower confidence limit is always equal to 0.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\examples{
  # Example: the birthwt data set from the MASS library
  # We want to find variables associated with low birth weight
  library(MASS)
  fullfit <- glm(low~age+lwt+race+smoke+ptl+ht+ui+ftv, family=binomial, data=birthwt)
  hypotheses <- c("age", "lwt", "race", "smoke", "ptl", "ht", "ui", "ftv")

  # Define the local test to be used in the closed testing procedure
  mytest <- function(hyps) {
    others <- setdiff(hypotheses, hyps)
    form <- formula(paste(c("low~",  paste(c("1", others), collapse="+"))))
    anov <- anova(glm(form, data=birthwt, family=binomial), fullfit, test="Chisq")
    res <- anov$"Pr("[2]                    # for R >= 2.14.0
    if (is.null(res)) res <- anov$"P("[2]   # earlier versions
    res
  }

  # perform the closed testing
  cl <- closed(mytest, hypotheses)
  summary(cl)

  # how many variables among a chosen set are associated with the response?
  pick(cl, c("ht", "lwt", "smoke", "ui"))
}

\keyword{htest}
