\name{Closed Testing}

\alias{closed}

\title{Closed Testing for Exploratory Research}

\description{Performs the closed testing procedure for user-specified local test.}

\usage{
closed (test, hypotheses, alpha = 0.05)
}

\arguments{
    \item{test}{A function that performs the local test. The function should accept a subvector of the \code{hypotheses} argument as input, and return a p-value.}
    \item{hypotheses}{Identifiers of the collection of elementary hypotheses.}
    \item{alpha}{The significance level of the test procedure.}
}

\details{The function \code{closed} performs the closed testing procedure on the collection \code{hypotheses}, testing all their intersection hypotheses and controlling the familywise error rate.}

\note{The number of intersection hypotheses is exponential in the number of elementary hypotheses. The number of elementary hypotheses is therefore limited to \code{log2(.Machine$integer.max+1)} (typically 31) for computational reasons.}

\value{The function \code{closed} returns an object of class \code{\link{closure}}.}

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
