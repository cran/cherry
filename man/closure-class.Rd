\name{Closure class}
\docType{class}
\alias{closure}
\alias{closure-class}
\alias{show,closure-method}
\alias{summary,closure-method}
\alias{defining}
\alias{defining,closure-method}
\alias{shortlist}
\alias{shortlist,closure-method}
\alias{hypotheses}
\alias{hypotheses,closure-method}
\alias{alpha,closure-method}
\alias{alpha<-,closure-method}
\alias{alpha<-}
\alias{alpha}

\title{Class "closure" for storing the result of the closed testing procedure}

\description{The class closure is the output of a call to \code{\link{closed}}. It stores the information needed to calculate confidence sets for the number of true and/or false hypotheses among a selected set.}

\section{Slots}{
  These slots are not meant to be directly accessed by the user.
  \describe{
    \item{\code{adjusted}:}{Object of class "numeric". Sores the adjusted p-values.}
    \item{\code{defining}:}{Object of class "integer". Stores the identifiers of the intersection hypotheses that are the defining rejections of the closed testing procedure. The identifiers should be read as binary booleans, i.e. 19 = 10011 (binary) is the intersection of the first, second and fifth elementary hypothesis.}
    \item{\code{hypotheses}:}{Object of class "character". Holds the names of the elementary hypotheses.}
    \item{\code{alpha}:}{Object of class "numeric". The type I error level chosen.}
    \item{\code{max.alpha}:}{Object of class "numeric". The largest value for which adjusted p-values have been calculated.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{(closure): Prints a brief description of the test results, including the upper bound of the number of true hypotheses and the corresponding lower bound to the number of false hypotheses among the full set.}
    \item{summary}{(closure): Prints the test results (as \code{show}) plus the defining rejections.}
    \item{defining}{(closure): Extracts the defining rejections as a list.}
    \item{shortlist}{(closure): Extracts the shortlist as a list.}
    \item{hypotheses}{(closure): Extracts the hypotheses slot.}
    }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{closed}}, \code{\link{pick}}.}

\keyword{methods}

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

  # perform the closed testing procedure
  cl <- closed(mytest, hypotheses)
  summary(cl)
  defining(cl)
  shortlist(cl)
  hypotheses(cl)

  # how many variables among a chosen set are associated with the response?
  pick(cl, c("ht", "lwt", "smoke", "ui"))
} 