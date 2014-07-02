\name{Select with shortcuts}

\alias{pickFisher}
\alias{pickSimes}
\alias{pickMeinshausen}
\alias{curveFisher}
\alias{curveSimes}
\alias{curveMeinshausen}


\title{Confidence limits for the number of true hypotheses, with shortcuts.}

\description{Calculates confidence limits for the number of true hypotheses on the basis of the closed testing procedure using specific local tests that allow shortcuts.}

\usage{
pickFisher (p, select = seq_along(p), alpha=0.05, silent=FALSE)

curveFisher (p, select = seq_along(p), order, alpha=0.05, plot = TRUE)

pickSimes (p, select = seq_along(p), alpha=0.05, hommel=FALSE, silent=FALSE)

curveSimes (p, select = seq_along(p), order, alpha=0.05, hommel=FALSE, plot = TRUE)

pickMeinshausen (p, PM, select = seq_along(p), alpha=0.05, silent=FALSE)

curveMeinshausen (p, PM, select = seq_along(p), order, alpha=0.05, plot = TRUE) 
}

\arguments{
    \item{p}{The vector of p-values for all tested hypotheses.}
    \item{PM}{A matrix of permutation p-values. Rows are hypotheses; columns are permutations.}
    \item{select}{The indexing vector of the p-values of the hypotheses to be selected. May be any type of appropriate indexing vector (integers, logical, or character).}
    \item{order}{The indexing vector specifying the order in which p-values of the hypotheses are to be rejected. May be integer or character.}
    \item{alpha}{The significance level of the test procedure.}
    \item{silent}{If \code{FALSE}, prints verbose result to the screen.}
    \item{plot}{If \code{TRUE}, plots the curve of correct rejections versus total rejections.}
    \item{hommel}{If \code{TRUE}, uses Hommel's (1983) variant of Simes' inequality. This variant is more conservative, but is valid under any form of dependence of the p-values.}
}

\details{The results of the \code{pickFisher}, \code{pickSimes} and \code{pickMeinshausen} functions are identical to applying \code{\link{closed}} and \code{\link{pick}}, for specific choices of the local test, but are computationally more efficient. \code{pickFisher} uses local tests based on Fisher combinations. This local test is only valid if p-values of true hypotheses are independent. \code{pickSimes} uses a local test based on Simes' inequality. It is valid if p-values of true hypotheses are independent but also under some forms of positive correlations. The Hommel variant of the Simes local test is valid under any dependence structure of the p-values. \code{pickMeinshausen} is a permutation-based variant of \code{pickSimes}. See the reference below.

In the curve functions, the user may specify either \code{select} or \code{order}. Specifying \code{order} fixes the precise order in which hypotheses are selected, whereas specifying \code{select} only specifies which hypotheses are candidates for selection, leaving the order to be chosen by the function to maximize the number of correct rejections.}

\value{For \code{pickFisher}, \code{pickSimes} and \code{pickMeinshausen}, the function returns the lower confidence limit for the number of false hypotheses (correct rejection) among the set \code{reject}. The upper confidence limit is always equal to the number of rejections made. \code{curveFisher} and \code{curveSimes} return the same confidence limit, but for selecting only the first 1,2,3,... hypotheses from the selected set.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Aldo Solari}

\references{
Goeman and Solari (2011) Statistical Science 26 (4) 584-597.

Meinshausen (2006) Scandinavian Journal of Statistics 33 (2), 227-237.
}


\examples{
  # Fisher's method
  data(NAEP)
  pickFisher(NAEP, c("NH", "NC", "IA"))
  pickFisher(NAEP, 1:7)
  curveFisher(NAEP)
  curveFisher(NAEP, order=7:1)


  # Simes method
  pickSimes(NAEP, c("NH", "NC", "IA"))
  pickSimes(NAEP, 1:7)
  curveSimes(NAEP)
  curveSimes(NAEP, select=1:7)
  
  # Meinshausen's method
  library(multtest)
  data(golub)
  smallglb <- golub[1:500,]
  TM<-sapply(1:nrow(smallglb), function(i) {
    mt.sample.teststat(smallglb[i,], golub.cl,test="t.equalvar",B=500)
  })
  PM<-2*(1-pt(abs(TM),df=length(golub.cl)-2)) # permutation matrix
  
  # p-values
  p<-apply(smallglb,1, function(z) t.test(z[golub.cl==0],z[golub.cl==1],var.equal=TRUE)$p.value) 
  
  pickMeinshausen(p, PM, select=1:100)
  pickMeinshausen(p, PM, select=sort.list(p)[1:100])
  curveMeinshausen(p,PM, select=1:200)
  curveMeinshausen(p,PM, order=1:200)
}

\keyword{htest}
