\name{regionpick}

\alias{regionpick}

\title{Confidence limits for the number of false hypotheses in a given region.}

\description{Calculates confidence limits for the number of false hypotheses on the basis of the region procedure within one or more regions.}

\usage{ regionpick (region, intervals, alpha, silent = FALSE, ignore_weights = TRUE) }

\arguments{
  \item{region}{An object of class \code{\link{region}}, typically created through a call to \code{\link{regionmethod}}.}
  \item{intervals}{A list containing one or more regions, specified by a left and rightbound.}
  \item{alpha}{For region objects with adjusted p-values, specifies the value of alpha for which confidence limits are to be calculated (optional).}
  \item{silent}{If \code{FALSE}, prints the result to the screen.}
  \item{ignore_weights}{If set to TRUE, a confidence interval for the number of false elementary hypotheses will be computed. If set to FALSE, a confidence interval for the combined weight of false elemenatary hypotheses will be computed.}
}

\details{The function \code{regionpick} calculates a confidence interval for the number (or weight) of false hypotheses among a selected set of elementary hypotheses.}

\value{The function returns the lower confidence limit for the number of false hypotheses (i.e. true findings) among the set of elementary hypotheses as specified by \code{intervals}. The upper confidence limit is always equal to the size of the set.}

\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}

\examples{

#generate data, where the response Y is associated with certain groups of covariates
#namely cov 3-6, 9-12, 15-18
set.seed(1)
n=100
p=20
X <- matrix(rnorm(n*p),n,p)
beta <- c(rep(0,2),rep(1,4),rep(0,2),rep(1,4),rep(0,2),rep(1,4),rep(0,2))
Y <- X \%*\% beta + rnorm(n)

# Define the local test to be used in the closed testing procedure
mytest <- function(left,right)
{
  X <- X[,(left:right),drop=FALSE]
  lm.out <- lm(Y ~ X)
  x <- summary(lm.out)
  return(pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE))  
}

# perform the region procedure
reg <- regionmethod(rep(1,p), mytest, isadjusted=TRUE)
summary(reg)

#what are the smallest regions that are found to be significant? 
implications(reg)

#how many covariates within the full region of length 20 are at least associated with the response?
regionpick(reg, list(c(1,p)), alpha=0.05)

#how many covariates within the two subregions, (1,5) and (16,20) 
#are at least associated with the response?
regionpick(reg, list(c(1,5),c(16,20)))

}
