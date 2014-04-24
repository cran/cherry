\name{Region class}
\docType{class}
\alias{region}
\alias{region-class}
\alias{show,region-method}
\alias{summary,region-method}
\alias{implications}
\alias{implications,region-method}
\alias{alpha,region-method}
\alias{pvalue}
\alias{pvalue,region-method}

\title{Class "region" for storing the result of the region procedure}

\description{The class region is the output of a call to \code{\link{regionmethod}}. It stores which region hypotheses have been rejected and can be used to calculate confidence sets for the number of true and/or false hypotheses among a selected region.}

\section{Slots}{
  These slots are not meant to be directly accessed by the user.
  \describe{
    \item{\code{weights}:}{Object of class "numeric". Stores the weights of elementary hypotheses.}
    \item{\code{isadjusted}:}{Object of class "logical". Stores whether adjusted p-values are calculated.}
    \item{\code{allpvalues}:}{Object of class "matrix". Stores (adjusted) p-values for all possible region hypotheses. Has value NA if adjusted p-value is larger than alpha. Corresponds to a 0x0 matrix if \code{all_pvalues} is set to FALSE (by default).}
    \item{\code{implications}:}{Object of class "matrix". Stores implications includinge (adjusted) p-values at chosen alpha-level}
    \item{\code{alpha}:}{Object of class "numeric". The type I error level chosen.}
    \item{\code{totalrejected}:}{Object of class "numeric". Stores the total number of rejected region hypotheses.}
  }
}


\section{Methods}{
  \describe{
    \item{show}{(region): Prints how many region hypotheses were tested and how many of them were rejected.}
    \item{summary}{(region): Prints the test results (as \code{show}).}
    \item{alpha}{(region): Retrieves the maximal alpha_value from the region object.}
    \item{implications}{(region,alpha): Retrieves the implications from a given region object. By default, the alpha-level is set to \code{alpha_max}, but the alpha level can be varied, if the supporting information is present in the region object.}
    \item{pvalue}{(region,left,right): Retrieves pvalues for all possible region hypotheses (indicated with left and right bound) from the region object. Only able to return the value if the allpvalues matrix is stored in the region object.}
  }
}

\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}

\seealso{\code{\link{regionmethod}}, \code{\link{regionpick}}.}

\keyword{methods}

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
reg <- regionmethod(rep(1,p), mytest, isadjusted=TRUE, all_pvalues=TRUE)
summary(reg)

#what are the smallest regions that are found to be significant? 
implications(reg)

#what are the smallest regions that are found to be significant at an alpha-level of 0.03? 
implications(reg, alpha=0.03)

#what is the adjusted p-value of the overall null-hypothesis 
#(corresponding to the region ranging from 1 to 20)?
pvalue(reg,1,20)

#how many covariates within the full region of length 20 are at least associated with the response?
regionpick(reg, list(c(1,p)), alpha=0.05)

#visualize the results by either plotting a polygon corresponding to the underlying graph
regionplot(reg)

#or by plotting the graph itself
regionplot2(reg)
} 