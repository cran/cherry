\name{Region Testing}

\alias{regionmethod}

\title{Testing of all possible region hypotheses}

\description{Tests all possible region hypotheses in a given region while controlling the FWER, using a user-specified local test.}

\usage{ 
regionmethod (weights, test, alpha_max = 0.05, all_pvalues = FALSE, 
              isadjusted = FALSE, verbose = FALSE) }

\arguments{
  \item{weights}{a vector that indicates the weight each elementary hypotheses should receive in the multiple testing procedure. The length of the vector should equal the number of elementary hypotheses. All values should be strictly positive.}
  \item{test}{A function that performs the local test. The function should have the left and rightbound of the region as input (as two separate numbers), and return a p-value.}
  \item{alpha_max}{The significance level of the test procedure.}
  \item{all_pvalues}{If set to TRUE, the procedure will return a matrix with the p-values of all tested region hypotheses. If set to FALSE, only the p-values of the implications will be returned.}
  \item{isadjusted}{If set to TRUE, adjusted p-values will be calculated. Otherwise, the p-values of all rejected hypotheses will equal alpha_max.}
  \item{verbose}{If set to TRUE, while running the method, a counter will indicate how many region hypotheses are already rejected.}
}

\details{The function \code{regionmethod} tests all possible region hypotheses within one main interval, while controlling the familywise error rate.}

\note{The number of region hypotheses is quadratic in the number of elementary hypotheses. The number of elementary hypotheses should for computational reasons not be too large. For values in between 1000 and 10.000 the region procedure can still be used, but one might consider calculating the individual raw p-values beforehand.
      
      If both \code{isadjusted} and \code{all_pvalues} are set to TRUE, afterwards implications can be obtained for all alpha-values smaller or equal to \code{alpha_max}.}

\value{The function \code{regionmethod} returns an object of class \code{\link{region}}.}

\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}

\references{
  Meijer, Krebs and Goeman (2015) Statistical Applications in Genetics and Molecular Biology 14 (1) 1-19.
}

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

#at least how many covariates within the full region of length 20 
#are associated with the response?
regionpick(reg, list(c(1,p)), alpha=0.05)

#visualize the results by either plotting a polygon corresponding to the underlying graph
regionplot(reg)

#or by plotting the graph itself
regionplot2(reg)

}