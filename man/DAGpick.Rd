\name{DAGpick}

\alias{DAGpick}

\title{Confidence limits for the number of false hypotheses in a given family of sets.}

\description{Calculates confidence limits for the number of false hypotheses on the basis of the DAG procedure within a family of sets.}

\usage{ DAGpick (DAG, indicators, optimization = "ILP")}
    
\arguments{
  \item{DAG}{DAG object, as returned by the function \code{DAGmethod}.}
  \item{indicators}{The names or indices of the sets (as specified in the DAGstructure object) for which you want to know the confidence limits. Note that, if there were duplicate sets in the original list, the index can be different from the one in the original list given to \code{construct}.}
  \item{optimization}{Can be, in ascending order of accuracy and computational costs: "LP" (linear programming) or "ILP" (integer linear programming).}
}
  
\value{The function \code{DAGpick} returns the lower bound of a 1-alpha confidence set for the number of false sets. For the moment, the function can only be used on DAG objects that correspond to a DAG with two-way logical relationships.}
  
\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}
  
\seealso{\code{\link{DAGmethod}}, \code{\link{DAG}}.}
  
\examples{

  #Generate data, where the response Y is associated with two (out of 4) covariates
  set.seed(1)
  n=100
  p=4
  X <- matrix(rnorm(n*p),n,p)
  beta <- c(0,0.5,0.5,0)
  Y <- X \%*\% beta + rnorm(n)
    
  # Let us assume we have the following sets that we want to test: 
  sets <- list(c(1,2,3,4), c(1,2), c(2,3,4), c(2,3), 1, 2, 3, 4)
  names(sets) <- c(1234, 12, 234, 23, 1, 2, 3, 4)
    
  # Start by making the corresponding graph structure
  struct <- construct(sets)

  # Define the local test to be used in the closed testing procedure. 
  # This test expects a set as input. 
  mytest <- function(set)
  { 
    X <- X[,set,drop=FALSE]
    lm.out <- lm(Y ~ X)
    x <- summary(lm.out)
    return(pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE))  
  }
    
  # Perform the DAG procedure (default is all-parents method).
  DAG <- DAGmethod(struct, mytest, isadjusted=TRUE)
  summary(DAG)
    
  # How many of the elementary hypotheses (the last 4 sets) have to be false 
  # with probability 1-alpha? 
  # Sets (don't have to be elementary hypotheses in general) must be specified 
  # by their index or name. 
  DAGpick(DAG, 5:8)
  DAGpick(DAG, c("1","2","3","4")) #as above, but while using names

}