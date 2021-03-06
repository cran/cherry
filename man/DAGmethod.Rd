\name{DAG Testing}

\alias{DAGmethod}

\title{Testing of hypotheses, forming a DAG.}

\description{Tests all hypotheses in a given DAG while controlling the FWER, using a user-specified local test.}

\usage{ DAGmethod (DAGstructure, test, alpha_max = 0.05, method = "all", isadjusted = FALSE, 
        optimization = "none", degree = "group", pvalues = NULL, verbose = FALSE)}

\arguments{
  \item{DAGstructure}{DAGstructure object, as returned by the function \code{construct}.}
  \item{test}{A function that performs the local test. The function should have a set as input and return a p-value.}
  \item{alpha_max}{The significance level of the test procedure.}
  \item{method}{Type of DAG procedure that is chosen. "all" gives the all-parents procedure, "any" the any-parent procedure.}
  \item{isadjusted}{If set to TRUE, adjusted p-values will be calculated. Otherwise, the p-values of all rejected hypotheses will equal alpha_max.}
  \item{optimization}{Can be, in ascending order of accuracy and computational costs: "none", "LP" (linear programming) or "ILP" (integer linear programming).}
  \item{degree}{Can be "group" or "individual". If "group" is chosen, optimization is done on a group level, otherwise on an individual level (more accurate, but more time-consuming).}
  \item{pvalues}{Optional (in case of stored p-values): a vector in which the raw p-values of the exact sets as found in the DAGstructure argument are stored (in the same order). If the test function is provided, this argument is not necessary.}
  \item{verbose}{If set to TRUE, while running the method, a counter will indicate how many hypotheses are already rejected.}
}

\details{The function \code{DAGmethod} tests all possible hypotheses within a given DAG structure, while controlling the familywise error rate.}

\value{The function \code{DAGmethod} returns an object of class \code{\link{DAG}}.}

\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}

\references{
  Meijer and Goeman (2015) Biometrical Journal 57 (1) 123-143.
}

\seealso{\code{\link{DAG}}, \code{\link{DAGstructure}}, \code{\link{construct}}, \code{\link{DAGpick}}.}

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
  
  # Check whether the DAG has toway logical relations: 
  istwoway(struct)
  
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
  
  # What are the smallest sets that are found to be significant? 
  implications(DAG)
  
  # What is the adjusted p-value of the null-hypothesis corresponding to the fourth set, 
  # which is set c(2,3)? 
  # To look up the pvalue, the function uses the index or name of the set 
  # in the list of sets stored in the DAGstructure. 
  # (Note that, if there were duplicate sets in the original list, this index can be different from 
  # the one in the original list given to \code{construct})
  pvalue(DAG,4)
  pvalue(DAG, "23") #as above, but while using names
  
  # How many of the elementary hypotheses (the last 4 sets) have to be false 
  # with probability 1-alpha? 
  # Sets (don't have to be elementary hypotheses in general) must be specified 
  # by their index or name. 
  DAGpick(DAG, 5:8)
  DAGpick(DAG, c("1","2","3","4")) #as above, but while using names
  
}