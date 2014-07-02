\name{DAG class}
\docType{class}
\alias{DAG}
\alias{DAG-class}
\alias{show,DAG-method}
\alias{summary,DAG-method}
\alias{implications,DAG-method}
\alias{alpha,DAG-method}
\alias{pvalue,DAG-method}

\title{Class "DAG" for storing the result of the DAG multiple testing method}

\description{The class DAG is the output of a call to \code{\link{DAGmethod}}. It stores which hypotheses have been rejected and can be used to calculate confidence sets for the number of true and/or false hypotheses among a selected set of hypotheses.}

\section{Slots}{
  These slots are not meant to be directly accessed by the user.
  \describe{
    \item{\code{sets}:}{Object of class "list". Stores unique original sets that are to be tested.}
    \item{\code{method}:}{Object of class "character". Stores whether the any-parent or all-parents method has been used.}
    \item{\code{isadjusted}:}{Object of class "logical". Stores whether adjusted p-values are calculated.}
    \item{\code{allpvalues}:}{Object of class "numeric". Stores (adjusted) p-values for all hypotheses. Has value NA if adjusted p-value is larger than alpha.}
    \item{\code{implications}:}{Object of class "logical". Stores whether hypotheses are implications at chosen alpha-level}
    \item{\code{alpha}:}{Object of class "numeric". The type I error level chosen.}
    \item{\code{rejected}:}{Object of class "logical". Stores for each hypothesis wether this hypothesis has been rejected}
    \item{\code{leaf_based_sets}:}{Object of class "list". Stores sets expressed in the indices of their corresponding leaf nodes (which are sets itself).}
    \item{\code{twoway}:}{Object of class "logical". Is TRUE if the final DAG structure has twoway logical relationships.}
  }
}
    

\section{Methods}{
  \describe{
    \item{show}{(DAG): Prints how many hypotheses there are in total and how many of them were rejected.}
    \item{summary}{(DAG): Prints the test results (as \code{show}).}
    \item{alpha}{(DAG): Retrieves the maximal alpha_value from the DAG object.}
    \item{implications}{(DAG): Retrieves the implications from a given DAG object.}
    \item{pvalue}{(DAG,indicator): Retrieves pvalues for all possible hypotheses (as specified by indicator) from the DAG object.}
  }
}

\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}

\seealso{\code{\link{DAGmethod}}, \code{\link{DAGpick}}.}

\keyword{methods}

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
  
  # What are the smallest sets that are found to be significant? If the sets have names, 
  # as in our example, the \code{implications} function will return the names 
  # of the implying sets, together with their (adjusted) p-value. 
  # If no names are provided, indices will be used instead of the names. 
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