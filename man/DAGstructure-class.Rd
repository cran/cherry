\name{DAGstructure class}
\docType{class}
\alias{DAGstructure}
\alias{DAGstructure-class}
\alias{istwoway}
\alias{istwoway,DAGstructure-method}

\title{Class "DAGstructure" for storing the result of the construct method that constructs the DAG.}

\description{The class DAGstructure is the output of a call to \code{\link{construct}}. It stores the DAG structure that is induced by the given sets.}

\section{Slots}{
  These slots are not meant to be directly accessed by the user.
  \describe{
    \item{\code{parents}:}{Object of class "list". Stores the parents of each set, indicated by indices.}
    \item{\code{children}:}{Object of class "list". Stores the children of each set, indicated by indices.}
    \item{\code{sets}:}{Object of class "list". Stores the sets that are used for the DAG construction.}
    \item{\code{twoway}:}{Object of class "logical". Is TRUE if the final DAG structure has twoway logical relationships.
    }
 }
}


\section{Methods}{
  \describe{
    \item{istwoway}{(DAGstructure): Indicates whether given DAGstructure has twoway relationships.}
  }
}

\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}

\seealso{\code{\link{DAGmethod}}}

\keyword{methods}

\examples{
  
  # Let us assume we have the following sets that we want to test: 
  sets <- list(c(1,2,3,4), c(1,2), c(2,3,4), c(2,3), 1, 2, 3, 4)
  # The sets need to have names in order to be able to look up their p-values later
  names(sets) <- c(1234, 12, 234, 23, 1, 2, 3, 4)
  
  # Start by making the corresponding graph structure
  struct <- construct(sets)
  
  # Check whether the DAG has toway logical relations: 
  istwoway(struct)
  
} 