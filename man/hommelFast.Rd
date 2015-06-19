\name{hommelFast}

\alias{hommelFast}

\title{Calculates adjusted p-values of Hommel's method efficiently.}

\description{Calculates adjusted p-values of Hommel's method efficiently.}

\usage{ hommelFast (pvalues, simes = TRUE) }

\arguments{
  \item{pvalues}{A vector of p-values.}
  \item{simes}{If \code{TRUE}, a Simes' test is used, if \code{FALSE} Hommel's test is used.}
}

\value{Returns a \code{\link{hommel}} object.}

\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}

\examples{
  
  #Generate a vector of pvalues
  set.seed(1)
  n <- 1000
  pvalues <- c(runif(0.50*n,0,0.001), runif(0.50*n,0,1))
  
  #Create an hommel object in which the adjusted p-values are stored, based on a Simes'test
  #(or Hommel's test, by choosing simes = FALSE):
  hom <- hommelFast(pvalues, simes = TRUE)
  
  #Retrieve the first 10 adjusted p-values by using the \code{pvalue} method. 
  pvalue(hom,1:10)
  
  
}
