\name{hommel class}
\docType{class}
\alias{hommel}
\alias{hommel-class}
\alias{pvalue,hommel-method}

\title{Class "hommel" for storing the adjusted p-values of Hommel's multiple testing method.}

\description{The class hommel is the output of a call to \code{\link{hommelFast}}. It stores the adjusted p-values of all hypotheses.}

\section{Slots}{
  These slots are not meant to be directly accessed by the user.
  \describe{
    \item{\code{pvalues}:}{Object of class "numeric". Stores original (unadjusted) p-values for all hypotheses.}
    \item{\code{jvalues}:}{Object of class "numeric". Stores "jvalues", where each jvalue is the size of the largest intersection hypothesis that cannot be rejected within the closed testing procedure up to the corresponding alpha level that is stored in jumpalpha.}
    \item{\code{jumpalpha}:}{Object of class "numeric". Stores the corresponding alpha value for each jvalue.}
    \item{\code{adjusted}:}{Object of class "numeric". Stores adjusted p-values for all hypotheses.}
    \item{\code{simes}:}{Object of class "logical". Is TRUE if the adjusted p-values are calculated based on Simes' test. Is FALSE if they are based on Hommel's test.}
  }
}


\section{Methods}{
  \describe{
    \item{pvalue}{(hommel,indicator): Retrieves adjusted pvalues for all possible hypotheses (as specified by indicator) from the hommel object. If no indicator is given, all adjusted p-values are returned.}
  }
}

\author{Rosa Meijer: \email{r.j.meijer@lumc.nl}}

\seealso{\code{\link{pickSimes}}, \code{\link{curveSimes}}.}

\keyword{methods}

\examples{
  
  #Generate a vector of pvalues
  set.seed(1)
  n <- 1000
  pvalues <- c(runif(0.50*n,0,0.001), runif(0.50*n,0,1))
  
  #Create an hommel object in which the adjusted p-values are stored, based on a Simes'test
  #(or Hommel's test, by choosing simes = FALSE):
  hom <- hommelFast(pvalues, simes = TRUE)
  
  #To check how many of the hypotheses within a certain subset have to be false with 95\% 
  #confidence, the pickSimes function can be used. If no subset is selected, the whole set of
  #hypotheses is used.  
  pickSimes(hom, alpha=0.05, silent=FALSE)
  
  #The same can be done for any subset, for example the first 10 hypotheses: 
  pickSimes(hom, select=1:10, alpha=0.05, silent=FALSE)
  
  #A plot can be made that displays the results of pickSimes (i.e. the lower confidence bounds)
  #but for selecting only the first 1,2,3,.... hypotheses from the selected set.
  curveSimes(hom, select=1:10, alpha=0.05, plot=FALSE)
    
} 