\name{Plot region objects}

\alias{regionplot}
\alias{regionplot2}

\title{Visualizing of the region hypotheses that could be rejected.}

\description{Visualizes region objects as created through a call to \code{\link{regionmethod}}.}

\usage{
  regionplot (region, alpha, color="red")
  
  regionplot2 (region, alpha, color_rej="red", color_unrej="grey") 
}

\arguments{
  \item{region}{An object of class \code{\link{region}}, typically created through a call to \code{\link{regionmethod}}.}
  \item{alpha}{For region objects with adjusted p-values, specifies the value of alpha for which rejections should be plotted (optional).}
  \item{color}{Color that is used to indicate rejected region hypotheses.}
  \item{color_rej}{Color that is used to indicate rejected region hypotheses.}
  \item{color_unrej}{Color that is used to indicate unrejected region hypotheses.}
}

\details{Both plot functions create a graph that visualizes all possible region hypotheses. Each region hypothesis is a node in the graph, and from each region hypothesis two edged connect the hypothesis with its child hypotheses. The \code{regionplot2} function visualized the graph with its nodes and edges. This function is especially useful for \code{\link{region}} objects with a limited number of elementary hypotheses. The \code{regionplot} function does not display the nodes and edges separately, but draws a polygon that follows the original graph structure. 
}

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

#visualize the results by either plotting a polygon corresponding to the underlying graph
regionplot(reg)

#or by plotting the graph itself
regionplot2(reg)


}
