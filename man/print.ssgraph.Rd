\name{print.ssgraph}
\alias{print.ssgraph}

\title{Print function for \code{S3} class \code{"ssgraph"} }

\description{
Prints the information about the selected graph which could be a graph with links for which their estimated posterior probabilities are greater than 0.5
or graph with the highest posterior probability.
It provides adjacency matrix, size and posterior probability of the selected graph.
}

\usage{\method{print}{ssgraph}( x, round = 2, ... )}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"ssgraph"}, from function \code{\link{ssgraph}}. }
  \item{round}{A value to round the probabilities to the specified number of decimal places.}
 
  \item{\dots}{System reserved (no specific usage). }
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 
}

\author{ Reza Mohammadi }

\seealso{\code{\link{ssgraph}}}

\examples{
\dontrun{
require( "BDgraph" )

# Generating multivariate normal data from a 'random' graph
dat <- bdgraph.sim( n = 50, p = 6, vis = TRUE )
  
out_ss <- ssgraph( data = dat )
   
print( out_ss )
}
}