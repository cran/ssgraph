\name{plot.ssgraph}
\alias{plot.ssgraph}

\title{ Plot function for \code{S3} class \code{"ssgraph"} }

\description{ Visualizes structure of the selected graphs which could be a graph with links for which their estimated posterior probabilities are greater than 0.5
or graph with the highest posterior probability. }

\usage{ \method{plot}{ssgraph}( x, cut = 0.5, layout = layout.circle, ... ) }

\arguments{
  \item{x}{ An object of \code{S3} class \code{"ssgraph"}, from function \code{\link{ssgraph}}. }
  
  \item{cut}{ Threshold for including the links in the selected graph based on the estimated posterior probabilities of the links; See the examples. } 

  \item{layout}{The vertex placement algorithm which is according to \pkg{igraph} package. }
  
  \item{\dots}{System reserved (no specific usage).}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 
}

\author{ Reza Mohammadi }

\seealso{\code{\link{ssgraph}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
dat <- bdgraph.sim( n = 50, p = 6, vis = TRUE )
  
out_ss <- ssgraph( data = dat )

plot( out_ss )

plot( out_ss, cut = 0.3 )
}
}
