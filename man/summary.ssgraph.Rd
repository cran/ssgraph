\name{summary.ssgraph}
\alias{summary.ssgraph}

\title{ Summary function for \code{S3} class \code{"ssgraph"} }

\description{
Provides a summary of the results for function \code{\link{ssgraph}}.  
}

\usage{
\method{summary}{ssgraph}( object, round = 2, vis = TRUE, ... )
}

\arguments{
  \item{object}{ An object of \code{S3} class \code{"ssgraph"}, from function \code{\link{ssgraph}}. }
  \item{round}{ A value for rounding all probabilities to the specified number of decimal places. }
  \item{vis}{ Visualize the results. }
  \item{\dots}{System reserved (no specific usage).}
}

\value{
	\item{selected_g}{The adjacency matrix corresponding to the selected graph which has the highest posterior probability.}
	\item{p_links}{An upper triangular matrix corresponding to the posterior probabilities of all possible links.}
	\item{K_hat}{The estimated precision matrix.}
}

\references{
Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and Wit, E. C. (2017). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108v5} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{ssgraph}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
   
ssgraph.obj <- ssgraph( data = data.sim, save = TRUE )
   
summary( ssgraph.obj )
}
}
