\name{predict.ssgraph}

\alias{predict.ssgraph}

\title{ Predict function for \code{S3} class "\code{bdgraph}" }

\description{
Provides predict values of the results for function \code{\link{ssgraph}}.  
}

\usage{
\method{predict}{ssgraph}( object, iter = 1, ... )
}

\arguments{
  \item{object}{ object of \code{S3} class "\code{ssgraph}", from function \code{\link{ssgraph}}. }
  \item{iter}{ number of predictions. }

  \item{\dots}{ additional parameters. }
}

\references{
Wang, H. (2015). Scaling it up: Stochastic search structure learning in graphical models, \emph{Bayesian Analysis}, 10(2):351-377

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{ssgraph}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
   
ssgraph.obj <- ssgraph( data = data.sim )
   
predict( ssgraph.obj, iter = 20 )
}
}
