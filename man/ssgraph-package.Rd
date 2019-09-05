\name{ssgraph-package}
\alias{ssgraph-package}

\alias{ copula_C }
\alias{ copula_S }

\docType{package}

\title{ Bayesian Graphical Estimation using Spike-and-Slab Priors }

\description{
The \code{R} package \pkg{ssgraph} is for Bayesian estimation of graphical models by using spike-and-slab priors.
The package is implemented the recent improvements in the Bayesian graphical models' literature, including Wang (2015).
To speed up the computations, the computationally intensive tasks of the package are implemented in \code{C++} in parallel using \pkg{OpenMP}. 
}

\section{How to cite this package}{
To cite \pkg{ssgraph} in publications use:
\preformatted{
Mohammadi R. (2019). ssgraph: Bayesian Graphical Estimation using 
Spike-and-Slab Priors, R package version 1.9, 
https://cran.r-project.org/package=ssgraph
}
}

\author{Reza Mohammadi \email{a.mohammadi@uva.nl} }

\references{
Wang, H. (2015). Scaling it up: Stochastic search structure learning in graphical models, \emph{Bayesian Analysis}, 10(2):351-377

George, E. I. and McCulloch, R. E. (1993). Variable selection via Gibbs sampling. \emph{Journal of the American Statistical Association}, 88(423):881-889

Griffin, J. E. and Brown, P. J. (2010) Inference with normal-gamma prior distributions in regression problems. \emph{Bayesian Analysis}, 5(1):171-188

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645

Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138
}

\examples{
\dontrun{

library( ssgraph )

# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 100, p = 8, size = 10, vis = TRUE )

# Running algorithm based on GGMs
ssgraph.obj <- ssgraph( data = data.sim, iter = 5000, save = TRUE )

summary( ssgraph.obj )

# To compare the result with true graph
compare( data.sim, ssgraph.obj, main = c( "Target", "ssgraph" ), vis = TRUE )
}
}

\keyword{package}
