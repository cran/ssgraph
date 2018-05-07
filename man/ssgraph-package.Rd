\name{ssgraph-package}
\alias{ssgraph-package}

\docType{package}

\title{ Bayesian Graphical Estimation using Spike-and-Slab Priors }

\description{
The \code{R} package \pkg{ssgraph} provides statistical tools for Bayesian estimation in graphical models by using spike-and-slab priors.
The package is implemented the recent improvements in the Bayesian graphical models' literature, including Wang (2015).
To speed up the computations, the computationally intensive tasks of the package are implemented in \code{C++} in parallel using \pkg{OpenMP}. 
}

\section{How to cite this package}{
Whenever using this package, please cite as
\preformatted{
Mohammadi R. (2018). ssgraph: Bayesian Graphical Estimation using Spike-and-Slab Priors, 
R package version 1.1, https://cran.r-project.org/package=ssgraph
}
}

\author{Reza Mohammadi <a.mohammadi@uva.nl>}

\references{
Wang, H. (2015). Scaling it up: Stochastic search structure learning in graphical models, \emph{Bayesian Analysis}, 10(2):351-377

George, E.I. and McCulloch, R.E. (1993). Variable selection via Gibbs sampling. \emph{Journal of the American Statistical Association}, 88(423):881-889

Griffin, J.E. and Brown, P.J. (2010) Inference with normal-gamma prior distributions in regression problems. \emph{Bayesian Analysis}, 5(1):171-188

Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 
}

\keyword{ Package, ssgraph, Bayesian Structure Learning, Gaussian Graphical Models, Gaussian Copula Graphical Models, Spike and Slab Priors }
