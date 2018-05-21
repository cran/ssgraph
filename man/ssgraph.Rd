\name{ssgraph}
\alias{ssgraph}

\title{	Algorithm for graphical models using spike-and-slab priors  }
\description{
	This function has a sampling algorithm for Bayesian model determination in undirected graphical models, based on spike-and-slab priors. 
}
\usage{
ssgraph( data, n = NULL, method = "ggm", iter = 5000, burnin = iter / 2, 
         var1 = 4e-04, var2 = 1, lambda = 1, g.prior = 0.5, g.start = "full", 
         sig.start = NULL, save.all = FALSE, print = 1000, cores = 2 )
}

\arguments{
	\item{data}{
	There are two options: (1) an (\eqn{n \times p}) matrix or a \code{data.frame} corresponding to the data, 
	(2) an (\eqn{p \times p}) covariance matrix as \eqn{S=X'X} which \eqn{X} is the data matrix 
	(\eqn{n} is the sample size and \eqn{p} is the number of variables). 
	It also could be an object of class \code{"sim"}, from function \code{bdgraph.sim}.
	The input matrix is automatically identified by checking the symmetry.
	}

	\item{n}{The number of observations. It is needed if the \code{"data"} is a covariance matrix.}

	\item{method}{
		A character with two options \code{"ggm"} (default) and \code{"gcgm"}. 
		Option \code{"ggm"} is for Gaussian graphical models based on Gaussianity assumption.
		Option \code{"gcgm"} is for Gaussian copula graphical models for the data that not follow Gaussianity assumption 
		(e.g. continuous non-Gaussian, discrete, or mixed dataset).
	}
	  
	\item{iter}{The number of iteration for the sampling algorithm.}
	\item{burnin}{The number of burn-in iteration for the sampling algorithm.}
	
	\item{var1}{ Value for the variance of the the prior of precision matrix for the places that there is no link in the graph. }
	\item{var2}{ Value for the variance of the the prior of precision matrix for the places that there is link in the graph. }
	\item{lambda}{ Value for the parameter of diagonal element of the prior of precision matrix. }

	\item{g.prior}{
		For determining the prior distribution of each edge in the graph. 
		There are two options: a single value between \eqn{0} and \eqn{1} (e.g. \eqn{0.5} as a noninformative prior) 
		or an (\eqn{p \times p}) matrix with elements between \eqn{0} and \eqn{1}.
	}

	\item{g.start}{
		Corresponds to a starting point of the graph. It could be an (\eqn{p \times p}) matrix, \code{"empty"} (default), or \code{"full"}. 
		Option \code{"empty"} means the initial graph is an empty graph and \code{"full"} means a full graph. 
		It also could be an object with \code{S3} class \code{"ssgraph"} or \code{"bdgraph"}; 
		this option can be used to run the sampling algorithm from the last objects of previous run (see examples).     
	} 
	
	\item{sig.start}{
		Corresponds to a starting point of the covariance matrix. It must be positive definite matrix.     
	} 

	\item{save.all}{
		Logical: if FALSE (default), the adjacency matrices are NOT saved. 
		If TRUE, the adjacency matrices after burn-in are saved.
	}

	\item{print}{ Value to see the number of iteration for the MCMC algorithm. }	
	
	\item{cores}{ The number of cores to use for parallel execution. 
	   The default is to use \code{2} CPU cores of the computer.
	   The case \code{cores="all"} means all CPU cores to use for parallel execution. 
	   }	
}

\value{
	An object with \code{S3} class \code{"ssgraph"} is returned:
	
	\item{p_links}{ An upper triangular matrix which corresponds the estimated posterior probabilities of all possible links. }
	
	\item{K_hat}{ The posterior estimation of the precision matrix. }

	For the case "save.all = TRUE" is returned:

	\item{sample_graphs}{ A vector of strings which includes the adjacency matrices of visited graphs after burn-in.}
	\item{graph_weights}{ A vector which includes the counted numbers of visited graphs after burn-in. }

	\item{all_graphs}{A vector which includes the identity of the adjacency matrices for all iterations after burn-in. 
		  It is needed for monitoring the convergence of the MCMC sampling algorithm.}

	\item{all_weights}{A vector which includes the waiting times for all iterations after burn-in. 
	                   It is needed for monitoring the convergence of the MCMC sampling algorithm.}
}

\references{
Wang, H. (2015). Scaling it up: Stochastic search structure learning in graphical models, \emph{Bayesian Analysis}, 10(2):351-377

George, E.I. and McCulloch, R.E. (1993). Variable selection via Gibbs sampling. \emph{Journal of the American Statistical Association}, 88(423):881-889

Griffin, J.E. and Brown, P.J. (2010) Inference with normal-gamma prior distributions in regression problems. \emph{Bayesian Analysis}, 5(1):171-188

Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 
}

\author{ Reza Mohammadi }

\examples{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 80, p = 7, prob = 0.5, vis = TRUE )

# Running algorithm based on GGMs
ssgraph.obj <- ssgraph( data = data.sim, iter = 1000 )

summary( ssgraph.obj )

# To compare the result with true graph
compare( data.sim, ssgraph.obj, colnames = c("True graph", "ssgraph"), vis = TRUE )

\dontrun{

# Running algorithm with starting points from previous run
ssgraph.obj2 <- ssgraph( data = data.sim, iter = 5000, g.start = ssgraph.obj )
    
compare( data.sim, ssgraph.obj, ssgraph.obj2, 
         colnames = c( "True graph", "Frist run", "Second run" ) )
}
}
