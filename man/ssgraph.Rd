\name{ssgraph}
\alias{ssgraph}

\title{	Algorithm for graphical models using spike-and-slab priors  }
\description{
	This function has a sampling algorithm for Bayesian model determination in undirected graphical models, based on sbike-and-slab prios. 
}
\usage{
ssgraph( data, n = NULL, iter = 5000, burnin = iter / 2, v1 = NULL, 
         v2 = NULL, lambda = 1, g.prior = 0.5, g.start = "full", 
         sig.start = NULL, print = 100 )
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
	  
	\item{iter}{The number of iteration for the sampling algorithm.}
	\item{burnin}{The number of burn-in iteration for the sampling algorithm.}
	
	\item{v1}{ Variance for the the prior of precision matrix for the places that there is no link in the graph. }
	\item{v2}{ Variance for the the prior of precision matrix for the places that there is link in the graph. }
	\item{lambda}{ Value for the diagonal element of the prior of precision matrix. }

	\item{g.prior}{
		For determining the prior distribution of each edge in the graph. 
		There are two options: a single value between \eqn{0} and \eqn{1} (e.g. \eqn{0.5} as a noninformative prior) 
		or an (\eqn{p \times p}) matrix with elements between \eqn{0} and \eqn{1}.
	}

	\item{g.start}{
		Corresponds to a starting point of the graph. It could be \code{"empty"} and \code{"full"}. 
		Option \code{"empty"} means the initial graph is an empty graph and \code{"full"} (default) means a full graph. 
		It also could be an object with \code{S3} class \code{"bdgraph"}; 
		with this option we could run the sampling algorithm from the last objects of previous run.     
	} 
	
	\item{sig.start}{
		Corresponds to a starting point of the covariance matrix. It must be positive definite matrix.     
	} 

	\item{print}{ Value to see the number of iteration for the MCMC algorithm. }	
}

\value{
	An object with \code{S3} class \code{"ssgraph"} is returned:
	
	\item{p_links}{ An upper triangular matrix which corresponds the estimated posterior probabilities of all possible links. }
	
	\item{K_hat}{ The posterior estimation of the precision matrix. }
}

\references{
Wang, H. (2015). Scaling it up: Stochastic search structure learning in graphical models, \emph{Bayesian Analysis}, 10(2)351-377

Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138
}

\author{ Reza Mohammadi }

\examples{
p <- 6
n <- 3 * p

K      <- toeplitz( c( 1, 0.5, rep( 0, p - 2 ) ) )
K[1,p] <- 0.4
K[p,1] <- 0.4

S      <- n * solve( K )

out_ss <- ssgraph( data = S, n = n, iter = 1000 )

round( out_ss $ p_links, 2 )
round( out_ss $ K_hat, 2 )
}
