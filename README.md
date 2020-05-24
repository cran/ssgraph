# **ssgraph** 
  
![](https://www.r-pkg.org/badges/version/ssgraph) ![](https://www.r-pkg.org/badges/last-release/ssgraph) ![](https://cranlogs.r-pkg.org/badges/ssgraph) 
![](https://cranlogs.r-pkg.org/badges/grand-total/ssgraph) 

The `R` package **ssgraph** is designed for Bayesian structure learning in graphical models using spike-and-slab priors. 
To speed up the computations, the computationally intensive tasks of the package are implemented in `C++` in parallel using **OpenMP**. 

## Installation

You can install the latest version from CRAN using:

``` r
install.packages( "ssgraph" )
```

``` r
require( "ssgraph" )
```

## Example

This is a simple example to see the preformance of the 

Frist, by using the function `bdgraph.sim` we simulate 60 observations (n = 60) from a multivariate
Gaussian distribution with 8 variables (p = 8) and “scale-free” graph structure, as follows:

``` r
data.sim = bdgraph.sim( n = 100, p = 8, graph = "scale-free", vis = TRUE )
round( head( data.sim $ data, 4 ), 2 )
```

Since the generated data are Gaussian, we run `ssgraph` function by choosing `method = "ggm"`, as follows:

``` r
ssgraph.obj <- ssgraph( data = data.sim, method = "ggm", iter = 5000, save = TRUE )

summary( ssgraph.obj )
```

To compare the result with true graph
``` r
compare( data.sim, ssgraph.obj, main = c( "Target", "ssgraph" ), vis = TRUE )
```
