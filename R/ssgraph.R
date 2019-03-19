## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2018 - 2019 Reza Mohammadi                                                     |
#                                                                                                  |
#     This file is part of ssgraph package.                                                        |
#                                                                                                  |
#     ssgraph is free software: you can redistribute it and/or modify it under                     |
#     the terms of the GNU General Public License as published by the Free                         |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                  |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  R code for Graphcial models based on spike and slab priors                                      |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

ssgraph = function( data, n = NULL, method = "ggm", not.cont = NULL, iter = 5000, 
                    burnin = iter / 2, var1 = 4e-04, var2 = 1, lambda = 1, g.prior = 0.5, 
                    g.start = "full", sig.start = NULL, save = FALSE, print = 1000, 
                    cores = NULL )
{
    if( iter < burnin ) stop( " Number of iteration must be more than number of burn-in" )
    if(  var1 <= 0    ) stop( " 'var1' must be more than 0" )
    if(  var2 <= 0    ) stop( " 'var2' must be more than 0" )
    
    num_machine_cores = BDgraph::detect_cores()
    if( is.null( cores ) ) cores = num_machine_cores - 1
    if( cores == "all" )   cores = num_machine_cores

    .C( "omp_set_num_cores", as.integer( cores ), PACKAGE = "ssgraph" )
    
    burnin <- floor( burnin )
    
    if( class( data ) == "sim" )
    {
        not.cont <- data $ not.cont
        data     <- data $ data
    }
    
    if( !is.matrix( data ) & !is.data.frame( data ) ) stop( " Data must be a matrix or dataframe" )
    if( is.data.frame( data ) ) data <- data.matrix( data )
    
    if( any( is.na( data ) ) ) 
    {
        if( method == "ggm" ) stop( " 'ggm' method does not deal with missing values. You could choose option method = gcgm" )	
        gcgm_NA = 1
    }else{
        gcgm_NA = 0
    }
    
    p <- ncol( data )
    if( p < 3 ) stop( " Number of variables/nodes ('p') must be more than or equal with 2" )
    if( is.null( n ) ) n <- nrow( data )

    if( is.data.frame( g.prior ) ) g.prior <- data.matrix( g.prior )
    if( class( g.prior ) == "dtCMatrix" ) g.prior = as.matrix( g.prior )
    if( ( class( g.prior ) == "bdgraph" ) | ( class( g.prior ) == "ssgraph" ) ) g.prior <- as.matrix( BDgraph::plinks( g.prior ) )

    if( !is.matrix( g.prior ) )
    {
        if( ( g.prior <= 0 ) | ( g.prior >= 1 ) ) stop( " 'g.prior' must be between 0 and 1" )
        g.prior = matrix( g.prior, p, p )
    }else{
        if( ( nrow( g.prior ) != p ) | ( ncol( g.prior ) != p ) ) stop( " 'g.prior' and 'data' have non-conforming size" )
        if( any( g.prior < 0 ) || any( g.prior > 1 ) ) stop( " Element of 'g.prior', as a matrix, must be between 0 and 1" )
    }
    
    g_prior = g.prior
    g_prior[ lower.tri( g_prior, diag = TRUE ) ] <- 0
    g_prior = g_prior + t( g_prior )
    
    if( method == "gcgm" )
    {
        if( isSymmetric( data ) ) stop( " method='gcgm' requires all data" )
        
        if( is.null( not.cont ) )
        {
	        not.cont = c( rep( 1, p ) )
	        for( j in 1:p )
	            if( length( unique( data[ , j ] ) ) > min( n / 2 ) ) not.cont[ j ] = 0
        }else{
            if( !is.vector( not.cont )  ) stop( " 'not.cont' must be a vector with length of number of variables" )
            if( length( not.cont ) != p ) stop( " 'not.cont' must be a vector with length of number of variables" )
            if( ( sum( not.cont == 0 ) + sum( not.cont == 1 ) ) != p ) stop( " Element of 'not.cont', as a vector, must be 0 or 1" )
        }
        
        R <- 0 * data
        for( j in 1:p )
            if( not.cont[ j ] )
                R[ , j ] = match( data[ , j ], sort( unique( data[ , j ] ) ) ) 
        R[ is.na( R ) ] = 0     # dealing with missing values	
        
        # copula for continuous non-Gaussian data
        if( gcgm_NA == 0 && min( apply( R, 2, max ) ) > ( n - 5 * n / 100 ) )
        {
            # copula transfer 
            data = stats::qnorm( apply( data, 2, rank ) / ( n + 1 ) )
            data = t( ( t( data ) - apply( data, 2, mean ) ) / apply( data, 2, stats::sd ) )
            
            method = "ggm"
        }else{	
            # for non-Gaussian data
            Z                  <- stats::qnorm( apply( data, 2, rank, ties.method = "random" ) / ( n + 1 ) )
            Zfill              <- matrix( stats::rnorm( n * p ), n, p )   # for missing values
            Z[ is.na( data ) ] <- Zfill[ is.na( data ) ]                  # for missing values
            Z                  <- t( ( t( Z ) - apply( Z, 2, mean ) ) / apply( Z, 2, stats::sd ) )
            S                  <- t( Z ) %*% Z
        }
    } 
    
    if( method == "ggm" ) 
    {
        if( isSymmetric( data ) )
        {
            if ( is.null( n ) ) stop( " Please specify the number of observations 'n'" )
            cat( " Input is identified as the covriance matrix. \n" )
            S <- data
        }else{
            S <- t( data ) %*% data
        }
    }
    
    if( ( class( g.start ) == "bdgraph" ) | ( class( g.start ) == "ssgraph" ) ) 
    {
        G <- g.start $ last_graph
        K <- g.start $ last_K
    } 

    if( class( g.start ) == "sim" ) 
    {
        G <- unclass( g.start $ G )
        K <- g.start $ K
    }
    
    if( class( g.start ) == "graph" ) G <- unclass( g.start )
    
    if( ( class( g.start ) == "character" ) && ( g.start == "empty" ) ) G = matrix( 0, p, p )
    if( ( class( g.start ) == "character" ) && ( g.start == "full"  ) ) G = matrix( 1, p, p )
    if( is.matrix( g.start ) ) 
    {
        if( ( sum( g.start == 0 ) + sum( g.start == 1 ) ) != ( p ^ 2 ) ) stop( " Element of 'g.start', as a matrix, must be 0 or 1" )
        G = g.start
    }
    
    if( ( nrow( G ) != p ) | ( ncol( G ) != p ) ) stop( " 'g.start' and 'data' have non-conforming size" )

    diag( G ) = 0    
    if( !isSymmetric( G ) )
    {
        G[ lower.tri( G, diag( TRUE ) ) ] <- 0
        G  = G + t( G )
    }

    if( ( class( g.start ) != "sim" ) & ( class( g.start ) != "bdgraph" ) & ( class( g.start ) != "ssgraph" ) )
    {
        if( is.null( sig.start ) ) sigma = S else sigma = sig.start
        K = solve( sigma )      # precision or concentration matrix (omega)
    }
    
    if( save == TRUE )
    {
        qp1           = ( p * ( p - 1 ) / 2 ) + 1
        string_g      = paste( c( rep( 0, qp1 ) ), collapse = '' )
        sample_graphs = c( rep ( string_g, iter - burnin ) )  # vector of numbers like "10100" 
        graph_weights = c( rep ( 0, iter - burnin ) )         # waiting time for every state
        all_graphs    = c( rep ( 0, iter - burnin ) )         # vector of numbers like "10100"
        all_weights   = c( rep ( 1, iter - burnin ) )         # waiting time for every state		
        size_sample_g = 0
    }
    
    if( ( save == TRUE ) && ( p > 50 & iter > 20000 ) )
    {
        cat( "  WARNING: Memory needs to run this function is around " )
        print( ( iter - burnin ) * utils::object.size( string_g ), units = "auto" ) 
    } 
    
    nmc     = iter - burnin
    p_links = matrix( 0, p, p )
    K_hat   = matrix( 0, p, p )

    mes <- paste( c( iter, " iteration is started.                    " ), collapse = "" )
    cat( mes, "\r" )
    
## - -  main BDMCMC algorithms implemented in C++ - - - - - - - - - - - - - - - - - - - - - - - - -|
    if( save == FALSE )
    { 
        if( method == "ggm" )
        {
             result = .C( "ggm_spike_slab_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), K = as.double(K), as.double(S), as.integer(p), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links), as.integer(n),
                         as.double(var1), as.double(var2), as.double(lambda), as.double(g_prior), as.integer(print), PACKAGE = "ssgraph" )
        }
        
        if( method == "gcgm" )
        {
            not_continuous = not.cont
            
             result = .C( "gcgm_spike_slab_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), K = as.double(K), as.double(S), as.integer(p), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links), as.integer(n),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(gcgm_NA),
                         as.double(var1), as.double(var2), as.double(lambda), as.double(g_prior), as.integer(print), PACKAGE = "ssgraph" )
        }
    }else{
        if( method == "ggm" )
        {
            result = .C( "ggm_spike_slab_map", as.integer(iter), as.integer(burnin), G = as.integer(G), K = as.double(K), as.double(S), as.integer(p), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links), as.integer(n),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.double(var1), as.double(var2), as.double(lambda), as.double(g_prior), as.integer(print), PACKAGE = "ssgraph" )
        }
        
        if( method == "gcgm" )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_spike_slab_map", as.integer(iter), as.integer(burnin), G = as.integer(G), K = as.double(K), as.double(S), as.integer(p), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links), as.integer(n),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(gcgm_NA),
                         as.double(var1), as.double(var2), as.double(lambda), as.double(g_prior), as.integer(print), PACKAGE = "ssgraph" )
        }
    }
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
    
    label      = colnames( data )
    p_links    = matrix( result $ p_links, p, p, dimnames = list( label, label ) ) 
    K_hat      = matrix( result $ K_hat  , p, p, dimnames = list( label, label ) ) 
    
    last_graph = matrix( result $ G      , p, p, dimnames = list( label, label ) )
    last_K     = matrix( result $ K      , p, p, dimnames = list( label, label ) )

    p_links[ lower.tri( p_links, diag = TRUE ) ] = 0
    p_links = p_links / nmc
    K_hat   = K_hat / nmc

    if( save == TRUE )
    {
        size_sample_g = result $ size_sample_g
        sample_graphs = result $ sample_graphs[ 1 : size_sample_g ]
        graph_weights = result $ graph_weights[ 1 : size_sample_g ]
        all_graphs    = result $ all_graphs + 1
        all_weights   = result $ all_weights

        output = list( p_links = p_links, K_hat = K_hat, last_graph = last_graph, last_K = last_K,
                       sample_graphs = sample_graphs, graph_weights = graph_weights, 
                       all_graphs = all_graphs, all_weights = all_weights )
    }else{
        output = list( p_links = p_links, K_hat = K_hat, last_graph = last_graph, last_K = last_K )
    }
    
    class( output ) = "ssgraph"
    return( output )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Summary for the ssgraph boject                                                                |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
summary.ssgraph = function( object, round = 2, vis = TRUE, ... )
{
    p_links    = object $ p_links
    selected_g = BDgraph::select( p_links, cut = 0.5 )   
    
    if( vis == TRUE )
    {
        if( !is.null( object $ graph_weights ) ) 
            op = graphics::par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) ) 

        p        = ncol( p_links )
        subGraph = "Selected graph with edge posterior probability = 0.5"
        if( p < 20 ) size = 15 else size = 2
        
        # - - - plot selected graph
        G  <- igraph::graph.adjacency( selected_g, mode = "undirected", diag = FALSE )
        igraph::plot.igraph( G, layout = igraph::layout.circle, main = "Selected graph", sub = subGraph, vertex.color = "white", vertex.size = size, vertex.label.color = 'black' )
        
        if( !is.null( object $ graph_weights ) )
        {
            sample_graphs = object $ sample_graphs
            graph_weights = object $ graph_weights
            sum_gWeights  = sum( graph_weights )

            # - - - plot posterior distribution of graph
            graph_prob = graph_weights / sum_gWeights
            graphics::plot( x = 1 : length( graph_weights ), y = graph_prob, type = "h", main = "Posterior probability of graphs",
                            ylab = "Pr(graph|data)", xlab = "graph", ylim = c( 0, max( graph_prob ) ) )

            # - - - plot posterior distribution of graph size
            sizesample_graphs = sapply( sample_graphs, function( x ) length( which( unlist( strsplit( as.character( x ), "" ) ) == 1 ) ) )
            xx       <- unique( sizesample_graphs )
            weightsg <- vector()
            
            for( i in 1 : length( xx ) ) weightsg[ i ] <- sum( graph_weights[ which( sizesample_graphs == xx[ i ] ) ] )
            
            graphics::plot( x = xx, y = weightsg / sum_gWeights, type = "h", main = "Posterior probability of graphs size", ylab = "Pr(graph size|data)", xlab = "Graph size" )
            
            # - - - plot trace of graph size
            all_graphs     = object $ all_graphs
            sizeall_graphs = sizesample_graphs[ all_graphs ]
            
            graphics::plot( x = 1 : length( all_graphs ), sizeall_graphs, type = "l", main = "Trace of graph size", ylab = "Graph size", xlab = "Iteration" )

            graphics::par( op )
        }
    }

    return( list( selected_g = selected_g, p_links = round( p_links, round ), K_hat = round( object $ K_hat, round ) ) )
}  
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Plot for the ssgraph boject                                                                   |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
plot.ssgraph = function( x, cut = 0.5, layout = layout.circle, ... )
{
    if( ( cut < 0 ) || ( cut > 1 ) ) stop( " Value of 'cut' must be between 0 and 1." )
    
    selected_g = BDgraph::select( x, cut = cut )
    
    G = igraph::graph.adjacency( selected_g, mode = "undirected", diag = FALSE )
    igraph::plot.igraph( G, layout = layout, sub = paste0( "Edge posterior probability = ", cut ), ... )	   		
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Print for the ssgraph boject                                                                  |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
print.ssgraph = function( x, ... )
{
    p_links = x $ p_links
    selected_g = BDgraph::select( p_links, cut = 0.5 )
    
    cat( paste( "\n Adjacency matrix of selected graph \n" ), fill = TRUE )
    print( selected_g )
    
    cat( paste( "\n Edge posterior probability of the links \n" ), fill = TRUE )
    print( round( p_links, 2 ) )
} 
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |



  
