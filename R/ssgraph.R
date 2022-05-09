## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2022  Reza Mohammadi                                |
#                                                                              |
#     This file is part of ssgraph package.                                    |
#                                                                              |
#     ssgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  R code for Graphical models based on spike and slab priors                  |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

ssgraph = function( data, n = NULL, method = "ggm", not.cont = NULL, iter = 5000, 
                    burnin = iter / 2, var1 = 4e-04, var2 = 1, lambda = 1, g.prior = 0.5, 
                    g.start = "full", sig.start = NULL, save = FALSE, print = 1000, 
                    cores = NULL )
{
    if( iter < burnin ) stop( "'iter' must be higher than 'burnin'" )
    if(  var1 <= 0    ) stop( "'var1' must be a positive value" )
    if(  var2 <= 0    ) stop( "'var2' must be a positive value" )
    
    burnin <- floor( burnin )
    if( print > iter ) print = iter
    
    cores = BDgraph::get_cores( cores = cores )
    
    list_S_n_p = BDgraph::get_S_n_p( data = data, method = method, n = n, not.cont = not.cont )
    
    S      = list_S_n_p $ S
    n      = list_S_n_p $ n
    p      = list_S_n_p $ p
    method = list_S_n_p $ method
    colnames_data = list_S_n_p $ colnames_data
    
    if( method == "gcgm" )
    {
        not.cont = list_S_n_p $ not.cont
        R        = list_S_n_p $ R
        Z        = list_S_n_p $ Z
        data     = list_S_n_p $ data
        gcgm_NA  = list_S_n_p $ gcgm_NA
    }
    
    g_prior = BDgraph::get_g_prior( g.prior = g.prior, p = p )
    G       = BDgraph::get_g_start( g.start = g.start, g_prior = g_prior, p = p )
    
    if( ( inherits( g.start, "bdgraph" ) ) | ( inherits( g.start, "ssgraph" ) ) ) 
        K <- g.start $ last_K

    if( inherits( g.start, "sim" ) ) 
        K <- g.start $ K
    
    if( ( !inherits( g.start, "sim" ) ) & ( !inherits( g.start, "bdgraph" ) ) & ( !inherits( g.start, "ssgraph" ) ) )
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
    
    p_links = matrix( 0, p, p )
    K_hat   = matrix( 0, p, p )

    cat( paste( c( iter, " MCMC sampling ... in progress: \n" ), collapse = "" ) ) 
    
## - -  main BDMCMC algorithms implemented in C++  - - - - - - - - - - - - - - |
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
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
    
    p_links    = matrix( result $ p_links, p, p, dimnames = list( colnames_data, colnames_data ) ) 
    K_hat      = matrix( result $ K_hat  , p, p, dimnames = list( colnames_data, colnames_data ) ) 
    
    last_graph = matrix( result $ G      , p, p, dimnames = list( colnames_data, colnames_data ) )
    last_K     = matrix( result $ K      , p, p, dimnames = list( colnames_data, colnames_data ) )

    p_links[ lower.tri( p_links, diag = TRUE ) ] = 0
    nmc     = iter - burnin
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
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Summary for the ssgraph object                                            |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

summary.ssgraph = function( object, round = 2, vis = TRUE, ... )
{
    p_links    = object $ p_links
    selected_g = BDgraph::select( p_links, cut = 0.5 )   
    
    if( vis == TRUE )
    {
        if( !is.null( object $ graph_weights ) ) 
            op = graphics::par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) ) 

        # - - - plot selected graph
        sub_g = "Graph with edge posterior probability > 0.5"
        BDgraph::plot.graph( selected_g, main = "Selected graph", sub = sub_g, ... )
        
        if( !is.null( object $ graph_weights ) )
        {
            sample_graphs = object $ sample_graphs
            graph_weights = object $ graph_weights
            sum_gWeights  = sum( graph_weights )

            # - - - plot posterior distribution of graph
            graph_prob = graph_weights / sum_gWeights
            graphics::plot( x = 1 : length( graph_weights ), y = graph_prob, type = "h", col = "gray60", 
                            main = "Posterior probability of graphs",
                            ylab = "Pr( graph | data )", xlab = "graph", ylim = c( 0, max( graph_prob ) ) )

            # - - - plot posterior distribution of graph size
            sizesample_graphs = sapply( sample_graphs, function( x ) length( which( unlist( strsplit( as.character( x ), "" ) ) == 1 ) ) )
            xx       <- unique( sizesample_graphs )
            weightsg <- vector()
            
            for( i in 1 : length( xx ) ) weightsg[ i ] <- sum( graph_weights[ which( sizesample_graphs == xx[ i ] ) ] )
            
            prob_zg = weightsg / sum_gWeights
            graphics::plot( x = xx, y = prob_zg, type = "h", col = "gray10", 
                            main = "Posterior probability of graphs size", 
                            ylab = "Pr( graph size | data )", xlab = "Graph size",
                            ylim = c( 0, max( prob_zg ) ) )
            
            # - - - plot trace of graph size
            all_graphs     = object $ all_graphs
            sizeall_graphs = sizesample_graphs[ all_graphs ]
            
            graphics::plot( x = 1 : length( all_graphs ), sizeall_graphs, type = "l", col = "gray40", 
                            main = "Trace of graph size", ylab = "Graph size", xlab = "Iteration" )

            graphics::par( op )
        }
    }

    return( list( selected_g = selected_g, p_links = round( p_links, round ), K_hat = round( object $ K_hat, round ) ) )
}  
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Plot for the ssgraph object                                               |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

plot.ssgraph = function( x, cut = 0.5, ... )
{
    BDgraph::plot.graph( x, cut = cut, sub = paste0( "Edge posterior probability = ", cut ), ... )    
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Print for the ssgraph object                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
print.ssgraph = function( x, ... )
{
    p_links = x $ p_links
    selected_g = BDgraph::select( p_links, cut = 0.5 )
    
    cat( paste( "\n Adjacency matrix of selected graph \n" ), fill = TRUE )
    print( selected_g )
    
    cat( paste( "\n Edge posterior probabilities of the links \n" ), fill = TRUE )
    print( round( p_links, 2 ) )
} 
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |



  
