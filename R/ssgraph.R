## ------------------------------------------------------------------------------------------------|
#     Copyright (C) 2018  Reza Mohammadi
#
#     This file is part of ssgraph package.
#
#     ssgraph is free software: you can redistribute it and/or modify it under 
#     the terms of the GNU General Public License as published by the Free 
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.
#
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>
## ------------------------------------------------------------------------------------------------|
#  R code for Graphcial models based on spike and slab priors
## ------------------------------------------------------------------------------------------------|
ssgraph = function( data, n = NULL, iter = 5000, burnin = iter / 2, 
                    v1 = NULL, v2 = NULL, lambda = 1, g.prior = 0.5, 
                    g.start = "full", sig.start = NULL, print = 100 )
{
    if( iter < burnin ) stop( "Number of iteration must be more than number of burn-in" )
   
    if( class( data ) == "sim" ) data <- data $ data
    
    if( !is.matrix( data ) & !is.data.frame( data ) ) stop( "Data must be a matrix or dataframe" )
    if( is.data.frame( data ) ) data <- data.matrix( data )
    if( any( is.na( data ) ) ) stop( "ssgraph does not deal with missing values" )	

    dimd <- dim( data )
    p    <- dimd[ 2 ]
    if( is.null( n ) ) n <- dimd[1]
    p1 = p - 1
    
    if( isSymmetric( data ) )
    {
        if ( is.null( n ) ) stop( "Please specify the number of observations 'n'" )
        cat( "Input is identified as the covriance matrix. \n" )
        S <- data
    }else{
        S <- t( data ) %*% data
    }
    
    if( class( g.start ) == "bdgraph" ) G <- g.start $ last_graph
    if( class( g.start ) == "sim"     ) G <- as.matrix( g.start $ G )
    if( class( g.start ) == "character" && g.start == "empty"  ) G = matrix( 0, p, p )
    if( class( g.start ) == "character" && g.start == "full"   ) G = matrix( 1, p, p )
    if( is.matrix( g.start ) ) G = g.start
    if( ( sum( G == 0 ) + sum( G == 1 ) ) != ( p ^ 2 ) ) stop( "Element of 'g.start', as a matrix, must have 0 or 1" )

    diag( G ) = 0    
    if( !isSymmetric( G ) ) G = G + t( G )
    
    if( is.null( sig.start ) ) sigma = S else sigma = sig.start
    K = solve( sigma )      # precision or concentration matrix (omega)
    
    if( is.null( v1 ) ) v1 = 0.02 ^ 2
    if( is.null( v2 ) ) v2 = ( 50 ^ 2 ) * v1
    pii   = g.prior

    V = matrix( v1, p, p )
    V[ G == 1 ] = v2
    diag( V ) = 0

    nmc   = iter - burnin
    p_hat = matrix( 0, p, p )
    K_hat = matrix( 0, p, p )
    
    a_gam   = 0.5 * n + 1 
    sqrt_v1 = sqrt( v1 )
    sqrt_v2 = sqrt( v2 )
    one_pii = 1 - pii
    
    ind_all = matrix( 0, nrow = p1, ncol = p ) 
    for( i in 1:p  ) ind_all[ , i ] = ( 1:p )[ -i ]
    
    for( i_mcmc in 1 : iter )    
    {
        if ( i_mcmc %% print == 0 )
        {
            cat( paste( c( " iter = ", i_mcmc ), collapse = "" ), "\r" )
            flush.console()	
        }    
        
        # sample sigma and C = inv( sigma )        
        for( i in 1:p  )
        {
            ind_noi = ind_all[ , i ]
            
            # --- updating precision matrix -------------------------------------------------------|
            V_12     = V[ ind_noi, i ]
            sigma_11 = sigma[ ind_noi, ind_noi ] 
            sigma_12 = sigma[ ind_noi, i       ]
            
            K_11_inv = sigma_11 - sigma_12 %*% t( sigma_12 ) / sigma[ i, i ]
            
            K_u = ( S[ i, i ] + lambda ) * K_11_inv + diag( 1 / V_12 )
            #K_u = ( K_u + t( K_u ) ) / 2   # check ?? 
            
            ## --- sampling beta from N( mean, Ci ) -----------------------------------------------|            
            # Ci = solve( K_u )
            chol_K_u     = chol( K_u )
            inv_chol_K_u = solve( chol_K_u )
            Ci           = inv_chol_K_u %*% t( inv_chol_K_u )
            #chol_Ci <- chol( Ci )
            chol_Ci  <- inv_chol_K_u
            
            mean    <- -Ci %*% S[ ind_noi, i ]
            z_N     <- matrix( rnorm( p1 ), nrow = p1, ncol = 1 )
            beta    <- t( chol_Ci ) %*% z_N + mean
            ## ------------------------------------------------------------------------------------|            
            
            K[ ind_noi, i       ] = beta
            K[ i,       ind_noi ] = beta
            
            b_gam = ( S[ i, i ] + lambda ) * 0.5
            gam   = rgamma( n = 1, shape = a_gam, scale = 1 / b_gam )
            
            K[ i, i ] = gam + t( beta ) %*% K_11_inv %*% beta
            
            # --- updating graph matrix -----------------------------------------------------------|
            w1 = one_pii * exp( - 0.5 * ( ( beta ^ 2 ) / v1 ) ) / sqrt_v1
            w2 = pii     * exp( - 0.5 * ( ( beta ^ 2 ) / v2 ) ) / sqrt_v2
            
            prob_e = w2 / ( w1 + w2 )
            
            G_12 = ( runif( p1 ) < prob_e ) * 1
            
            G[ ind_noi, i       ] = G_12;
            G[ i,       ind_noi ] = G_12;
            
            # --- updating matrix V according to one-column change of graph matrix ----------------|
            v              = G_12 * v2
            v[ G_12 == 0 ] = v1
            
            V[ ind_noi, i       ] = v;        
            V[ i,       ind_noi ] = v;
            
            # --- updating Covariance matrix according to one-column change of precision matrix ---|
            K_11_inv_X_beta = K_11_inv %*% beta
            
            sigma[ ind_noi, ind_noi ] = K_11_inv + K_11_inv_X_beta %*% t( K_11_inv_X_beta ) / gam; 
            sigma_12                  = -K_11_inv_X_beta / gam     
            sigma[ ind_noi, i ]       = sigma_12     
            sigma[ i, ind_noi ]       = t( sigma_12 )
            sigma[ i, i ]             = 1 / gam
            # -------------------------------------------------------------------------------------|
        }
        
        if( i_mcmc > burnin )
        {
            p_hat = p_hat + G
            K_hat = K_hat + K
        }
    }
    
    diag( p_hat ) = nmc
    p_hat[ lower.tri( p_hat ) ] = 0
    output = list( p_links = p_hat / nmc, K_hat = K_hat / nmc )
    
    class( output ) = "ssgraph"
    return( output )
}
   
