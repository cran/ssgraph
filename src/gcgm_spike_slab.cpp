// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2018 - 2019 Reza Mohammadi                                                    |
//                                                                                                 |
//     This file is part of ssgraph package.                                                       |
//                                                                                                 |
//     ssgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

#include "matrix.h"
#include "copula.h"

extern "C" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// MCMC sampling algorithm for Gaussian copula graphical models using spike-and-slab priors  
// it is for Bayesian model averaging (MA)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_spike_slab_ma( int *iter, int *burnin, int G[], double K[], double S[], int *p, 
            double K_hat[], double p_links[], int *n,
            double Z[], int R[], int not_continuous[], int *gcgm, 
            double *v1, double *v2, double *lambda, double g_prior[], int *print )
{
    double lambda_c = *lambda;
    double a_gam    = 0.5 * *n + 1.0; 
    
    double sqrt_v1  = sqrt( *v1 );
    double sqrt_v2  = sqrt( *v2 );
    double inv_v1   = 1.0 / *v1;
    double inv_v2   = 1.0 / *v2;

    double alpha = 1.0, beta = 0.0, alpha1 = -1.0;
    double w1, w2, prob_e, inv_V_ij, sigmaii_inv, Sii_lambda, b_gam, gam, K0ii;
    
    int print_c = *print, iteration = *iter, burn_in = *burnin;
    int info, one = 1, dim = *p, pxp = dim * dim;
    int p1 = dim - 1, p1xp1 = p1 * p1, i, row_i, ii, k, ik, k_p1, G_ij;
    
    char transT = 'T', transN = 'N', sideU = 'U', diagN = 'N';																	
    
    vector<double> sigma( pxp ); 
    vector<double> copyK( pxp ); 
    memcpy( &copyK[0], K, sizeof( double ) * pxp );
    inverse( &copyK[0], &sigma[0], &dim );			
    
    vector<double> sigma_11( p1xp1 );       // sigma_11 = sigma[ -i, -i ] 
    vector<double> sigma_12( p1 );          // sigma_12 = sigma[ -i, i  ]
    vector<double> K_11_inv( p1xp1 );
    vector<double> K_u( p1xp1 );            // K_u  = Sii_lambda * K_11_inv
    vector<double> chol_K_u( p1xp1 );       // chol_K_u = chol( K_u )
    vector<double> inv_chol_K_u( p1xp1 );   // inv_chol_K_u = solve( chol_K_u )
    vector<double> sig_K_12( p1xp1, 0.0 );  // sig_K_12 = inv_chol_K_u %*% t( inv_chol_K_u )
    vector<double> mean( p1 );              // mean = -sig_K_12 %*% S[ -i, i ]
    vector<double> S_12( p1 );
    vector<double> K_12( p1 );      // K_12 = ssgraph::rmvnorm( n = 1, mean = mean, sigma = sig_K_12 )
    // vector<double> K_11_inv_X_K_12( p1 );   // K_11_inv_X_K_12 = K_11_inv %*% K_12;
    
    // -- Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - | 
    GetRNGstate();
    for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
    {
        if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

        // - - - STEP 1: copula - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		
        
        get_S( K, Z, R, not_continuous, S, gcgm, n, &dim );
        
        // - - - STEP 2: updating graph G and precision matrix K, row by row - - - - - - - - - - - |		
       
        for( row_i = 0; row_i < dim; row_i++ )
        {
            ii = row_i * dim + row_i;
            
            // - - updating precision matrix - - - - - - - - - - - - - - - - - - - - - - - - - - - |
            sub_matrices1( &sigma[0], &sigma_12[0], &sigma_11[0], &row_i, &dim );
            
            // K_11_inv = sigma_11 - sigma_12 %*% t( sigma_12 ) / sigma[ i, i ]
            sigmaii_inv = - 1.0 / sigma[ ii ];
            memcpy( &K_11_inv[0], &sigma_11[0], sizeof( double ) * p1xp1 );
            F77_NAME(dger)( &p1, &p1, &sigmaii_inv, &sigma_12[0], &one, &sigma_12[0], &one, &K_11_inv[0], &p1 );
            
            Sii_lambda = S[ ii ] + lambda_c;  // Sii_lambda = S[ i, i ] + lambda
            for( k = 0; k < p1xp1; k++ ) 
                K_u[ k ] = Sii_lambda * K_11_inv[ k ];
            
            for( k = 0; k < row_i; k++ )
            {
                inv_V_ij = ( G[ k * dim + row_i ] == 0 ) ? inv_v1 : inv_v2;
                K_u[ k * p1 + k ] += inv_V_ij; 
            }
            
            for( k = row_i + 1; k < dim; k++ )
            {
                k_p1 = k - 1;
                inv_V_ij = ( G[ k * dim + row_i ] == 0 ) ? inv_v1 : inv_v2;
                K_u[ k_p1 * p1 + k_p1 ] += inv_V_ij; 
            }
            
            // - - sampling K_12 from N( mean, sig_K_12 ) - - - - - - - - - - - - - - - - - - - - -|            
            cholesky( &K_u[0], &chol_K_u[0], &p1 );        
            
            memcpy( &inv_chol_K_u[0], &chol_K_u[0], sizeof( double ) * p1xp1 );
            F77_NAME(dtrtri)( &sideU, &diagN, &p1, &inv_chol_K_u[0], &p1, &info );
            
            F77_NAME(dgemm)( &transN, &transT, &p1, &p1, &p1, &alpha, &inv_chol_K_u[0], &p1, &inv_chol_K_u[0], &p1, &beta, &sig_K_12[0], &p1 );
            
            sub_row_mins( S, &S_12[0], &row_i, &dim );   // S_12 = S[ -i, i ]
            
            F77_NAME(dgemv)( &transN, &p1, &p1, &alpha1, &sig_K_12[0], &p1, &S_12[0], &one, &beta, &mean[0], &one );
            
            rmvnorm_chol( &K_12[0], &mean[0], &inv_chol_K_u[0], &p1 );
            
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |            
            
            for( k = 0; k < row_i; k++ )
            {
                K[ k * dim + row_i ] = K_12[ k ];   // K[ -i, i ] = K_12
                K[ row_i * dim + k ] = K_12[ k ];   // K[ i, -i ] = K_12
            }
            
            for( k = row_i + 1; k < dim; k++ )
            {
                k_p1 = k - 1;
                K[ k * dim + row_i ] = K_12[ k_p1 ];
                K[ row_i * dim + k ] = K_12[ k_p1 ];
            }
            
            b_gam = Sii_lambda * 0.5; // b_gam = ( S[ i, i ] + lambda ) * 0.5
            gam = Rf_rgamma( a_gam, 1.0 / b_gam );  // gam = rgamma( n = 1, shape = a_gam, scale = 1 / b_gam )
            
            vector<double> K_11_inv_X_K_12( p1 );   // K_11_inv_X_K_12 = K_11_inv %*% K_12;
            F77_NAME(dgemv)( &transN, &p1, &p1, &alpha, &K_11_inv[0], &p1, &K_12[0], &one, &beta, &K_11_inv_X_K_12[0], &one );
            
            K0ii = F77_NAME(ddot)( &p1, &K_12[0], &one, &K_11_inv_X_K_12[0], &one );
            K[ ii ] = gam + K0ii;   // K[ i, i ] = gam + t( K_12 ) %*% K_11_inv_X_K_12
            
            // - - updating graph matrix - - - - - - - - - - - - - - - - - - - - - - - - - - - - --|
            
            for( k = 0; k < row_i; k++ )
            {
                ik = k * dim + row_i;
                w1 = ( 1.0 - g_prior[ ik ] ) * exp( - 0.5 * K_12[ k ] * K_12[ k ] * inv_v1 ) / sqrt_v1;
                w2 = g_prior[ ik ]           * exp( - 0.5 * K_12[ k ] * K_12[ k ] * inv_v2 ) / sqrt_v2;
                
                prob_e = w2 / ( w1 + w2 );
                G_ij = ( unif_rand() < prob_e ) ? 1 : 0;
                
                G[ row_i * dim + k ] = G_ij;   // G[ k, i ] = G_ij;
                G[ k * dim + row_i ] = G_ij;   // G[ i, k ] = G_ij;
            }
            
            for( k = row_i + 1; k < dim; k++ )
            {
                k_p1 = k - 1;
                ik = k * dim + row_i;
                w1 = ( 1.0 - g_prior[ ik ] ) * exp( - 0.5 * K_12[ k_p1 ] * K_12[ k_p1 ] * inv_v1 ) / sqrt_v1;
                w2 = g_prior[ ik ]           * exp( - 0.5 * K_12[ k_p1 ] * K_12[ k_p1 ] * inv_v2 ) / sqrt_v2;
                
                prob_e = w2 / ( w1 + w2 );
                G_ij = ( unif_rand() < prob_e ) ? 1 : 0;
                
                G[ row_i * dim + k ] = G_ij;   // G[ k, i ] = G_ij;
                G[ k * dim + row_i ] = G_ij;   // G[ i, k ] = G_ij;
            }
            
            // - - updating Covariance matrix according to one-column change of precision matrix --|
            
            update_sigma( &sigma[0], &row_i, &K_11_inv[0], &K_11_inv_X_K_12[0], &gam, &dim );
            
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |            
        }
        
        // - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
        if( i_mcmc >= burn_in )
            for( i = 0; i < pxp ; i++ )
            {
                K_hat[   i ] += K[ i ];
                p_links[ i ] += G[ i ];
            }	
        // - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
    }  
    PutRNGstate();
    // -- End of main loop for MCMC algorithm - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// MCMC sampling algorithm for Gaussian copula graphical models using spike-and-slab priors  
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_spike_slab_map( int *iter, int *burnin, int G[], double K[], double S[], int *p, 
                         double K_hat[], double p_links[], int *n,
                         int all_graphs[], double all_weights[], 
                         char *sample_graphs[], double graph_weights[], int *size_sample_g,
                         double Z[], int R[], int not_continuous[], int *gcgm, 
                         double *v1, double *v2, double *lambda, double g_prior[], int *print )
{
    double lambda_c = *lambda;
    double a_gam    = 0.5 * *n + 1.0; 
    
    double sqrt_v1  = sqrt( *v1 );
    double sqrt_v2  = sqrt( *v2 );
    double inv_v1   = 1.0 / *v1;
    double inv_v2   = 1.0 / *v2;
    
    double alpha = 1.0, beta = 0.0, alpha1 = -1.0;
    double w1, w2, prob_e, inv_V_ij, sigmaii_inv, Sii_lambda, b_gam, gam, K0ii;
    
    int print_c = *print, iteration = *iter, burn_in = *burnin;
    int info, one = 1, dim = *p, pxp = dim * dim;
    int p1 = dim - 1, p1xp1 = p1 * p1, i, j, row_i, ii, k, ik, k_p1, G_ij;
    
    char transT = 'T', transN = 'N', sideU = 'U', diagN = 'N';																	
    
    bool this_one;
    int count_all_g = 0, qp = dim * ( dim - 1 ) / 2;
    int size_sample_graph = *size_sample_g, counter;
    vector<char> char_g( qp );              // char string_g[pp];
    string string_g;
    vector<string> sample_graphs_C( iteration - burn_in );
    
    vector<double> sigma( pxp ); 
    vector<double> copyK( pxp ); 
    memcpy( &copyK[0], K, sizeof( double ) * pxp );
    inverse( &copyK[0], &sigma[0], &dim );			
    
    vector<double> sigma_11( p1xp1 );       // sigma_11 = sigma[ -i, -i ] 
    vector<double> sigma_12( p1 );          // sigma_12 = sigma[ -i, i  ]
    vector<double> K_11_inv( p1xp1 );
    vector<double> K_u( p1xp1 );            // K_u  = Sii_lambda * K_11_inv
    vector<double> chol_K_u( p1xp1 );       // chol_K_u = chol( K_u )
    vector<double> inv_chol_K_u( p1xp1 );   // inv_chol_K_u = solve( chol_K_u )
    vector<double> sig_K_12( p1xp1, 0.0 );  // sig_K_12 = inv_chol_K_u %*% t( inv_chol_K_u )
    vector<double> mean( p1 );              // mean = -sig_K_12 %*% S[ -i, i ]
    vector<double> S_12( p1 );
    vector<double> K_12( p1 );      // K_12 = ssgraph::rmvnorm( n = 1, mean = mean, sigma = sig_K_12 )
    // vector<double> K_11_inv_X_K_12( p1 );   // K_11_inv_X_K_12 = K_11_inv %*% K_12;
    
    //-- Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
    GetRNGstate();
    for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
    {
        if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
        
        //- - - STEP 1: copula - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |		
        
        get_S( K, Z, R, not_continuous, S, gcgm, n, &dim );
        
        //- - - STEP 2: updating graph G and precision matrix K, row by row - - - - - - - - - - ---|		
             
        for( row_i = 0; row_i < dim; row_i++ )
        {
            ii = row_i * dim + row_i;
            
            // --- updating precision matrix - - - - - - - - - - - - - - - - - - - - - - - - - - --|
            sub_matrices1( &sigma[0], &sigma_12[0], &sigma_11[0], &row_i, &dim );
            
            // K_11_inv = sigma_11 - sigma_12 %*% t( sigma_12 ) / sigma[ i, i ]
            sigmaii_inv = - 1.0 / sigma[ ii ];
            memcpy( &K_11_inv[0], &sigma_11[0], sizeof( double ) * p1xp1 );
            F77_NAME(dger)( &p1, &p1, &sigmaii_inv, &sigma_12[0], &one, &sigma_12[0], &one, &K_11_inv[0], &p1 );
            
            Sii_lambda = S[ ii ] + lambda_c;  // Sii_lambda = S[ i, i ] + lambda
            for( k = 0; k < p1xp1; k++ ) 
                K_u[ k ] = Sii_lambda * K_11_inv[ k ];
            
            for( k = 0; k < row_i; k++ )
            {
                inv_V_ij = ( G[ k * dim + row_i ] == 0 ) ? inv_v1 : inv_v2;
                K_u[ k * p1 + k ] += inv_V_ij; 
            }
            
            for( k = row_i + 1; k < dim; k++ )
            {
                k_p1 = k - 1;
                inv_V_ij = ( G[ k * dim + row_i ] == 0 ) ? inv_v1 : inv_v2;
                K_u[ k_p1 * p1 + k_p1 ] += inv_V_ij; 
            }
            
            // --- sampling K_12 from N( mean, sig_K_12 ) - - - - - - - - - - - - - - - - - - - - -|            
            cholesky( &K_u[0], &chol_K_u[0], &p1 );        
            
            memcpy( &inv_chol_K_u[0], &chol_K_u[0], sizeof( double ) * p1xp1 );
            F77_NAME(dtrtri)( &sideU, &diagN, &p1, &inv_chol_K_u[0], &p1, &info );
            
            F77_NAME(dgemm)( &transN, &transT, &p1, &p1, &p1, &alpha, &inv_chol_K_u[0], &p1, &inv_chol_K_u[0], &p1, &beta, &sig_K_12[0], &p1 );
            
            sub_row_mins( S, &S_12[0], &row_i, &dim );   // S_12 = S[ -i, i ]
            
            F77_NAME(dgemv)( &transN, &p1, &p1, &alpha1, &sig_K_12[0], &p1, &S_12[0], &one, &beta, &mean[0], &one );
            
            rmvnorm_chol( &K_12[0], &mean[0], &inv_chol_K_u[0], &p1 );
            
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |            
            
            for( k = 0; k < row_i; k++ )
            {
                K[ k * dim + row_i ] = K_12[ k ];   // K[ -i, i ] = K_12
                K[ row_i * dim + k ] = K_12[ k ];   // K[ i, -i ] = K_12
            }
            
            for( k = row_i + 1; k < dim; k++ )
            {
                k_p1 = k - 1;
                K[ k * dim + row_i ] = K_12[ k_p1 ];
                K[ row_i * dim + k ] = K_12[ k_p1 ];
            }
            
            b_gam = Sii_lambda * 0.5; // b_gam = ( S[ i, i ] + lambda ) * 0.5
            gam   = Rf_rgamma( a_gam, 1.0 / b_gam );  // gam = rgamma( n = 1, shape = a_gam, scale = 1 / b_gam )
            
            vector<double> K_11_inv_X_K_12( p1 );   // K_11_inv_X_K_12 = K_11_inv %*% K_12;
            F77_NAME(dgemv)( &transN, &p1, &p1, &alpha, &K_11_inv[0], &p1, &K_12[0], &one, &beta, &K_11_inv_X_K_12[0], &one );
            
            K0ii = F77_NAME(ddot)( &p1, &K_12[0], &one, &K_11_inv_X_K_12[0], &one );
            K[ ii ] = gam + K0ii;   // K[ i, i ] = gam + t( K_12 ) %*% K_11_inv_X_K_12
            
            // - - updating graph matrix - - - - - - - - - - - - - - - - - - - - - - - - - - - - --|
            
            for( k = 0; k < row_i; k++ )
            {
                ik = k * dim + row_i;
                w1 = ( 1.0 - g_prior[ ik ] ) * exp( - 0.5 * K_12[ k ] * K_12[ k ] * inv_v1 ) / sqrt_v1;
                w2 = g_prior[ ik ]           * exp( - 0.5 * K_12[ k ] * K_12[ k ] * inv_v2 ) / sqrt_v2;
                
                prob_e = w2 / ( w1 + w2 );
                G_ij = ( unif_rand() < prob_e ) ? 1 : 0;
                
                G[ row_i * dim + k ] = G_ij;   // G[ k, i ] = G_ij;
                G[ k * dim + row_i ] = G_ij;   // G[ i, k ] = G_ij;
            }
            
            for( k = row_i + 1; k < dim; k++ )
            {
                k_p1 = k - 1;
                ik = k * dim + row_i;
                w1 = ( 1.0 - g_prior[ ik ] ) * exp( - 0.5 * K_12[ k_p1 ] * K_12[ k_p1 ] * inv_v1 ) / sqrt_v1;
                w2 = g_prior[ ik ]           * exp( - 0.5 * K_12[ k_p1 ] * K_12[ k_p1 ] * inv_v2 ) / sqrt_v2;
                
                prob_e = w2 / ( w1 + w2 );
                G_ij = ( unif_rand() < prob_e ) ? 1 : 0;
                
                G[ row_i * dim + k ] = G_ij;   // G[ k, i ] = G_ij;
                G[ k * dim + row_i ] = G_ij;   // G[ i, k ] = G_ij;
            }
            
            // - - updating Covariance matrix according to one-column change of precision matrix --|
            
            update_sigma( &sigma[0], &row_i, &K_11_inv[0], &K_11_inv_X_K_12[0], &gam, &dim );
            
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |            
        }
        
        // - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
        if( i_mcmc >= burn_in )
        {
            counter = 0;	
            for( j = 1; j < dim; j++ )
                for( i = 0; i < j; i++ )
                    char_g[ counter++ ] = G[ j * dim + i ] + '0'; 
            
            for( i = 0; i < pxp ; i++ )
            {
                K_hat[   i ] += K[ i ];
                p_links[ i ] += G[ i ];
            }	
            
            string_g = string( char_g.begin(), char_g.end() );	
            
            this_one = false;
            for( i = 0; i < size_sample_graph; i++ )
                if( sample_graphs_C[ i ] == string_g )
                {
                    graph_weights[ i ]++;           // += all_weights[count_all_g];
                    all_graphs[ count_all_g ] = i;
                    this_one = true;
                    break;
                } 
                
                if( !this_one || size_sample_graph == 0 )
                {
                    sample_graphs_C[ size_sample_graph ] = string_g;
                    graph_weights[ size_sample_graph ]   = all_weights[ count_all_g ];
                    all_graphs[ count_all_g ]            = size_sample_graph; 
                    size_sample_graph++;				
                }
                
                count_all_g++; 
        } 
        // - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
    }  
    PutRNGstate();
    // -- End of main loop for MCMC algorithm - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
    
    for( int i = 0; i < ( iteration - burn_in ); i++ ) 
    {
        sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
        sample_graphs[ i ][ qp ] = '\0';
    }
    
    *size_sample_g = size_sample_graph;
    
}

} // End of exturn "C"
