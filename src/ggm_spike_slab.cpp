// ------------------------------------------------------------------------------------------------|
//     Copyright (C) 2018  Reza Mohammadi                                                          |
//                                                                                                 |
//     This file is part of ssgraph package.                                                       |
//                                                                                                 |
//     ssgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// ------------------------------------------------------------------------------------------------|

#include "matrix.h"
#include "rmvnorm.h"

extern "C" {
// ------------------------------------------------------------------------------------------------|
// MCMC sampling algorithm for Gaussian Graphical models using spike-and-slab priors  
// it is for Bayesian model averaging (MA)
// ------------------------------------------------------------------------------------------------|
void ggm_spike_slab_ma( int *iter, int *burnin, int G[], double K[], double S[], int *p, 
            double K_hat[], double p_links[], int *n,
            double *v1, double *v2, double *lambda, double *g_prior, int *print )
{
    double var_1 = *v1;
    double var_2 = *v2;
    double lambda_c = *lambda;
    
    double a_gam = 0.5 * *n + 1.0; 
    double sqrt_v1     = sqrt( var_1 );
    double sqrt_v2     = sqrt( var_2 );
    double g_prior_c     = *g_prior;
    double one_g_prior_c = 1.0 - g_prior_c;
    
    //omp_set_num_threads( 2 );
	int print_c = *print, iteration = *iter, burn_in = *burnin;
	int one = 1, dim = *p, pxp = dim * dim;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
	char transT = 'T', transN = 'N', sideU = 'U';  // , transN = 'N'																	

	double w1, w2, prob_e;
	int G_ij;
	double inv_V_ij;
	int k_p1;
	int p1    = dim - 1;
	int p1xp1 = p1 * p1;

//-- Main loop for birth-death MCMC ---------------------------------------------------------------| 
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
        // sample graph G and precision matrix K        
		for( int i = 0; i < dim; i++ )
		{
		    int ii = i * dim + i;
		    
		    // --- updating precision matrix ------------------------------------------------------|
		    
		    //K_11_inv = solve( K[ -i, -i ] )
		    vector<double> K_11( p1xp1 );
		    sub_matrix_22( &K[0], &K_11[0], &i, &dim );
		    
		    vector<double> K_11_inv( p1xp1 );
		    inverse( &K_11[0], &K_11_inv[0], &p1 );
		    
		    double Sii_lambda = S[ ii ] + lambda_c;  // Sii_lambda = S[ i, i ] + lambda
	        // K_u  = Sii_lambda * K_11_inv
	        vector<double> K_u( p1xp1 );
	        for( int k = 0; k < p1xp1; k++ ) 
	            K_u[ k ] = Sii_lambda * K_11_inv[ k ];
	        
            for( int k = 0; k < i; k++ )
		    {
		        // inv_V_ij = ifelse( G[ i, k ] == 0, 1.0 / v1, 1.0 / v2 );
		        // K_u[ k, k ] = K_u [ k, k ] + inv_V_ij; 
		        inv_V_ij = ( G[ k * dim + i ] == 0 ) ? ( 1.0 / var_1 ) : ( 1.0 / var_2 );
		        K_u[ k * p1 + k ] += inv_V_ij; 
		    }
	    
		    for( int k = i + 1; k < dim; k++ )
		    {
		        k_p1 = k - 1;
		        inv_V_ij = ( G[ k * dim + i ] == 0 ) ? ( 1.0 / var_1 ) : ( 1.0 / var_2 );
		        K_u[ k_p1 * p1 + k_p1 ] += inv_V_ij; 
		    }
		    
            // --- sampling K_12 from N( mean, sig_K_12 ) -----------------------------------------|            
            //sig_K_12 = solve( K_u )
		    vector<double> sig_K_12( p1xp1 );
		    vector<double> copy_K_u( p1xp1 );
		    memcpy( &copy_K_u[0], &K_u[0], sizeof( double ) * p1xp1 );
		    inverse( &copy_K_u[0], &sig_K_12[0], &p1 );	
		    
            // mean = -sig_K_12 %*% S[ -i, i ]
            vector<double> mean( p1 );
            vector<double> S_12( p1 );
            sub_row_mins( S, &S_12[0], &i, &dim );   // S_12 = S[ -i, i ]
            
            F77_NAME(dgemv)( &transN, &p1, &p1, &alpha1, &sig_K_12[0], &p1, &S_12[0], &one, &beta, &mean[0], &one );

            // K_12 = ssgraph::rmvnorm( n = 1, mean = mean, sigma = sig_K_12 )
            vector<double> K_12( p1 );
            rmvnorm_c( &K_12[0], &mean[0], &sig_K_12[0], &p1 );
            
            // ------------------------------------------------------------------------------------|            
		            
            // K[ -i, i ] = K_12
            // K[ i, -i ] = K_12
            for( int k = 0; k < i; k++ )
            {
                K[ k * dim + i ] = K_12[ k ];
                K[ i * dim + k ] = K_12[ k ];
            }
            
            for( int k = i + 1; k < dim; k++ )
            {
                k_p1 = k - 1;
                K[ k * dim + i ] = K_12[ k_p1 ];
                K[ i * dim + k ] = K_12[ k_p1 ];
            }
                
            double b_gam = Sii_lambda * 0.5; // b_gam = ( S[ i, i ] + lambda ) * 0.5
            // gam   = rgamma( n = 1, shape = a_gam, scale = 1 / b_gam )
            double gam = Rf_rgamma( a_gam, 1.0 / b_gam );  // ?? to check
            
            // K_11_inv_X_K_12 = K_11_inv %*% K_12;
            vector<double> K_11_inv_X_K_12( p1 );
            F77_NAME(dgemv)( &transN, &p1, &p1, &alpha, &K_11_inv[0], &p1, &K_12[0], &one, &beta, &K_11_inv_X_K_12[0], &one );

            // K[ i, i ] = gam + t( K_12 ) %*% K_11_inv_X_K_12
            // double K0ii = F77_NAME(ddot)( &p1, &K_12[0], &one, &K_11_inv_X_K_12[0], &one );
            double K0ii = 0;
            for( int k = 0; k < p1; k++ ) K0ii += ( K_12[ k ] * K_11_inv_X_K_12[ k ] ); // ?? to check
            K[ ii ] = gam + K0ii;

            // --- updating graph matrix ----------------------------------------------------------|
            
            for( int k = 0; k < i; k++ )
            {
                // w1 = one_g_prior * exp( - 0.5 * ( ( ( K_12[ k, 1 ] ) ^ 2 ) / v1 ) ) / sqrt_v1;
                // w2 = g_prior     * exp( - 0.5 * ( ( ( K_12[ k, 1 ] ) ^ 2 ) / v2 ) ) / sqrt_v2;
                w1 = one_g_prior_c * exp( - 0.5 * K_12[ k ] * K_12[ k ] / var_1 ) / sqrt_v1;
                w2 = g_prior_c     * exp( - 0.5 * K_12[ k ] * K_12[ k ] / var_2 ) / sqrt_v2;
                
                prob_e = w2 / ( w1 + w2 );
                G_ij = ( unif_rand() < prob_e ) ? 1 : 0;
                
                G[ i * dim + k ] = G_ij;   // G[ k, i ] = G_ij;
                G[ k * dim + i ] = G_ij;   // G[ i, k ] = G_ij;
            }
            
            for( int k = i + 1; k < dim; k++ )
            {
                k_p1 = k - 1;
                // w1 = one_g.prior * exp( - 0.5 * ( ( ( K_12[ k_p1, 1 ] ) ^ 2 ) / v1 ) ) / sqrt_v1;
                // w2 = g.prior     * exp( - 0.5 * ( ( ( K_12[ k_p1, 1 ] ) ^ 2 ) / v2 ) ) / sqrt_v2;
                w1 = one_g_prior_c * exp( - 0.5 * K_12[ k_p1 ] * K_12[ k_p1 ] / var_1 ) / sqrt_v1;
                w2 = g_prior_c     * exp( - 0.5 * K_12[ k_p1 ] * K_12[ k_p1 ] / var_2 ) / sqrt_v2;
                
                prob_e = w2 / ( w1 + w2 );
                G_ij = ( unif_rand() < prob_e ) ? 1 : 0;
                
                G[ i * dim + k ] = G_ij;   // G[ k, i ] = G_ij;
                G[ k * dim + i ] = G_ij;   // G[ i, k ] = G_ij;
            }
		}
		
        //----- saving result ---------------------------------------------------------------------|	
        if( i_mcmc >= burn_in )
            for( int i = 0; i < pxp ; i++ )
            {
                K_hat[i] += K[i];
                p_links[i] += G[i];
            }	
        //----- End of saving result --------------------------------------------------------------|	
	}  
	PutRNGstate();
//-- End of main loop for MCMC algorithm ----------------------------------------------------------| 
}
       
} // End of exturn "C"
