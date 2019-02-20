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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves vector sub_A which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
    int subj = *sub, pdim = *p, subxp = subj * pdim;
    
    memcpy( sub_A       , A + subxp           , sizeof( double ) * subj );		
    memcpy( sub_A + subj, A + subxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes symmatric matrix A (p x p) and 
// retrieves A12(1x(p-1)) and A22((p-1)x(p-1))
// Like A12=A[j, -j], and A22=A[-j, -j] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
    int pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;
    int size_psub  = sizeof( double ) * psub;
    int size_mpsub = sizeof( double ) * mpsub;
    
    memcpy( A12,        A + subxp,            size_psub );	
    memcpy( A12 + psub, A + subxp + psub + 1, size_mpsub );	

    #pragma omp parallel
    {	
        int i, ixpdim, ixp1;
        
        #pragma omp for
        for( i = 0; i < psub; i++ )
        {	
            ixpdim = i * pdim;
            ixp1   = i * p1;
            
            memcpy( A22 + ixp1       , A + ixpdim           , size_psub );
            memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, size_mpsub );
        }
    }
    
    #pragma omp parallel
    {	
        int i, ixpdim, ixp1;
        
        #pragma omp for
        for( i = psub + 1; i < pdim; i++ )
        {
            ixpdim = i * pdim;
            ixp1   = ( i - 1 ) * p1;
            
            memcpy( A22 + ixp1       , A + ixpdim           , size_psub );
            memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, size_mpsub );
        }
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// inverse function for symmetric positive-definite matrices (p x p)
// WARNING: Matrix you pass is overwritten with the result
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void inverse( double A[], double A_inv[], int *p )
{
    int info, dim = *p;
    char uplo = 'U';
    
    // creating an identity matrix
    #pragma omp parallel for
    for( int i = 0; i < dim; i++ )
        for( int j = 0; j < dim; j++ )
            A_inv[ j * dim + i ] = ( i == j );
    
    // LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
    F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Cholesky decomposition of symmetric positive-definite matrix
// A = U' %*% U
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void cholesky( double A[], double U[], int *p )
{
    char uplo = 'U';
    int info, dim = *p;
    
    memcpy( U, A, sizeof( double ) * dim * dim );	
    
    F77_NAME(dpotrf)( &uplo, &dim, U, &dim, &info );	
    
    #pragma omp parallel for
    for( int i = 0; i < dim; i++ )
        for( int j = 0; j < i; j++ )
            U[ j * dim + i ] = 0.0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// update sigma for gm_spike_slab functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void update_sigma( double sigma[], int *sub, double K_11_inv[], double K_11_inv_X_K_12[], double *gam, int *p )
{
    int i = *sub, dim = *p, p1 = dim - 1; //, one, p1xp1 = p1 * p1;
    
    double alpha_gam = 1.0 / *gam;
    double alpha_ij  = - alpha_gam;

    // sigma[ -i, -i ] = K_11_inv + K_11_inv_X_K_12 %*% t( K_11_inv_X_K_12 ) / gam;
    #pragma omp parallel for
    for( int col = 0; col < i; col++ )
    {
		int col_x_p  = col * dim;
		int col_x_p1 = col * p1;
		
        for( int row = 0; row < i;  row++ ) 
            sigma[ col_x_p + row ] = K_11_inv[ col_x_p1 + row ] + K_11_inv_X_K_12[ col ] * K_11_inv_X_K_12[ row ] * alpha_gam;
        
        for( int row = i; row < p1; row++ ) 
            sigma[ col_x_p + ( row + 1 ) ] = K_11_inv[ col_x_p1 + row ] + K_11_inv_X_K_12[ col ] * K_11_inv_X_K_12[ row ] * alpha_gam;
    }    
    
    #pragma omp parallel for
    for( int col = i; col < p1; col++ )
    {
		int col1_x_p = ( col + 1 ) * dim;
		int col_x_p1 = col * p1;
		
        for( int row = 0; row < i ; row++ ) 
            sigma[ col1_x_p + row ] = K_11_inv[ col_x_p1 + row ] + K_11_inv_X_K_12[ col ] * K_11_inv_X_K_12[ row ] * alpha_gam;
        
        for( int row = i; row < p1; row++ ) 
            sigma[ col1_x_p + ( row + 1 ) ] = K_11_inv[ col_x_p1 + row ] + K_11_inv_X_K_12[ col ] * K_11_inv_X_K_12[ row ] * alpha_gam;
    }
/*
    vector<double> copy_sigma_11( p1xp1 );
    memcpy( &copy_sigma_11[0], &K_11_inv[0], sizeof( double ) * p1xp1 );
    
    // K_11_inv := K_11_inv + K_11_inv_X_K_12 %*% t( K_11_inv_X_K_12 ) / gam;
    F77_NAME(dger)( &p1, &p1, &alpha_gam, &K_11_inv_X_K_12[0], &one, &K_11_inv_X_K_12[0], &one, &copy_sigma_11[0], &p1 );
    
    // sigma[ -i, -i ] = sigma_11
    #pragma omp parallel for
    for( int col = 0; col < i; col++ )
    {
        for( int row = 0; row < i;  row++ ) sigma[ col * dim + row         ] = copy_sigma_11[ col * p1 + row ];
        for( int row = i; row < p1; row++ ) sigma[ col * dim + ( row + 1 ) ] = copy_sigma_11[ col * p1 + row ];
        //memcpy( &sigma[0] + col * dim        , &copy_sigma_11[0] + col * p1    , sizeof( double ) * i          );
        //memcpy( &sigma[0] + col * dim + i + 1, &copy_sigma_11[0] + col * p1 + i, sizeof( double ) * ( p1 - i ) );
    }    
    
    #pragma omp parallel for
    for( int col = i; col < p1; col++ )
    {
        for( int row = 0; row < i ; row++ ) sigma[ ( col + 1 ) * dim + row         ] = copy_sigma_11[ col * p1 + row ];
        for( int row = i; row < p1; row++ ) sigma[ ( col + 1 ) * dim + ( row + 1 ) ] = copy_sigma_11[ col * p1 + row ];
        //memcpy( &sigma[0] + ( col + 1 ) * dim        , &copy_sigma_11[0] + col * p1    , sizeof( double ) * i          );
        //memcpy( &sigma[0] + ( col + 1 ) * dim + i + 1, &copy_sigma_11[0] + col * p1 + i, sizeof( double ) * ( p1 - i ) );
    }
*/
    #pragma omp parallel for
    for( int k = 0; k < i; k++ )
    {
        sigma[ i * dim + k ] = K_11_inv_X_K_12[ k ] * alpha_ij;
        sigma[ k * dim + i ] = K_11_inv_X_K_12[ k ] * alpha_ij;
    }
    
    #pragma omp parallel for
    for( int k = ( i + 1 ); k < dim; k++ )
    {
        int k_p1 = k - 1;
        sigma[ k * dim + i ] = K_11_inv_X_K_12[ k_p1 ] * alpha_ij; 
        sigma[ i * dim + k ] = K_11_inv_X_K_12[ k_p1 ] * alpha_ij; 
    }
    
/*        
    F77_NAME(dscal)( &p1, &alpha_ij, &K_11_inv_X_K_12[0], &one );
    memcpy( &sigma[0] + i * dim        , &K_11_inv_X_K_12[0]    , sizeof( double ) * i          );
    memcpy( &sigma[0] + i * dim + i + 1, &K_11_inv_X_K_12[0] + i, sizeof( double ) * ( p1 - i ) );
    for( int k = 0; k < i; k++ ) sigma[ k * dim + i ]  = K_11_inv_X_K_12[ k ];
    for( int k = ( i + 1 ); k < dim; k++ ) sigma[ k * dim + i ]  = K_11_inv_X_K_12[ k - 1 ]; 
*/    
    sigma[ i * dim + i ]  = alpha_gam;   // sigma[ i, i ]  = 1 / gam;
}    

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Simulate one sample from multivarate normal distribution R ~ N_p( mu, sig )
// where chol_sig = chol( sig )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rmvnorm_chol( double sample[], double mean[], double chol_sig[], int *p )
{
  // GetRNGstate();
   int p1 = *p, one = 1;
   char transN = 'N';
   double alpha = 1.0, beta1 = 1.0;
   
   //    vector<double> chol_sig( p1 * p1 );
   //    cholesky( &sig[0], &chol_sig[0], &p1 );        
   //for( int row = 0; row < p1; row++ )
   //    for( int col = 0; col < row; col++ )
   //       chol_sig[ row * p1 + col ] = 0.0;
   
   vector<double> z_N( p1 );
   for( int row = 0; row < p1; row++ ) z_N[ row ] = norm_rand();
   
   memcpy( &sample[0], &mean[0], sizeof( double ) * p1 );
   F77_NAME(dgemv)( &transN, &p1, &p1, &alpha, &chol_sig[0], &p1, &z_N[0], &one, &beta1, &sample[0], &one );
  // PutRNGstate();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

  
