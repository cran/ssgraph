// ------------------------------------------------------------------------------------------------|
//     Copyright (C) 2018 Reza Mohammadi                                                           |
//                                                                                                 |
//     This file is part of ssgraph package.                                                       |
//                                                                                                 |
//     BDgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// ------------------------------------------------------------------------------------------------|

#include "matrix.h"

// ------------------------------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves square sub_matrix B (p_sub x p_sub), dictated by vector sub
void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int i, j, ixp, subixp, psub = *p_sub, pdim = *p;
	
	for( i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[i] * pdim;
		
		for( j = 0; j < psub; j++ )
			sub_A[ixp + j] = A[subixp + sub[j]]; 
	}
}
   
// ------------------------------------------------------------------------------------------------|
// Takes symmetric matrix A (p x p) and 
// retrieves upper part of sub_matrix B (p_sub x p_sub), dictated by vector sub
void sub_matrix_upper( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int i, j, ixp, subixp, psub = *p_sub, pdim = *p;
			
	for( i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[i] * pdim;
			
		for( j = 0; j <= i; j++ )
			sub_A[ixp + j] = A[subixp + sub[j]]; 
	}
}
   
// ------------------------------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves vector sub_A which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
void sub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int subj = *sub, pdim = *p, subxp = subj * pdim;

	memcpy( sub_A       , A + subxp           , sizeof( double ) * subj );		
	memcpy( sub_A + subj, A + subxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// ------------------------------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(2 x p-2) which is sub rows of matrix A, minus two elements
// Likes A[(i,j), -(i,j)] in R ONLY FOR SYMMETRIC MATRICES
void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col, sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}
}

// ------------------------------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(p-2 x 2) which is sub cols of matrix A, minus two elements
// Likes A[-(i,j), (i,j)] in R 
void sub_cols_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int subi = *row, subj = *col, pdim = *p, p2 = pdim - 2, subixp = subi * pdim, subjxp = subj * pdim;

	memcpy( sub_A           , A + subixp           , sizeof( double ) * subi );		
	memcpy( sub_A + subi    , A + subixp + subi + 1, sizeof( double ) * ( subj - subi - 1 ) );	
	memcpy( sub_A + subj - 1, A + subixp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	

	memcpy( sub_A + p2           , A + subjxp           , sizeof( double ) * subi );		
	memcpy( sub_A + p2 + subi    , A + subjxp + subi + 1, sizeof( double ) * ( subj - subi - 1 ) );	
	memcpy( sub_A + p2 + subj - 1, A + subjxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// ------------------------------------------------------------------------------------------------|
// Takes symmatric matrix A ( p x p ) and 
// retrieves A12 ( 1 x ( p - 1 ) ) and A22 ( ( p - 1 ) x ( p - 1 ) )
// Like A12 = A[ j, -j ], and A22 = A[ -j, -j ] in R
void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, ixpdim, ixp1, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;

	memcpy( A12,        A + subxp,            sizeof( double ) * psub );	
	memcpy( A12 + psub, A + subxp + psub + 1, sizeof( double ) * ( pdim - psub - 1 ) );	

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		ixp1   = i * p1;
		
		memcpy( A22 + ixp1       , A + ixpdim           , sizeof( double ) * psub );
		memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
		ixp1   = ( i - 1 ) * p1;
		
		memcpy( A22 + ixp1       , A + ixpdim           , sizeof( double ) * psub);
		memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}
}
    
// ------------------------------------------------------------------------------------------------|
// Takes symmatric matrix A ( p x p ) and 
// retrieves A22 ( ( p - 1 ) x ( p - 1 ) )
// Like A22 = A[ -j, -j ] in R
void sub_matrix_22( double A[], double A22[], int *sub, int *p )
{
    int i, ixpdim, ixp1, pdim = *p, p1 = pdim - 1, psub = *sub, mpsub = pdim - psub - 1;
    
    for( i = 0; i < psub; i++ )
    {	
        ixpdim = i * pdim;
        ixp1   = i * p1;
        
        memcpy( A22 + ixp1       , A + ixpdim           , sizeof( double ) * psub );
        memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
    }
    
    for( i = psub + 1; i < pdim; i++ )
    {
        ixpdim = i * pdim;
        ixp1   = ( i - 1 ) * p1;
        
        memcpy( A22 + ixp1       , A + ixpdim           , sizeof( double ) * psub);
        memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
    }
}

// ------------------------------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves A11(2x2), A12(2x(p-2)), and A22((p-2)x(p-2))
// Like A11=A[e, e], A12=A[e, -e], and A22=A[-e, -e] in R
void sub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, j, ixp, ij, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[0] = A[sub0 * pdim + sub0];
	A11[1] = A[sub0 * pdim + sub1];
	A11[2] = A11[1];                   // for symmetric matrices
	A11[3] = A[sub1 * pdim + sub1];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[i + i]     = A[ixp + sub0];
		A12[i + i + 1] = A[ixp + sub1];
	
		for( j = 0; j < sub0; j++ )
			A22[j * p2 + i] = A[ixp + j];

		for( j = sub0 + 1; j < sub1; j++ )
		{
			ij = ixp + j;
			A22[(j - 1) * p2 + i] = A[ij];
			A22[i * p2 + j - 1]   = A[ij];
		}
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[(j - 2) * p2 + i] = A[ij];
			A22[i * p2 + j - 2]   = A[ij];
		}
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		
		A12[i + i - 2] = A[ixp + sub0];
		A12[i + i - 1] = A[ixp + sub1];
	
		for( j = sub0 + 1; j < sub1; j++ )
			A22[(j - 1) * p2 + i - 1] = A[ixp + j];
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[(j - 2) * p2 + i - 1] = A[ij];
			A22[(i - 1) * p2 + j - 2] = A[ij];
		}
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
				
		A12[i + i - 4] = A[ixp + sub0];
		A12[i + i - 3] = A[ixp + sub1];
		
		for( j = sub1 + 1; j < pdim; j++ )
			A22[(j - 2) * p2 + i - 2] = A[ixp + j];
	}
}
   
// ------------------------------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves A11_inv(2x2), A21((p-2)x2), and A22((p-2)x(p-2))
// Like A11_inv=inv( A[e, e] ), A21=A[-e, e], and A22=A[-e, -e] in R
void sub_matrices_inv( double A[], double A11_inv[], double A21[], double A22[], int *row, int *col, int *p )
{
	int i, ixp, ixp2, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;
	int sub0xp = sub0 * pdim, sub1xp = sub1 * pdim, sub0_plus = sub0 + 1, sub1_plus = sub1 + 1;
	
	double a11 = A[ sub0 * pdim + sub0 ];
	double a12 = A[ sub0 * pdim + sub1 ];
	double a22 = A[ sub1 * pdim + sub1 ];

	double det_A11 = a11 * a22 - a12 * a12;
	A11_inv[0]     = a22 / det_A11;
	A11_inv[1]     = - a12 / det_A11;
	A11_inv[2]     = A11_inv[1];
	A11_inv[3]     = a11 / det_A11;

	memcpy( A21           , A + sub0xp            , sizeof( double ) * sub0 );		
	memcpy( A21 + sub0    , A + sub0xp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );	
	memcpy( A21 + sub1 - 1, A + sub0xp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );	

	memcpy( A21 + p2           , A + sub1xp            , sizeof( double ) * sub0 );		
	memcpy( A21 + p2 + sub0    , A + sub1xp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );	
	memcpy( A21 + p2 + sub1 - 1, A + sub1xp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );	
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp  = i * pdim;
		ixp2 = i * p2;

		memcpy( A22 + ixp2           , A + ixp            , sizeof( double ) * sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );	
	}
 
	for( i = sub0_plus; i < sub1; i++ )
	{
		ixp  = i * pdim;
		ixp2 = ( i - 1 ) * p2;

		memcpy( A22 + ixp2           , A + ixp            , sizeof( double ) * sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );	
	}
	
	for( i = sub1_plus; i < pdim; i++ )
	{
		ixp  = i * pdim;
		ixp2 = ( i - 2 ) * p2;
				
		memcpy( A22 + ixp2           , A + ixp            , sizeof( double ) * sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );		
	}
}
      
// ------------------------------------------------------------------------------------------------|
// inverse function for symmetric positive-definite matrices (p x p)
// WARNING: Matrix you pass is overwritten with the result
void inverse( double A[], double A_inv[], int *p )
{
	int info, dim = *p;
	char uplo = 'U';

	// creating an identity matrix
	#pragma omp parallel for
	for( int i = 0; i < dim; i++ )
		for( int j = 0; j < dim; j++ )
			A_inv[j * dim + i] = ( i == j );
	
	// LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
	F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info );
}
    
// ------------------------------------------------------------------------------------------------|
// Cholesky decomposition of symmetric positive-definite matrix
// A = U' %*% U
void cholesky( double A[], double U[], int *p )
{
	char uplo = 'U';
	int info, dim = *p, pxp = dim * dim;
	
	memcpy( U, A, sizeof( double ) * pxp );	
	
	F77_NAME(dpotrf)( &uplo, &dim, U, &dim, &info );	
	
	//#pragma omp parallel for
	for( int i = 0; i < dim; i++ )
		for( int j = 0; j < i; j++ )
			U[ j * dim + i ] = 0.0;
}
  
