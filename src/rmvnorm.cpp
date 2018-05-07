// ------------------------------------------------------------------------------------------------|
//     Copyright (C) 2018 Reza Mohammadi                                                           |
//                                                                                                 |
//     This file is part of ssgraph package.                                                       |
//                                                                                                 |
//     ssgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// ------------------------------------------------------------------------------------------------|

#include "rmvnorm.h"

// ------------------------------------------------------------------------------------------------|
// Simulate one sample from multivarate normal distribution R ~ N_p( mu, sig )
// ------------------------------------------------------------------------------------------------|
void rmvnorm_c( double sample[], double mu[], double sig[], int *p )
{
    GetRNGstate();
    int dim = *p, one = 1;
    char transT = 'T';
    double alpha = 1.0, beta = 1.0;
    
    vector<double> chol_sig( dim * dim );
    cholesky( &sig[0], &chol_sig[0], &dim );        
    
    vector<double> z_N( dim );
    for( int i = 0; i < dim; i++ ) z_N[ i ] = norm_rand();
    
    memcpy( &sample[0], &mu[0], sizeof( double ) * dim );
    F77_NAME(dgemv)( &transT, &dim, &dim, &alpha, &chol_sig[0], &dim, &z_N[0], &one, &beta, &sample[0], &one );
    PutRNGstate();
}
   
