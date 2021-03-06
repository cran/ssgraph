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

#ifndef matrix_H
#define matrix_H

#include "util.h"

extern "C" {
	void sub_row_mins( double A[], double sub_A[], int *sub, int *p );

	void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p );
	
	void sub_matrix_22( double A[], double A22[], int *sub, int *p );

	void inverse( double A[], double A_inv[], int *p );

	void cholesky( double A[], double U[], int *p );
	
	void update_sigma( double sigma[], int *sub, double K_11_inv[], double K_11_inv_X_K_12[], double *gam, int *p );
	
	void rmvnorm_chol( double sample[], double mean[], double chol_sig[], int *p );
	    
}

#endif
