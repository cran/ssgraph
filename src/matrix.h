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

#ifndef matrix_H
#define matrix_H

#include "util.h"

extern "C" {
	void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p  );

	void sub_matrix_upper( double A[], double sub_A[], int sub[], int *p_sub, int *p  );

	void sub_row_mins( double A[], double sub_A[], int *sub, int *p );

	void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p );

	void sub_cols_mins( double A[], double sub_A[], int *row, int *col, int *p );

	void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p );
	
	void sub_matrix_22( double A[], double A22[], int *sub, int *p );

	void sub_matrices( double A[], double A11[], double A21[], double A22[], int *row, int *col, int *p );

	void sub_matrices_inv( double A[], double A11_inv[], double A21[], double A22[], int *row, int *col, int *p );
	
	void inverse( double A[], double A_inv[], int *p );

	void cholesky( double A[], double U[], int *p );
}

#endif
