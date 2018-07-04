// ------------------------------------------------------------------------------------------------|
//     Copyright (C) 2018 Reza Mohammadi                                                      |
//                                                                                                 |
//     This file is part of BDgraph package.                                                       |
//                                                                                                 |
//     ssgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// ------------------------------------------------------------------------------------------------|
  
#ifndef copula_H
#define copula_H

#include "matrix.h"

extern "C" {
	void get_mean( double Z[], double K[], double *mu_ij, double *sigma, int *i, int *j, int *n, int *p );

	void get_bounds( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n );

	void copula( double Z[], double K[], int R[], int is_discrete[], int *n, int *p );
	 
	void get_bounds_NA( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n );

	void copula_NA( double Z[], double K[], int R[], int is_discrete[], int *n, int *p );

	void get_S( double K[], double Z[], int R[], int is_discrete[], double S[], int *gcgm, int *n, int *p );
}

#endif
