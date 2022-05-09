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

#include "util.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

extern "C" {
	void omp_set_num_cores( int *cores ) 
	{
	    #ifdef _OPENMP
	        omp_set_num_threads( *cores );
	    #else
	        Rprintf( "  This OS does not support multi-threading for the ssgraph package  \n" ); 
	    #endif
	}
}
