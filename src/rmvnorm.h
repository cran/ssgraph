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

#ifndef rmvnorm_H
#define rmvnorm_H

#include "matrix.h"

extern "C" {
    void rmvnorm_c( double R[], double mu[], double sig[], int *p );
}

#endif
