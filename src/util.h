#ifndef UTIL_H
#define UTIL_H

#include <Rconfig.h>

#ifdef _OPENMP
    #include <omp.h>
#endif


#define USE_FC_LEN_T  // For Fortran character strings
#include <Rconfig.h>  // included by R.h, so define USE_FC_LEN_T early
#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>

#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <R_ext/Arith.h>     // for the special values like NA, NaN  

#include <R_ext/Visibility.h>

#include <sstream>
#include <string>            // std::string, std::to_string
#include <vector>            // for using vector

#include <math.h>            // isinf, sqrt
#include <limits>            // for std::numeric_limits<double>::max()
#include <algorithm>         // for transform function
#include <functional>        // for transform function

#ifndef FCONE
# define FCONE
#endif

using namespace std;

#endif
