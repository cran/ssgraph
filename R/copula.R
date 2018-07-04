## ------------------------------------------------------------------------------------------------|
# copula funiction based on C++: IT DOES NOT WORKING WITH MESSING DATA
copula_C = function( K, Z, R, is.discrete, S, n, p )
{
    #S = 0 * S
    #Z = 0 * Z
    gcgm = 0
    
    is_discrete = is.discrete
    
    # void copula( double Z[], double K[], int R[], int *n, int *p )
    result = .C( "copula", Z_new = as.double(Z), as.double(K), as.integer(R), as.integer(is_discrete), as.integer(n), as.integer(p), PACKAGE = "ssgraph" )
    Z_new = matrix( result $ Z_new, nrow = n, ncol = p )
    #    S_new = t( Z ) %*% Z
    
    # get_S( double K[], double Z[], int R[], double S[], int *gcgm, int *n, int *p )
    #result = .C( "get_S", as.double(K), as.double(Z), as.integer(R), as.integer(is_discrete), S_new = as.double(S), as.integer(gcgm), as.integer(n), as.integer(p), PACKAGE = "ssgraph" )
    #S_new = matrix( result $ S_new, p, p )
    return( Z_new )
}

## ------------------------------------------------------------------------------------------------|
# copula funiction based on C++: IT DOES NOT WORKING WITH MESSING DATA
copula_S = function( K, Z, R, is.discrete, S, n, p )
{

    #S = 0 * S
    #Z = 0 * Z
    gcgm = 0
    
    is_discrete = is.discrete
    
    # void get_S( double K[], double Z[], int R[], double S[], int *gcgm, int *n, int *p )
    result = .C( "get_S", as.double(K), as.double(Z), as.integer(R), as.integer(is_discrete), S_new = as.double(S), as.integer(gcgm), as.integer(n), as.integer(p), PACKAGE = "ssgraph" )
    S_new = matrix( result $ S_new, p, p )
    return( S_new )
}

