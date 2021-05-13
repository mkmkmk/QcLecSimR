# Modular exponentiation
# https://en.wikipedia.org/wiki/Modular_exponentiation
#
#
#   Author:  Mariusz Krej
#
library(Rcpp)

if (FALSE)
{
    modpows = function(base, exponent, modulus)
    {
        if (modulus == 1) 
            return(0)
        c = 1
        eprime = 0
        while (eprime < exponent) 
        {
            c = (c * base) %% modulus
            eprime = eprime + 1
        }    
        return(c)
    }
}

if (FALSE)
{
    cfun = 
        "
    // [[Rcpp::export]]
    double modpows(int base, int exponent, long long modulus)
    {
    if (modulus <= 1)
    return 0;
    long long c = 1LL;
    for (int i = 0; i < exponent; ++i)
    c = (c * base) % modulus;
    return c;
    }
    "
    sourceCpp(code = cfun)
} 

if (TRUE)
{
    # repeated squaring power modulo
    cfun = 
        "
            // repeated squaring power modulo
            // [[Rcpp::export]]
            double modpows(int base, int exponent, long long modulus)
            {
                if (modulus <= 1)
                    return 0;
                if (exponent < 1)
                    return 1;
                long long res = 1LL;
                for (int i = 0; i < 32; ++i)
                {
                    if (exponent & 1)
                    {
                        long long ires = base;
                        for (int k = 0; k < i; ++k)
                            ires = (ires * ires) % modulus;
                        res = (res * ires) % modulus;
                    }
                    exponent >>= 1;
                    if (!exponent)
                        break;
                }
                return res;
            }
        "
    sourceCpp(code = cfun)
}

modpow = function(base, exponent, modulus) 
    sapply(exponent, function(it) modpows(base, it, modulus))

# modpow tests
stopifnot(c(modpow(5, 0, 13) == 1))
stopifnot(c(modpow(5, 1, 13) == 5))
stopifnot(c(modpow(777, 1, 13) == 777 %% 13))
stopifnot(c(modpow(5, 3, 13) == 5^3 %% 13))
stopifnot(c(modpow(2, 100, 7) == 2))
stopifnot(c(modpow(3, 100, 7) == 4))
stopifnot(c(modpow(10, 333, 5) == 0))
stopifnot(c(modpow(10, 333, 10) == 0))
stopifnot(c(modpow(10, 333, 2) == 0))
stopifnot(c(modpow(666, 0, 666) == 1))


