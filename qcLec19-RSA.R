# ----------------------
# ----------     AD S. Aaronson, Introduction to Quantum Information Science Lecture Notes
# ----------     https://www.scottaaronson.com/qclec.pdf
# ----------------------
# 
# This script is designed to run line by line --> RStudio and Ctrl + Enter
#
#   Author:  Mariusz Krej
#
# ----------------------
# ad RSA encryption
#
# https://www.scottaaronson.com/democritus/lec8.html
# https://www.scottaaronson.com/qclec.pdf
#
#
#   Author:  Mariusz Krej
#

library(pracma)

if("rstudioapi" %in% rownames(installed.packages()) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

if (FALSE)
    rm(list=ls())


source("modpow.R")

epsilon = 1e-12

c = 11
c = 7
c = 5
c = 3

mx = 50
mx = 100 
mx = 1000


while(1)
{
    p = trunc(runif(1, c, mx))
    q = trunc(runif(1, c, mx))
    if (p==q)
        next
    if (!all(abs(q / 2:(q-1) - trunc(q / 2:(q-1))) > epsilon))
        next
    if (!all(abs(p / 2:(p-1) - trunc(p / 2:(p-1))) > epsilon))
        next
    Npq = (p-1) * (q-1)
    if (Npq/c - trunc(Npq/c) != 0)
    {
        x = trunc(runif(1, p*q/2, p*q))
        cat(sprintf("p=%g q=%g x=%g\n", p, q, x))
        break
    }
}    

N = p * q
N   

Npq = (p-1)*(q-1)
Npq
Npq / c

stopifnot(x < N)
stopifnot(Npq / c - trunc(Npq / c) != 0)

if (FALSE)
{
    modpow(x, 0:20, N)
    modpow(x, 0:20 + Npq, N)
    modpow(x, 1, N)
    modpow(x, c, N)
}
retval = modpow(x, c, N)
retval

if(FALSE)
{
    modpow(x, c, N)
    modpow(retval, 1, N)
}

# by raising the retval to the next powers, it jumps over c positions in a row
# but k is such that (c * k) stands next to the end of the period (p-1)(q-1), right where x sits
# (after the end of the period there is 1, because it is a copy of x^0, the next is x^1)
k = which((c*(1:Npq)) %% Npq == 1)[1]
stopifnot(!is.na(k))
if (FALSE)
{
    k
    (c*k) %% Npq
    (c*(1:Npq)) %% Npq
}
# this is x !! 
res = modpow(retval, k, N)
stopifnot(res == x)
cat("OK!\n")


# -------
# + it looks like the period is a few shorter than Npq - that's interesting
# ok, here it is confirmed:
# https://www.scottaaronson.com/qclec.pdf, lecture 19
# "the period of f might not equal φ(N), the most we can say is that the period divides φ(N)"
#
if (FALSE)
{
    which(modpow(retval, 1:Npq, N) == x)
    which(modpow(retval, 1:Npq, N) == 1)
    modpow(retval, k, N)
    plot(1:(Npq),modpow(x, 1:(Npq), N), type = 'l', col = 'blue')
    diff(which(modpow(x, 1:(Npq), N) == 1))

    diff(which(modpow(retval, 1:Npq, N) == x))
    f = diff(which(modpow(retval, 1:Npq, N) == x))[1]
    print(Npq/f)
    modpow(x, 0:k, N)
    modpow(retval, 0:(k+10), N)
}

# ------- period attack (classical)

# (you can pick your own x until it meets "we get lucky" conditions)
if (
    abs(x/p - trunc(x/p)) > epsilon &&
    abs(x/q - trunc(x/q)) > epsilon
   ) # we get lucky
{  
    cat("period search...")
    period = diff(which(modpow(x, 1:N, N) == 1))
    stopifnot(abs(max(period) - min(period)) < epsilon)
    period = period[1]
    cat("ok\n")
        
    rec = NA
    
    if (period %% 2 == 0) # we get lucky (again)
    {
        xs2 = modpow(x, period/2, N)
        if (xs2 > 1 && xs2 < N - 1) # we get lucky (again)
        {
            if (FALSE)
            {
                ((xs2 + 1) * (xs2 - 1)) %% N
                gcd(xs2 - 1, N)
                gcd(xs2 + 1, N)
                c((xs2 - 1) / p, (xs2 - 1) / q, (xs2 + 1) / p, (xs2 + 1) / q)
                pmul = (xs2 - 1) / p
                qmul = (xs2 - 1) / q
                gcd(p*pmul, p*q)
                gcd(q*qmul, p*q)
            }
            cat("we get lucky\n")        
            rec = max(gcd(xs2 - 1, N), gcd(xs2 + 1, N))
        }
    }

    # less lucky
    if (is.na(rec)) 
    {
        # here is a little cheat because Npq should not be known,
        # but maybe this can be circumvented e.g. by checking all pre-matching multipliers?
        mul = Npq / period 
        paq = N - period * mul + 1
        # p^2 - paq*p + N = 0 --> quadratic equation
        d = paq^2 - 4*N  # ---> b^2 - 4*a*c
        stopifnot(d >= 0)
        rec = max( (paq + sqrt(d)) / 2, (paq - sqrt(d)) / 2) # ---> (-b +/- sqrt(d))/(2*a)
    }
    
    cat(sprintf("recovered (p or q) = %d\n", rec))
    stopifnot(any(abs(c(q, p) - rec) < epsilon))

}


