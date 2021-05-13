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

if("rstudioapi" %in% rownames(installed.packages()) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

if (FALSE)
    rm(list=ls())

source("quSimLib.R")
source("modpow.R")




# ----- AD lecture 18

# Bernstein-Vazirani, (18.4)
n = 8
n = 4
#allbits = expand.grid(rep(list(0:1), n))
for(si in 1:2^n - 1)
{
    # s = unlist(allbits[1 + si,], use.names = F)
    s = as.bits(si, n)
    ket_s = as.qubit(si, n)
    inner = c()
    for(xrow in 1:2^n - 1)
    {
        x = as.bits(xrow, n)
        inner = c(inner, (-1)^sum(s*x))
    }
    inner = t(t(inner)) / sqrt(2^n)
    Hn = multikron(rep(list(H), n))
    stopifnot(max(Mod(Hn %*% inner - ket_s)) < epsilon)
}


# ad Birthday paradox 
# if there were N days in the year, then you’d need about √N 
# people in the room for a likely birthday collision
# N - number of days, n - number of people
N = 365
n = 23
# ≈≈
(1 - 1/N)^(n*(n-1)/2)
.5
# ≈≈
(1-1/N)^(n^2/2)
.5
# ≈≈
(1-1/N)^(n^2)
(.5)^2
# ≈≈
log((1-1/N)^(n^2))
-2*log(2)
# ≈≈
n^2*log(1-1/N)
-2*log(2)
# ≈≈
n^2
-2*log(2)/log(1-1/N)
# ≈≈
n^2
2*log(2)/log(N/(N-1))
# ≈≈
# N --> N+1
n^2
2*log(2)/log((N+1)/N)
# ≈≈
n^2
2*log(2)/log(1+1/N)
# ≈≈
# when z is small: log(1 + z) ≈ z
# (1/N) is small => log(1+1/N) ≈ 1/N
n^2
2*log(2)*N
# ≈≈
n
sqrt(2*log(2)) * sqrt(N)
# ≈≈
n
sqrt(N)


#---------
# https://qiskit.org/textbook/ch-algorithms/bernstein-vazirani.html
#
# bits: o210
CX03l = swap(2, 0, 4) %*% multikron(RCX, I, I) %*% swap(2, 0, 4)
CX13l = swap(1, 2, 4) %*% multikron(RCX, I, I) %*% swap(1, 2, 4)
CX23l = multikron(RCX, I, I)

all(Z %*% H %*% q0 == qm)
bv_func = multikron(I, I, I, I) # 000
bv_func = CX23l %*% CX13l %*% CX03l # 111
bv_func = CX23l %*% CX03l # 101
bv_func = CX23l %*% CX13l # 110
bv_func = CX13l %*% CX03l # 011
{
    full_bv = multikron(I, H, H, H) %*% bv_func %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    # meas bits 210
    res = matrix(0, 8, 1)
    for (val in 0:7)
        res[1+val, 1] = Mod(t(Conj(as.qubit(val, 4))) %*% full_bv)^2 + Mod(t(Conj(as.qubit(8+val, 4))) %*% full_bv)^2
    round(res, 10)
    as.bits(which.max(res) - 1, 3)
}

# ad Simon
if (FALSE)
{
    n=10
    t=1:2^(n/2)
    plot(t, (1-(t*(t-1)/2/(2^n-1)))^t, type='l', col ='blue')
}


