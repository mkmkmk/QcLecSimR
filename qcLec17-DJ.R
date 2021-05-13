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



# ----- AD lecture 17

#
#---- (17.3) 

x = q0;  f = 0
kronecker(x, NOT %*% as.qubit(f)) 
x = q0;  f = 1
kronecker(x, NOT %*% as.qubit(f)) 
x = q1;  f = 0
kronecker(x, NOT %*% as.qubit(f)) 
x = q1;  f = 1
kronecker(x, NOT %*% as.qubit(f)) 

# (17.1) -> (17.2)  &  (17.1) -> (17.3)
n = 2
for(ix in 1:2^n - 1)
{
    row = as.bits(ix, 2)
    f = row[1]
    x = as.qubit(row[2])
    ok = all(requal(
        2^-.5 * (kronecker(x, as.qubit((0 + f)%%2)) - kronecker(x, as.qubit((1 + f)%%2))),
        (-1)^f * kronecker(x, qm)
    ))
    stopifnot(ok)
}    
# if the second register is placed in the |+〉state, then nothing happens
n = 2
for(ix in 1:2^n - 1)
{
    row = as.bits(ix, 2)
    f = row[1]
    x = as.qubit(row[2])
    ok = all(requal(
        2^-.5 * (kronecker(x, as.qubit((0 + f)%%2)) + kronecker(x, as.qubit((1 + f)%%2))),
        kronecker(x, qp)
    ))
    stopifnot(ok)
}    

x = q1
x = q0
kronecker(x, qm)
kronecker(x, (q0 - q1) / 2^.5)
kronecker(x, q0)/2^.5 - kronecker(x, q1)/2^.5




kronecker(x, qm)


kronecker(x, NOT %*% H %*% as.qubit(f)) 

x = q0;  f = 1
kronecker(H %*% x, NOT %*% H %*% as.qubit(f)) 
kronecker(H %*% x, NOT %*% H %*% as.qubit(f)) 
x = q1;  f = 1
kronecker(H %*% x, NOT %*% H %*% as.qubit(f)) 

x = q0; f = 0


# Deutch
#---- (17.5) 
f = c(0, 1)
(-1)^(1-f)
(-1)^(1) * (-1)^(-f)
-(-1)^(-f)
-(-1)^f
(-1)^(1-f) == -(-1)^f


# ---- (17.9) Deutsch-Josza
multikron(H %*% q0, H %*% q0, H %*% q0, H %*% q0)
multikron(H %*% q0, H %*% q0, H %*% q0, H %*% q1)
multikron(H %*% q0, H %*% q0, H %*% q1, H %*% q0)
multikron(H %*% q0, H %*% q0, H %*% q1, H %*% q1)
multikron(H %*% q0, H %*% q1, H %*% q0, H %*% q0)
multikron(H %*% q0, H %*% q1, H %*% q0, H %*% q1)
multikron(H %*% q0, H %*% q1, H %*% q1, H %*% q0)
multikron(H %*% q0, H %*% q1, H %*% q1, H %*% q1)
multikron(H %*% q1, H %*% q0, H %*% q0, H %*% q0)
multikron(H %*% q1, H %*% q0, H %*% q0, H %*% q1)
multikron(H %*% q1, H %*% q0, H %*% q1, H %*% q0)
multikron(H %*% q1, H %*% q0, H %*% q1, H %*% q1)
multikron(H %*% q1, H %*% q1, H %*% q0, H %*% q0)
multikron(H %*% q1, H %*% q1, H %*% q0, H %*% q1)
multikron(H %*% q1, H %*% q1, H %*% q1, H %*% q0)
multikron(H %*% q1, H %*% q1, H %*% q1, H %*% q1)

# left and right sides of eq. (17.9)
n = 10
n = 8
n = 6
n = 4
for (ix in 1:2^n - 1)
{
    x = as.bits(ix, n)
    right = c()
    for(iy in 1:2^n - 1)
    {
        y = as.bits(iy, n)
        right = c(right, (-1)^sum(x*y))
    }
    right = t(t(right)) / sqrt(2^n)
    
    left = list()
    for(i in 1:length(x))
        left[[i]] = if (x[i] == 0) H %*% q0 else H %*% q1
    
    left = multikron(left)
    stopifnot(max(Mod(left - right)) < epsilon)
}



# ---------------
# ad https://qiskit.org/textbook/ch-algorithms/deutsch-jozsa.html#3.-Creating-Quantum-balanceds--
# Deutsch-Jozsa
# "One of the ways we can guarantee our circuit is balanced is by performing a CNOT for each qubit..."
#

# oracle bit - leftmost bit (bit0, bit ids: 3210)
CX03l = multikron(I, I, SW) %*% multikron(I, SW, I) %*% multikron(SW %*% CX %*% SW, I, I) %*% multikron(I, SW, I) %*% multikron(I, I, SW) 
CX13l = multikron(I, SW, I) %*% multikron(SW %*% CX %*% SW, I, I) %*% multikron(I, SW, I)
CX23l = multikron(SW %*% CX %*% SW, I, I)
balanced_left = CX23l %*% CX13l %*% CX03l


# oracle bit - rightmost bit (bit3, bit ids: 3210)
CX03r = multikron(SW, I, I) %*% multikron(I, SW, I) %*% multikron(I, I, CX) %*% multikron(I, SW, I) %*% multikron(SW, I, I) 
CX13r = multikron(I, SW, I) %*% multikron(I, I, CX) %*% multikron(I, SW, I)
CX23r = multikron(I, I, CX)
balanced_right = CX23r %*% CX13r %*% CX03r 

stopifnot(balanced_left == REV4 %*% balanced_right %*% REV4)
stopifnot(CX03r == swap(1, 3, 4) %*% multikron(I, I, CX) %*% swap(1, 3, 4))
stopifnot(CX13r == swap(1, 2, 4) %*% multikron(I, I, CX) %*% swap(1, 2, 4))
stopifnot(CX23l == multikron(RCX, I, I))

stopifnot(CX03l == swap(2, 0, 4) %*% multikron(RCX, I, I) %*% swap(2, 0, 4))
stopifnot(CX13l == swap(1, 2, 4) %*% multikron(RCX, I, I) %*% swap(1, 2, 4))
stopifnot(CX23l == multikron(RCX, I, I))

stopifnot(all(multikron(SW %*% CX %*% SW, I, I) == multikron(SW, I, I) %*% multikron(CX, I, I) %*% multikron(SW, I, I) ))

balanced_left %*% as.qubit(0, 4)
balanced_left %*% as.qubit(3, 4)
balanced_left %*% as.qubit(5, 4)
balanced_left %*% as.qubit(6, 4)

balanced_left %*% as.qubit(1, 4)
balanced_left %*% as.qubit(2, 4)
balanced_left %*% as.qubit(4, 4)
balanced_left %*% as.qubit(7, 4)

wrap = multikron(I, X, I, I)
wrap = multikron(I, X, I, I)
wrap = multikron(X, X, X, X)

# left vs right
if (F)
{
    # right
    balanced = balanced_right
    balancedwrp = wrap %*% balanced %*% wrap
    const0 = multikron(I, I, I, I)
    const1 = multikron(X, X, X, I)
    
    full_dj = multikron(H, H, H, I) %*% balanced %*% multikron(H, H, H, H %*% X) %*% as.qubit(0, 4)
    full_dj = multikron(H, H, H, I) %*% balancedwrp %*% multikron(H, H, H, H %*% X) %*% as.qubit(0, 4)
    full_dj = multikron(H, H, H, I) %*% const0 %*% multikron(H, H, H, H %*% X) %*% as.qubit(0, 4)
    full_dj = multikron(H, H, H, I) %*% const1 %*% multikron(H, H, H, H %*% X) %*% as.qubit(0, 4)
    
    # meas 
    # P( |000〉)
    Mod(t(Conj(as.qubit(0, 4))) %*% full_dj)^2 + Mod(t(Conj(as.qubit(1, 4))) %*% full_dj)^2
    # P( |111〉)
    Mod(t(Conj(as.qubit(14, 4))) %*% full_dj)^2 + Mod(t(Conj(as.qubit(15, 4))) %*% full_dj)^2
    
} else
{   
    # left
    balanced = balanced_left
    balancedwrp = wrap %*% balanced %*% wrap
    const0 = multikron(I, I, I, I)
    const1 = multikron(I, X, X, X)
    
    full_dj = multikron(I, H, H, H) %*% balanced %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    full_dj = multikron(I, H, H, H) %*% balancedwrp %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    full_dj = multikron(I, H, H, H) %*% const0 %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    full_dj = multikron(I, H, H, H) %*% const1 %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    
    # meas
    # P( |000〉)
    Mod(t(Conj(as.qubit(0, 4))) %*% full_dj)^2 + Mod(t(Conj(as.qubit(8, 4))) %*% full_dj)^2
    # P( |111〉)
    Mod(t(Conj(as.qubit(7, 4))) %*% full_dj)^2 + Mod(t(Conj(as.qubit(15, 4))) %*% full_dj)^2
    
}

all(kronecker(X, I) %*% q10 == q00)
all(kronecker(I, X) %*% q10 == q11)


