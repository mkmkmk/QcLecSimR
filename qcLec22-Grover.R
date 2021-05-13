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


# ----- AD lecture 22
# Grover's alg
# (22.3)
DIFF2 = multikron(H,H) %*% ctrl(Z) %*% multikron(Z,Z) %*% multikron(H,H)
DIFF2 * 2
DIFF3 = -multikron(H,H,H) %*% multikron(X,X,X) %*% ctrl(ctrl(Z)) %*% multikron(X,X,X) %*% multikron(H,H,H)
DIFF3 * 4

all(REV3 %*% DIFF3 %*% REV3 == DIFF3)
all(REV3 %*% ctrl(ctrl(Z)) %*% REV3 == ctrl(ctrl(Z)) )

swap(0,1,3) %*% multikron(I, ctrl(Z)) %*% swap(0,1,3) == multikron(I, ctrl(Z))
swap(1,2,3) %*% multikron(I, ctrl(Z)) %*% swap(1,2,3) == multikron(I, ctrl(Z))


sqrt(8) * multikron(H,H,H) %*% as.qubit(0, 3) 

DIFF3 %*% multikron(ctrl(Z), I) %*% swap(1,2,3) %*% multikron(I, ctrl(Z)) %*% swap(1,2,3) %*% multikron(H,H,H) %*% as.qubit(0, 3) 

multikron(ctrl(Z), I) %*% swap(1,2,3) %*% multikron(I, ctrl(Z)) %*% swap(1,2,3)

round(multikron(H, H) %*% CX %*% multikron(H, H), 10)

as.qubit(0,2) %*% Conj(t(as.qubit(0,2)))


2 * (q00 %*% Conj(t(q00))) - II

# (22.8)
all(multikron(H, H) %*% (2*(q00 %*% Conj(t(q00))) - II) %*% multikron(H, H) - DIFF2 < epsilon)
all(2*multikron(H, H) %*% (q00 %*% Conj(t(q00))) %*% multikron(H, H) - II - DIFF2 < epsilon)
all(2 * (multikron(H, H) %*% q00) %*% Conj(t(multikron(H, H) %*% q00)) - II - DIFF2 < epsilon)


# (22.12) N=8, K=3
(t(Conj(multikron(H, H, H) %*% as.qubit(0, 3))) %*% (as.qubit(2, 3) + as.qubit(4, 3)  + as.qubit(6, 3) ) / sqrt(3) - sqrt(3/8)) < epsilon

t(Conj(multikron(H, H, H) %*% as.qubit(0, 3))) %*% (as.qubit(2, 3) + as.qubit(4, 3)  + as.qubit(6, 3))
matrix(1,1,8) %*% (as.qubit(2, 3) + as.qubit(4, 3)  + as.qubit(6, 3) )

if (FALSE)
{
    per = 100
    t = 0:100
    plot(t, sin(2*pi*t/per)^2, type = 'l', col='blue')
}