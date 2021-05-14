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


# ----- AD lecture 7

# --------- ad Bloch sphere
if (FALSE)
{
    plot(exp(1i*pi/2), xlim = c(-1,1), ylim = c(-1,1))
    
    th = 0:99/100*pi
    plot(th, cos(th/2), type = 'l', col = "blue")
    lines(th, sin(th/2), type = 'l', col = "red")
}

th  = c(0, pi, pi/2, pi/2, pi/2, pi/2)
phi = c(0, 0,  0,    pi,   pi/2, -pi/2) 

bloch = toStateTb(th, phi)
bloch
all(Mod(bloch[,1] - q0) < 1e-10)
all(Mod(bloch[,2] - q1) < 1e-10)
all(Mod(bloch[,3] - qp) < 1e-10)
all(Mod(bloch[,4] - qm) < 1e-10)

toBloch(q0)
toBloch(q1)
toBloch(qp)
toBloch(qm)
toBloch((q1 - q0) / sqrt(2))

matrix(c(1, 0, 0, -1), 2, 2)
Z 
matrix(c(0, 1i, -1i, 0), 2, 2)
Y 
X 
NOT
S 
matrix(c(1, 0, 0, 1i), 2, 2)
Tg 
matrix(c(1, 0, 0, exp(1i*pi/4)), 2, 2)

all(Mod(S - (Tg%*%Tg)) < 1e-9)

toBloch(Z %*% qp)
toBloch(qm)
all(Mod(toState(toBloch(Z %*% qp) * pi) - qm) < 1e-9)

toBloch(X %*% q0) == toBloch(q1)

toBloch(S %*% qp)
toBloch(S %*% qm)
S %*% S %*% qm == qp

toState(toBloch(S %*% qm))

# p.7.2 cloning
kronecker(.33*q0+.44*q1, q0)

