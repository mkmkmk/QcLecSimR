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

#if (FALSE)
    rm(list=ls())

source("quSimLib.R")
source("modpow.R")


# ----- AD lecture 28
stopifnot(Mod(X %*% Z - (-1i) * Y) < epsilon)

X %*% qp == qp
-X %*% qm == qm


X %*% q0 == q1
X %*% q1 == q0
Z %*% qp == qm
Z %*% qm == qp
multikron(X, X) %*% bell == bell 


stopifnot(multikron(multikron(rep(list(I), 3)), Z) %*% as.qubit(0, 4) == as.qubit(0, 4))
stopifnot(multikron(Z, multikron(rep(list(I), 3))) %*% as.qubit(0, 4) == as.qubit(0, 4))
stopifnot(multikron(multikron(rep(list(I), 3)), Z) %*% multikron(Z, multikron(rep(list(I), 3)))%*% as.qubit(0, 4) == as.qubit(0, 4))

stopifnot(Mod(multikron(X,X) %*% multikron(X,X) - multikron(I,I)) < epsilon)
stopifnot(Mod(multikron(Z,Z) %*% multikron(Z,Z) - multikron(I,I)) < epsilon)
stopifnot(Mod(multikron(X,X) %*% multikron(Z,Z) + multikron(Y,Y)) < epsilon)
stopifnot(Mod(multikron(Z,Z) %*% multikron(I,I) - multikron(Z,Z)) < epsilon)


# (28.8)
state = multikron(S,I) %*% CX %*% multikron(H, I) %*% q00
stopifnot(multikron(Y,X) %*% state == state)
stopifnot(multikron(Z,Z) %*% state == state)


# Fig. 28.1 --> missing exp(-1i*pi/4) ??
# Magic State 
# https://arxiv.org/pdf/1501.01298.pdf, Fig. 2b
# https://en.wikipedia.org/wiki/Quantum_logic_gate#Rotation_operator_gates
# ROTZ(θ) is none of the gates in Table 4.1
# ROTZ(pi/2) is not a S-gate
# ROTZ(pi/2) is exp(-1i*pi/4)*S 
stopifnot(Mod(ROTZ(pi/2) - exp(-1i*pi/4) * S) < epsilon)
magic = Tg %*% qp
stopifnot(Mod((q0 + exp(1i*pi/4) * q1) / sqrt(2) - magic) < epsilon)
inpBit = t(t(exp(2i*pi*runif(2)) / sqrt(2)))
state = RCX %*% multikron(magic, inpBit)
stopifnot(Mod(multikron(q0, Tg %*% inpBit) - postmeas(state, q00, q01)) < epsilon)
stopifnot(Mod(multikron(q1, Tg %*% inpBit) - multikron(I, ROTZ(pi/2)) %*% postmeas(state, q10, q11)) < epsilon)
# or
stopifnot(Mod(multikron(q0, Tg %*% inpBit) - ctrl(ROTZ(pi/2)) %*% postmeas(state, q00, q01)) < epsilon)
stopifnot(Mod(multikron(q1, Tg %*% inpBit) - ctrl(ROTZ(pi/2)) %*% postmeas(state, q10, q11)) < epsilon)






