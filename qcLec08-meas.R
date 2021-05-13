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

#------------------------
# Nimish Mishra, Understanding the basics of measurements in Quantum Computation
# https://towardsdatascience.com/understanding-basics-of-measurements-in-quantum-computation-4c885879eba0
M0 = q0 %*% t(Conj(q0))

M00 = q00 %*% t(Conj(q00))
M01 = q01 %*% t(Conj(q01))
M10 = q10 %*% t(Conj(q10))
M11 = q11 %*% t(Conj(q11))


N = 10
step1_dud = kronecker(ROT(pi/2/N) %*% q0, q0)
step1_bomb = CNOT %*% step1_dud
state = step1_bomb 
state = step1_dud 
meas = M00 + M10
meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]

postmeas2 = function(state, measBasis)
{
    meas = 0
    for (it in measBasis)
        meas = meas + it %*% t(Conj(it))
    meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]
}

round(postmeas2(step1_bomb, list(q00, q10)), 3)
round(postmeas2(step1_dud, list(q00, q10)), 3)

qmm %*% t(Conj(qmm)) + qpm %*% t(Conj(qpm)) + qmp %*% t(Conj(qmp)) + qpp %*% t(Conj(qpp)) 
state = CNOT %*% kronecker(ROT(pi/2/N) %*% q0, qm)
postmeas(state , list(qmm, qpm))
meas = qmm %*% t(Conj(qmm)) +  qpm %*% t(Conj(qpm))
meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]

postmeas(state , list(qmm, qpm))


# ad flip does nothing to |+〉
NOT %*% qp == qp
CNOT %*% qpp == qpp 

# CNOT*CNOT == I, no tak --> crypto: c ⊕ k = m ⊕ k ⊕ k = m
all(CNOT %*% CNOT == diag(1,4)) 

