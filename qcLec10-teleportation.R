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



# ----- AD lecture 10

# quantum teleportation
# eq. 10.2
alfa = .666
alfa = runif(1)
beta = (1-alfa^2)^.5
eq102 = kronecker(alfa*q0 + beta*q1, bell) 
eq102 * sqrt(2) # ok

# eq. (10.3)
eq103 = kronecker(CNOT, I) %*% eq102 # ok
eq103 * sqrt(2) # ok

# eq. (10.4)
eq104 =  kronecker(HI, I) %*% eq103 
2 * eq104 # ok!

# Table 10.1
# if Alice sees:
ifAliceSees = q00
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
stopifnot(all(Mod(postm - multikron(ifAliceSees, alfa*q0 + beta*q1)) < epsilon))

ifAliceSees = q01
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm = kronecker(II, X) %*% postm
stopifnot(all(Mod(postm - multikron(ifAliceSees, alfa*q0 + beta*q1)) < epsilon))

ifAliceSees = q10
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm = kronecker(II, Z) %*% postm
stopifnot(all(Mod(postm - multikron(ifAliceSees, alfa*q0 + beta*q1)) < epsilon))

ifAliceSees = q11
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm = kronecker(II, Z) %*% kronecker(II, X) %*% postm
stopifnot(all(Mod(postm - multikron(ifAliceSees, alfa*q0 + beta*q1)) < epsilon))
# ok!!!!

CNOTI = kronecker(CNOT, I)
HII = multikron(H, I, I)
for(ii in 1:200)
{
    amp = exp(2i * pi * runif(2)) / sqrt(2)
    inputState = amp[1] * q0 + amp[2] * q1
    transCirc = HII %*% CNOTI %*% kronecker(inputState, bell) 
    # (in fact, Alice's measurement results do not have equal probabilities)
    ifAliceSeesInt = sample(0:3, 1)
    ifAliceSees = as.qubit(ifAliceSeesInt, 2)
    postm = postmeas(transCirc, kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1))
    if (ifAliceSeesInt %% 2 == 1)
        postm = kronecker(II, X) %*% postm
    if (trunc(ifAliceSeesInt / 2) == 1)
        postm = kronecker(II, Z) %*% postm
    stopifnot(all(Mod(postm - multikron(ifAliceSees, inputState)) < 1e-9))
}

# -------
# teleportation of half of another Bell's pair
# bit 0 is an additional bit, bits 1 2 3 are in Fig. 10.1
state = multikron(I, H, I, I) %*% multikron(I, CNOT, I) %*% kronecker(bell, bell) 
ifAliceSees = q00
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
all(postm == multikron(SWAP, SWAP) %*% multikron(q0, bell, q0))
all(multikron(I, SWAP, I) %*% multikron(SWAP, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q01
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, X) %*% postm
postm
all(postm == multikron(SWAP, SWAP) %*% multikron(q0, bell, q1))
all(multikron(I, SWAP, I) %*% multikron(SWAP, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q10
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, Z) %*% postm
postm
all(postm == multikron(SWAP, SWAP) %*% multikron(q1, bell, q0))
all(multikron(I, SWAP, I) %*% multikron(SWAP, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q11
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, Z) %*% multikron(I,I,I, X) %*% postm
postm
all(postm == multikron(SWAP, SWAP) %*% multikron(q1, bell, q1))
all(multikron(I, SWAP, I) %*% multikron(SWAP, I, I) %*% postm == multikron(ifAliceSees, bell))

#   ---x--- a
#      |
# a ---x---
#                  
# b ---x---
#      |
#   ---x--- b 

# a ---x-------
#      |          
#   ---x---x---
#          |
#   -------x--- a
#              
# b ----------- b

# -----------------
# ad GHZ
GHZ

ifISee = q0
postmeas(GHZ, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))

ifISee = q1
postmeas(GHZ, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))

someEn = (multikron(q1, q0, q0) + multikron(q0, q1, q0) + multikron(q0, q0, q1)) / sqrt(3)
ifISee = q0
postmeas(someEn, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))
ifISee = q1
postmeas(someEn, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))


