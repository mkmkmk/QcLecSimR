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

# ----- AD lecture 4 

# eq. (4.2) Hadamard
matrix(c(1,1,1,-1), nrow=2) / sqrt(2)
H
t(t(1:0))  # |0〉
q0

t(t(0:1))  # |1〉
q1

H %*% q0
round(H %*% (H %*% q0), 10)

H %*% q1
round(H %*% (H %*% q1), 10)
H %*% (H %*% (H %*% q1))


# --------- ad bomb
N = 10
step1_dud = kronecker(ROT(pi/2/N) %*% q0, q0)
step1_bomb = CNOT %*% step1_dud
step1_bomb
step1_dud
all(step1_dud == CNOT %*% CNOT %*% step1_dud) # $$$

# after step 1
#
# |C〉 - control bit, |P〉 - probe bit
#
#     | before meas          | after meas, if P=|0〉 eq. (5.1)
#     |----------------------|---------------------------------
#     | |CP〉| bomb  | dud   |   | bomb         | dud           
#     |------|-------|-------|---|--------------|-------------
#     | |00〉| 0.988 | 0.988 | * | 0.988/0.988  | 0.988/1     
#     | |01〉| 0.000 | 0.000 |   | 0            | 0         
#     | |10〉| 0.000 | 0.156 | * | 0.000/0.988  | 0.156/1 
#     | |11〉| 0.156 | 0.000 |   | 0            | 0         

round(postmeas(step1_bomb, q00, q10), 3)
round(postmeas(step1_dud, q00, q10), 3)
round(postmeas(step1_bomb, list(q00, q10)), 3)

# bomb "simulation"
bomb = CNOT %*% kronecker(ROT(pi/2/N), diag(1,2))
dud = kronecker(ROT(pi/2/N), diag(1,2))
qcirc = dud
qcirc = bomb
{
    ctl = q0
    probe = q0
    state = kronecker(ctl, probe)
    for(i in 1:N)
    {
        state = qcirc %*% state
        # we assume that the bomb did not explode:
        state = postmeas(state, q00, q10)
    }
    round(state, 3)
}


# -----------------
# -- bomb - swaped bits - ok
#  |αC〉
#  |00〉
#  |01〉
#  |10〉
#  |11〉
# error: swCNOT = kronecker (q0 %*% t(q0), NOT) + kronecker (q1 %*% t(q1),  diag(1,2))
swCNOT = SWAP %*% CNOT %*% SWAP
all(swCNOT == kronecker (diag(1,2), q0 %*% t(q0)) + kronecker (NOT, q1 %*% t(q1)) )
qcirc = swCNOT %*% kronecker(diag(1,2), ROT(pi/2/N))
# --
alfa = q0
meas = list(q00, q01)
meas = list(kronecker(alfa, q0), kronecker(alfa, q1))
# --
alfa = q1
meas = list(q10, q11)
# --
{
    ctl = q0
    state = kronecker(alfa, ctl)
    for(i in 1:N)
        state = postmeas(qcirc %*% state, meas)
    state
    all(state == kronecker(alfa, q0))
}
sum(Mod(state)^2)


