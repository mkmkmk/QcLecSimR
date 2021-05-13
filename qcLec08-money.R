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


# ----- AD lecture 8

# quantum money bomb attack
# ad eq. (8.1)

round(ROT(pi/2), 10)
round(ROT(pi/4), 10)
round(ROT(0), 10)
N=10

kronecker(ROT(pi/2/N) %*% q0, q0)

CNOT %*% kronecker(ROT(pi/2/N) %*% q0, q0)  # ok

#  |Cα〉
#  |00〉
#  |01〉
#  |10〉
#  |11〉
postmeas(CNOT %*% kronecker(ROT(pi/2/N) %*% q0, q0), q00, q10) # |00〉-> ok
postmeas(CNOT %*% kronecker(ROT(pi/2/N) %*% q0, q1), q01, q11) # |01〉-> ok

#      |Cα
qpp #  |++〉
qpm #  |+-〉
qmp #  |-+〉
qmm #  |--〉

postmeas(CNOT %*% kronecker(ROT(pi/2/N) %*% q0, qm), qpm, qmm)

# -----------------
qcirc = CNOT %*% kronecker(ROT(pi/2/N), diag(1,2))
hbase = kronecker(diag(1,2), H)
H
# --
alfa = q0
meas = list(q00, q10)
meas = list(kronecker(q0, alfa), kronecker(q1, alfa))
# --
alfa = q1
meas = list(q01, q11)
meas = list(kronecker(q0, alfa), kronecker(q1, alfa))
# --
# for |+〉and |-〉you rotate the basis before measuring 
# so that the measurement is in the { |+〉, |-〉}
# then |+〉<--> |0〉and |-〉<--> |1〉
qcirc = hbase %*% CNOT %*% kronecker(ROT(pi/2/N), diag(1,2))
# --
#alfa = H %*% qm
alfa = qm
#meas = list(qpm, qmm)
meas = list(q01, q11)
meas = list(kronecker(q0, H %*% alfa), kronecker(q1, H %*% alfa))
round(meas[[2]], 3)
# ok - even this little jumping is mentioned!
# --
#alfa = H %*% qp
alfa = qp
#meas = list(qpp, qmp)
meas = list(q00, q10)
meas = list(kronecker(q0, H %*% alfa), kronecker(q1, H %*% alfa))
# --
# ok, |-〉rotates
-CNOT
cminusNOT = kronecker (q0 %*% t(q0), diag(1,2)) + kronecker (q1 %*% t(q1), -NOT)
qcirc = hbase %*% cminusNOT %*% kronecker(ROT(pi/2/N), diag(1,2))
alfa = qm
meas = list(q01, q11)
{
    ctl = q0
    state = kronecker(ctl, alfa)
    for(i in 1:N)
    {
        state = postmeas(qcirc %*% state, meas)
        # state = meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]
        round(state, 3)
    }
    round(state, 3)
}
sum(Mod(state)^2)

