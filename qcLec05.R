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

# ----- AD lecture 5 

# eq. (5.5)
(1 - diag(1, 2))
NOT 

kronecker(NOT, diag(1, 2))
NOTI

q00
q01
q10
q11

NOTI %*% q00   # => q10
NOTI %*% q01   # => q11
NOTI %*% q10   # => q00
NOTI %*% q11   # => q01


# eq. (5.6)
(1 - 2 * (0:1) %*% t(0:1)) / 2^.5
H
matrix(c(rep(1, 3), -1), 2, 2) / 2^.5
H
H %*% 1:0
H %*% 0:1

kronecker(H, H)
H2
H2 %*% q00   # |++〉
H2 %*% q01   # |+-〉
H2 %*% q10   # |-+〉
H2 %*% q11   # |--〉


sum((H2 %*% q00)^2)
sum((H2 %*% q01)^2)
sum((H2 %*% q10)^2)
sum((H2 %*% q11)^2)

round(H2 %*% H2 %*% q00, 10)
round(H2 %*% H2 %*% q01, 10)
round(H2 %*% H2 %*% q10, 10)
round(H2 %*% H2 %*% q11, 10)


# ad order
kronecker(H, diag(1, 2)) * 2^.5
kronecker(diag(1, 2), H) * 2^.5

# eq. (5.7) (5.8) 
kronecker (q0 %*% t(q0), diag(1,2)) + kronecker (q1 %*% t(q1), NOT)
CNOT
NOT
diag(1, 2)
kronecker(H, I)
HI

# (H ⊗ I) (q0 ⊗ q0) == (H q0) ⊗ q0
kronecker(H, I) %*% kronecker(q0, q0) == kronecker(H %*% q0, q0)
HI %*% q00 == kronecker(H %*% q0, q0)

CNOT %*% kronecker(H %*% q0 , q0) # eq. (5.7)
CNOT %*% HI %*% q00

# ok, you can independently enter the sum on the CNOT, and then guess the result
CNOT %*% HI %*% q00
CNOT %*% kronecker(qp, q0)
CNOT %*% (q00 + q10) / sqrt(2)
(CNOT %*% q00 + CNOT %*% q10) / sqrt(2)
(q00 + q11) / sqrt(2)

CNOT %*% HI %*% q00
bell

H %*% q0   #  |+〉
(q00 + q10) / 2^.5
kronecker(H %*% q0, q0)
(q00 + q11) / 2^.5

HI %*% bell

# If Alice sees |0〉
# wg. eq. (5.1)
alfa = .5 # == beta
scale = (alfa^2+alfa^2)^.5
alfa/scale
1/2^.5

# Bell pair, if (after meas) one of pair is |0〉then ...
postmeas(bell, q00, q10) # ... both are |0〉
postmeas(bell, q00, q01) # ... both are |0〉
# Bell pair, if one of pair is |1〉then ...
postmeas(bell, q11, q01) # ... both are |1〉
postmeas(bell, q10, q11) # ... both are |1〉

ifOneIs = q0 # then ...
postmeas(bell, kronecker(ifOneIs, q0), kronecker(ifOneIs, q1))
postmeas(bell, kronecker(q0, ifOneIs), kronecker(q1, ifOneIs))

ifOneIs = q1 # then ...
postmeas(bell, kronecker(q0, ifOneIs), kronecker(q1, ifOneIs))
postmeas(bell, kronecker(ifOneIs, q0), kronecker(ifOneIs, q1))

