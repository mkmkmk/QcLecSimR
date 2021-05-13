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



# ----- AD lecture 6

# mixed states eq. (6.2)

# outer product with itself
q0 %*% t(Conj(q0))
q1 %*% t(Conj(q1))

# mixed state eq. (6.3), ok!
(q0 %*% t(Conj(q0)) + q1 %*% t(Conj(q1))) / 2 


bell %*%  t(Conj(bell)) 

# eq. (6.4), ok!
qp = H %*% q0  # |+〉
qm = H %*% q1  # |-〉

(qp %*% t(Conj(qp)) + qm %*% t(Conj(qm))) / 2 

# ad order
qp %*% t(Conj(qp))
t(Conj(qp)) %*% qp  

# ad p. 6.1.2 〈+|M|+〉= 19/2
M = matrix(c(1/2,-10, -10,1/2 ), 2, 2)
t(Conj(qp)) %*% M %*% qp  #〈+|M|+〉
#〈+|M|+〉== (H ρ H^†)[1,1]
HMH = round(H %*% M %*% t(Conj(H)), 10)
HMH
t(Conj(qp)) %*% M %*% qp  #〈+|M|+〉
t(Conj(qm)) %*% M %*% qm  #〈-|M|-〉
t(Conj(qm)) %*% M %*% qp  #〈-|M|+〉
t(Conj(qp)) %*% M %*% qm  #〈+|M|-〉

# ad rank
qr(q1 %*% t(Conj(q1)))$rank
qr((q0 %*% t(Conj(q0)) + q1 %*% t(Conj(q1))) / 2 )$rank
qr(qp %*% t(Conj(qp)))$rank

# eq. (6.9) --> ok but ⊗ is missing
(q00 + q01 + q10) / sqrt(3)
kronecker(q0, qp) * sqrt(2/3) + kronecker(q1, q0) * sqrt(1/3)  # ok!
eq69 = (q00 + q01 + q10) / sqrt(3)

# If Alice sees |0〉Bob sees 
bobA0 = sqrt(2/3) * qp
bobA0
# If Alice sees |1〉Bob sees 
bobA1 = sqrt(1/3) * q0
bobA1

# density mx
bobA0 %*% t(Conj(bobA0)) + bobA1 %*% t(Conj(bobA1))

# ad eq (6.10) tracing out - ok
# (ρB)jj′=∑i αij (αij′)*
aij = matrix(eq69, 2, 2)
rhoB = matrix(0, 2, 2)
for(j in 1+0:1)
    for(jp in 1+0:1)
        rhoB[j, jp] = aij[1, j] * Conj(aij[1, jp]) + aij[2, j] * Conj(aij[2, jp])
rhoB
# αij * αij^†
aij %*% t(Conj(aij))
