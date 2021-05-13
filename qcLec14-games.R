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



# ----- AD lecture 14

ifOneIs = q1 # then ...
postmeas(bell, kronecker(q0, ifOneIs), kronecker(q1, ifOneIs))

ifOneIs = q0 # then ...
postmeas(bell, kronecker(q0, ifOneIs), kronecker(q1, ifOneIs))

Mod(Conj(t(q00)) %*% bell)^2
Mod(Conj(t(q11)) %*% bell)^2

all(requal(ROT(pi * 1/8) %*% q1, ROT(pi * 5/8) %*% q0, 9))


# ad eq. (14.3)
Mod(Conj(t(q0)) %*% ROT(pi/8) %*% q0)^2 == cos(pi/8)^2
Mod(Conj(t(q0)) %*% ROT(pi/8) %*% q1)^2 == sin(pi/8)^2
# ad eq. (14.4)
Mod(Conj(t(q1)) %*% ROT(pi/8) %*% q0)^2 == sin(pi/8)^2
Mod(Conj(t(q1)) %*% ROT(pi/8) %*% q1)^2 == cos(pi/8)^2


# x=0 y=0 -> xy==0
Mod(Conj(t(q0)) %*% ROT(pi/8) %*% q0)^2

meas = q0
meas = q1
Mod(Conj(t(q0)) %*% meas)^2
Mod(Conj(t(q0)) %*% ROT(pi/8) %*% meas)^2

round(H2 %*% bell, 9)
round((qpp + qmm) / sqrt(2), 9)
round(H2 %*% bell, 9) == round((qpp + qmm) / sqrt(2), 9)

bell = CNOT %*% HI %*% q00

SWAP %*% kronecker(H %*% q0, q0)
kronecker(H %*% q0, q0)

kronecker(H, ROT(pi/8)) %*% bell 

# a + b = xy

# x=0 y=0  
# if x=0 Alice measures in {|0〉,|1〉}
# if y=0 Bob measures in {|π/8〉,|5π/8〉}
# success when a=0 b=0 or a=1 b=1 ---> Table 13.1
# P(a=0)P(b=0|a=0) + P(a=1)P(b=1|a=1)
# P(a=0) == P. Alice sees |0〉
pa0 = as.numeric(Mod(Conj(t(q00)) %*% bell)^2)
pa0
# P(b=0|a=0) == P. Bob sees |π/8〉when Alice sees |0〉
pb0 = Mod(Conj(t(q0)) %*% ROT(pi/8) %*% q0)^2     # (14.3)
pb0
# P(a=1) == P. Alice sees |1〉
pa1 = as.numeric(Mod(Conj(t(q11)) %*% bell)^2)
pa1
# P(b=1|a=1) == P. Bob sees |5π/8〉when Alice sees |1〉
pb1 = as.numeric(Mod(Conj(t(q1)) %*% ROT(pi * 5/8) %*% q0)^2)      # (14.4)
pb1

pa0*pb0 + pa1*pb1
cos(pi/8)^2

# x=1 y=0  
# success when a=0 b=0 or a=1 b=1

# if x=1 Alice measures in {|+〉,|−〉}
# if y=0 Bob measures in {|π/8〉,|5π/8〉}

# P(a=0) == P. Alice sees |+〉
pa0 = as.numeric(Mod(Conj(t(q0)) %*% qp)^2)
pa0 
# P(b=0|a=0) == P. Bob sees |π/8〉when Alice sees |+〉
pb0 = Mod(Conj(t(qp)) %*% ROT(pi/8) %*% q0)^2    
pb0
# P(a=1) == P. Alice sees |-〉
pa1 = as.numeric(Mod(Conj(t(q1)) %*% qm)^2)
pa1 
# P(b=1|a=1) == P. Bob sees |5π/8〉when Alice sees |+〉
pb1 = as.numeric(Mod(Conj(t(q1)) %*% ROT(pi * 5/8) %*% q0)^2)     
pb1

Mod(Conj(t(ROT(pi/8) %*% q0)) %*% qp)^2    

requal((qpp + qmm) / sqrt(2), (q00 + q11) / sqrt(2), 9)
requal(H2 %*% bell, bell, 9)
requal(H2 %*% H2 %*% H2 %*% H2 %*% H2 %*% bell, bell, 9)
requal(H2, H2 %*% H2 %*% H2 %*% H2 %*% H2, 9)

# ---- H --- M
# ---- R --- M

postmeas(bell, q00, q10)

Mod(Conj(t(q00)) %*% postmeas(IH %*% bell, q00, q10))^2
Mod(Conj(t(q00)) %*% postmeas(HI %*% bell, q00, q10))^2


