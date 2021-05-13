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



# ----- AD lecture 27

amp = exp(2i * pi * runif(2)) / sqrt(2)
inpBit = amp[1] * q0 + amp[2] * q1

inpState = multikron(q0, q0, inpBit)
bitEncoder = swap(0, 1, 3) %*% multikron(RCX, I) %*% swap(0, 1, 3) %*% multikron(I, RCX)

# Fig. 27.2
bitEncoder %*% inpState 
stopifnot(Mod(bitEncoder %*% inpState - (amp[1] * as.qubit(0, 3) + amp[2] * as.qubit(7, 3))) < epsilon)

bitDecoder = 
    swap(2, 3, 5) %*% multikron(RCX, I, I, I) %*% swap(2, 3, 5) %*%
    swap(1, 3, 5) %*% multikron(RCX, I, I, I) %*% swap(1, 3, 5) %*%
    swap(1, 2, 5) %*% multikron(I, RCX, I, I) %*% swap(1, 2, 5) %*%
    swap(0, 2, 5) %*% multikron(I, RCX, I, I) %*% swap(0, 2, 5) 


# bit flip correction
# Fig. 27.5
for (ii in 1:20)
{
    amp = exp(2i * pi * runif(2)) / sqrt(2)
    inpBit = amp[1] * q0 + amp[2] * q1
    
    randFlip = sample(0:3, 1)
    
    # no error state
    if(randFlip == 0)
        encState = (amp[1] * multikron(q0, q0, q0) + amp[2] * multikron(q1, q1, q1))
    
    # one bit flips
    if (randFlip == 1)
        encState = (amp[1] * multikron(q1, q0, q0) + amp[2] * multikron(q0, q1, q1))
    if (randFlip == 2)
        encState = (amp[1] * multikron(q0, q0, q1) + amp[2] * multikron(q1, q1, q0))
    if (randFlip == 3)
        encState = (amp[1] * multikron(q0, q1, q0) + amp[2] * multikron(q1, q0, q1))
    
    # 2 bit flip
    if (FALSE)
        encState = (amp[1] * multikron(q1, q1, q0) + amp[2] * multikron(q0, q0, q1))
    
    inpState = multikron(q0, q0, encState)
    decState = bitDecoder %*% inpState 
    
    meastb = sapply(0:3, function(m) sum(sapply(0:7, function(k) Mod((t(as.qubit(m*2^3 + k, 5)) %*% decState))^2)))
    meastb
    meas = which.max(meastb) - 1
    meas
    
    stopifnot(abs(meastb - as.qubit(meas, 2)) < epsilon)
    
    mbase = list()
    for(k in 1:2^3)
        mbase[[k]] = multikron(as.qubit(meas, 2), as.qubit(k - 1, 3))
    postm = postmeas(decState, mbase)
    postm
    
    if (meas == 1) # (meas == 2) # (!x && y)
        postm = multikron(I, I, I, I, X) %*% postm
    if (meas == 3) # (x && y)
        postm = multikron(I, I, I, X, I) %*% postm
    if (meas == 2) # (meas == 1) # (x && !y)
        postm = multikron(I, I, X, I, I) %*% postm
    
    refState = multikron(as.qubit(meas, 2), bitEncoder %*% multikron(q0, q0, inpBit))
    stopifnot(Mod(postm - refState) < epsilon)
}
