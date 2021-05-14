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


# ----- AD lecture 9 
HI = kronecker(H,I)
ZI = kronecker(Z,I)
XI = kronecker(X,I)

IH = kronecker(I,H)
IZ = kronecker(I,Z)
IX = kronecker(I,X)


# superdense coding

# eq. (9.1), ok, but I think you need to change the order in the last one
XI
ZI
II

# for the 2nd bit there is always I, so we don't touch it
II %*% bell
XI %*% bell
ZI %*% bell
(ZI %*% XI) %*% bell

I %*% q1 

XI %*% ZI %*% bell

# eq. (9.2) --> ???
bobT = matrix(c(1,1,0,0, 0,0,-1,1, 0,0,1,-1, 1,-1,0,0), 4) / sqrt(2)

HI %*% RCX == bobT
HI %*% CX == bobT
IH %*% CX == bobT
IH %*% RCX == bobT

rBobT = HI %*% CNOT
lBobT = SWAP %*% HI %*% CNOT %*% SWAP 

sqrt(2) * bobT
sqrt(2) * lBobT
sqrt(2) * rBobT
sqrt(2) * IH %*% SWAP %*% CNOT %*% SWAP  
sqrt(2) * IH %*% RCX

# rBobT - decode - ok
stopifnot(q00 == round(lBobT %*% bell, 9))
stopifnot(q01 == round(lBobT %*% ZI %*% bell, 9))
stopifnot(q10 == round(lBobT %*% XI %*% bell, 9))
stopifnot(q11 == round(lBobT %*% XI %*% ZI %*% bell, 9))

# rBobT - decode - ok
stopifnot(q00 == round(rBobT %*% bell, 9))
stopifnot(q10 == round(rBobT %*% IZ %*% bell, 9))
stopifnot(q01 == round(rBobT %*% IX %*% bell, 9))
stopifnot(q11 == round(rBobT %*% IX %*% IZ %*% bell, 9))

# bobT - not ok
round(bobT %*% bell, 9)
round(bobT %*% ZI %*% bell, 9)
round(bobT %*% XI %*% bell, 9)
round(bobT %*% ZI %*% XI %*% bell, 9)
round(bobT %*% bell, 9)
round(bobT %*% IZ %*% bell, 9)
round(bobT %*% IX %*% bell, 9)
round(bobT %*% IZ %*% IX %*% bell, 9)

# unitarity test
stopifnot(Mod(lBobT %*% t(Conj(lBobT)) - II) < epsilon)
stopifnot(Mod(lBobT %*% t(Conj(lBobT)) - II) < epsilon)
all(Mod(bobT %*% t(Conj(bobT)) - II) < epsilon) # (fail)


