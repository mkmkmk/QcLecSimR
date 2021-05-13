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

XI %*% ZI %*% bell # mistake?


# eq. (9.2) --> mistake?
kronecker(I, H)
IH
errBobT = matrix(c(1,1,0,0, 0,0,-1,1, 0,0,1,-1, 1,-1,0,0), 4) / sqrt(2)
bobT = SWAP %*% HI %*% CNOT %*% SWAP 

bobT == errBobT

sqrt(2) * errBobT
sqrt(2) * bobT
sqrt(2) * IH %*% SWAP %*% CNOT %*% SWAP  
sqrt(2) * IH %*% RCX

#decode - ok
round(bobT %*% bell, 9)
round(bobT %*% ZI %*% bell, 9)
round(bobT %*% XI %*% bell, 9)
round(bobT %*% XI %*% ZI %*% bell, 9)

if (FALSE)
{
    round(errBobT %*% XI %*% bell, 9)
    round(errBobT %*% ZI %*% bell, 9)
    round(errBobT %*% ZI %*% XI %*% bell, 9)
}




