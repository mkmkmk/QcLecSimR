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

# ------ ad (2.8)
{
    r = runif(1)
    r
    matrix(c(1/2,1/2,1/2,1/2), nrow=2) %*% t(t(c(r,1-r)))
}


# ------ ad fig. 3.2, error
U = matrix(c(1,1,-1,1), nrow=2) / sqrt(2)
U

t(t(c(1,1))) / sqrt(2) # |+〉
qp  

round(U %*% qp, 10)

t(t(1:0))  # |0〉
q0
U %*% q0


t(t(0:1))  # |1〉
q1
U %*% q1

round(U %*% (U %*% q0), 10)

