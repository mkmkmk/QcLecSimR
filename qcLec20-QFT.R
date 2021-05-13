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

# ----- AD lecture 20

if (FALSE)
{
    Q = 2^6
    Q = 4
    plot(exp(2i*pi/Q)^(0:(Q-1)), xlim = c(-1,1), ylim = c(-1,1), col = "blue")
}


B2 = matrix(c(1,0,0,1i), 2)
Mod(diag(exp(2i*pi/4)^(0:1), 2) - B2) < epsilon

B3 = diag(exp(2i*pi/4)^(0:2), 3)
B4 = diag(exp(2i*pi/8)^(0:3), 4)

exp(0)
exp(2i*pi/2) 

QFT4 = QFT(4)
QFT8 = QFT(8) 


# (20.6)
stopifnot(Mod(1/sqrt(2) * rbind(cbind(H, B2 %*% H), cbind(H, -B2 %*% H)) %*% SW - QFT4) < epsilon)

1/sqrt(2) * rbind(cbind(QFT4, B4 %*% QFT4), cbind(QFT4, -B4 %*% QFT4)) %*% REV3

QFT8
preF8 = 2^-.5 * rbind(cbind(QFT4, B4 %*% QFT4), cbind(QFT4, -B4 %*% QFT4))

# abc
# acb
# cab
ROT3 = swap(2, 1, 3) %*% swap(0, 1, 3)
all(Mod(preF8 %*% ROT3 - QFT8) < epsilon)

#abc
#cab
REV3
REV3 - multikron(I, SW) %*% multikron(SW, I) %*% multikron(I, SW)
REV3 - swap(0, 1, 3) %*% swap(2, 1, 3) %*% swap(0, 1, 3)
REV3 - swap(0, 1, 3) %*% ROT3
swap(0, 1, 3) %*% REV3 - swap(0, 1, 3) %*% swap(0, 1, 3) %*% ROT3
swap(0, 1, 3) %*% REV3 - ROT3

ZR4 = diag(0, 4)

# (20.12) - mistake?
all(Mod(multikron(QFT4, I) - rbind(cbind(QFT4, ZR4), cbind(ZR4, QFT4))) < epsilon)
all(Mod(multikron(I, QFT4) - rbind(cbind(QFT4, ZR4), cbind(ZR4, QFT4))) < epsilon)


# (20.14)
# controlled B4
CB4 = ctrl(gate = B4)

stopifnot(Mod(CB4 %*% rbind(cbind(QFT4, ZR4), cbind(ZR4, QFT4)) - rbind(cbind(QFT4, ZR4), cbind(ZR4, B4%*%QFT4))) < epsilon)


# (20.14), (20.12)
stopifnot(Mod(CB4 %*% multikron(I, QFT4) - rbind(cbind(QFT4, ZR4), cbind(ZR4, B4%*%QFT4))) < epsilon)


# (20.15)
multikron(H, I) * sqrt(2)
stopifnot(multikron(H, II) * sqrt(2) == rbind(cbind(II, II), cbind(II, -II)))


# (20.16), (20.15)
stopifnot(Mod(multikron(H, II) %*% rbind(cbind(QFT4, ZR4), cbind(ZR4, B4%*%QFT4)) - preF8) < epsilon)


# (20.16), (20.15), (20.12)
stopifnot(Mod(multikron(H, II) %*% rbind(cbind(QFT4, ZR4), cbind(ZR4, B4%*%QFT4)) - preF8) < epsilon)


# (20.16), (20.12), (20.14)
stopifnot(Mod(multikron(H, II) %*% CB4 %*% multikron(I, QFT4) - preF8) < epsilon)
stopifnot(Mod(multikron(H, II) %*% CB4 %*% multikron(I, QFT4) - preF8) < epsilon)


# ----- QFT4
rQFT4 = SW %*% QFT4 %*% SW
ZR2 = diag(0, 2)
CB3 = ctrl(gate = B2)

stopifnot(Mod(multikron(H, I) %*% CB3 %*% multikron(I, H) %*% SW - QFT4) < epsilon)
stopifnot(Mod(SW %*% multikron(I, H) %*% CB3 %*% multikron(H, I) - QFT4) < epsilon)
stopifnot(Mod(SW %*% (multikron(I, H) %*% CB3 %*% multikron(H, I) %*% SW) %*% SW - QFT4) < epsilon)
stopifnot(Mod(multikron(I, H) %*% CB3 %*% multikron(H, I) %*% SW - rQFT4) < epsilon)
stopifnot(Mod(multikron(I, H) %*% CB3 %*% multikron(H, I) - multikron(I, H) %*% CB3 %*% multikron(H, I)) < epsilon)
stopifnot(Mod(multikron(I, H) %*% CB3 %*% multikron(H, I) %*% SW - rQFT4) < epsilon)


# ----- QFT8
rQFT8 = REV3 %*% QFT8 %*% REV3
stopifnot(Mod(multikron(H, II) %*% CB4 %*% multikron(I, QFT4) %*% ROT3 - QFT8) < epsilon)
stopifnot(Mod(multikron(H, II) %*% CB4 %*% multikron(I, QFT4) %*% swap(0, 1, 3) %*% REV3 - QFT8) < epsilon)

# qclec QFT8 with un-reversed QFT4
stopifnot(Mod(multikron(H, II) %*% CB4 %*% multikron(I, QFT4 %*% SW) %*% REV3 - QFT8) < epsilon)

stopifnot(Mod(REV3 %*% multikron(H, II) %*% CB4 %*% multikron(I, QFT4) %*% swap(0, 1, 3) - rQFT8) < epsilon)
stopifnot(Mod(REV3 %*% multikron(H, II) %*% CB4 %*% multikron(I, QFT4 %*% SW) - rQFT8) < epsilon)

# qiskit QFT8 with un-reversed QFT4
stopifnot(Mod(REV3 %*% multikron(I, SW %*% QFT4) %*% CB4 %*% multikron(H, II) - QFT8) < epsilon)


# ----- QFT16
QFT16 = QFT(16)
B8 = diag(exp(2i*pi/16)^(0:8), 8)
III = multikron(rep(list(I), 3))
# ZR8 = diag(0, 8)
CB16 = ctrl(gate = B8)

stopifnot(Mod(ctrl(gate = B4) - CB4) < epsilon)
stopifnot(Mod(ctrl(gate = B8) - CB16) < epsilon)


# qclec QFT16 with un-reversed QFT8
stopifnot(Mod(multikron(H, III) %*% CB16 %*% multikron(I, QFT8 %*% REV3) %*% REV4 - QFT16) < epsilon)
# qiskit QFT16 with un-reversed QFT8
stopifnot(Mod(REV4 %*% multikron(I, REV3 %*% QFT8) %*% CB16 %*% multikron(H, III) - QFT16) < epsilon)

#-----
all(Mod(QFT4 %*% QFT4 - SW %*% CNOT %*% SW) < epsilon)
all(Mod(rQFT4 %*% rQFT4 - CNOT) < epsilon)
all(Mod(QFT4 %*% t(Conj(QFT4)) - II) < epsilon)

