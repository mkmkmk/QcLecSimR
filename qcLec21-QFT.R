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

# ----- AD lecture 21

III = multikron(rep(list(I), 3))


# --- ad Continued Fraction
rfun1 = function(a, i) if(i < 100) 1 + 1 / rfun1(a, i+1) else a
rfun1(1, 1) - (1 + sqrt(5)) / 2

cfract = function(r, iter = 50)
{
    i = 1
    tbi = c()
    for (ii in 1:iter)
    {
        i = trunc(r)
        r = r - i
        tbi = c(tbi, i)
        if (abs(r) < epsilon)
            break
        r = 1 / r
    }
    res = "-"
    for(i in rev(tbi))
        res = if (res != "-") sprintf("%g+1/(%s)", i, res) else sprintf("%g", i)
    res
}

cfract(pi, 4)
eval(parse(text = cfract(pi, 4))) - pi
3+16/113 - pi
355/113 - pi

cfract(pi)
abs(eval(parse(text = cfract(pi))) - pi) < 1e-100

# ---------------- 
# https://qiskit.org/textbook/ch-algorithms/shor.html#3.-Qiskit-Implementation

# ---- 7^x mod 15 circuit
n = 4
mul7m15 = multikron(rep(list(X), n)) %*% swap(0, 1, n) %*% swap(1, 2, n) %*% swap(2, 3, n) 
for (power in 1:255)
    stopifnot(all(Reduce("%*%", rep(list(mul7m15), power)) %*% as.qubit(8, n) == as.qubit(sum(2^(1:n-1) * as.bits(modpow(7, power, 15), n)), n)))

# ---- 11^x mod 15 circuit
mul11m15 = multikron(rep(list(X), n)) %*% swap(0, 2, n) %*% swap(1, 3, n)
for (power in 1:255)
    stopifnot(all(Reduce("%*%", rep(list(mul11m15), power)) %*% as.qubit(8, n) == as.qubit(sum(2^(1:n-1) * as.bits(modpow(11, power, 15), n)), n)))


# ---- 4-bit Shor’s Algorithm
#      ┌───┐                                                            ┌──────┐┌─┐         
# q_0: ┤ H ├───────■────────────────────────────────────────────────────┤0     ├┤M├─────────
#      ├───┤       │                                                    │      │└╥┘┌─┐      
# q_1: ┤ H ├───────┼──────────────■─────────────────────────────────────┤1     ├─╫─┤M├──────
#      ├───┤       │              │                                     │  QFT │ ║ └╥┘┌─┐   
# q_2: ┤ H ├───────┼──────────────┼──────────────■──────────────────────┤2     ├─╫──╫─┤M├───
#      ├───┤       │              │              │                      │      │ ║  ║ └╥┘┌─┐
# q_3: ┤ H ├───────┼──────────────┼──────────────┼──────────────■───────┤3     ├─╫──╫──╫─┤M├
#      └───┘┌──────┴──────┐┌──────┴──────┐┌──────┴──────┐┌──────┴──────┐└──────┘ ║  ║  ║ └╥┘
# q_4: ─────┤0            ├┤0            ├┤0            ├┤0            ├─────────╫──╫──╫──╫─
#           │             ││             ││             ││             │         ║  ║  ║  ║ 
# q_5: ─────┤1            ├┤1            ├┤1            ├┤1            ├─────────╫──╫──╫──╫─
#           │  7^1 mod 15 ││  7^2 mod 15 ││  7^4 mod 15 ││  7^8 mod 15 │         ║  ║  ║  ║ 
# q_6: ─────┤2            ├┤2            ├┤2            ├┤2            ├─────────╫──╫──╫──╫─
#      ┌───┐│             ││             ││             ││             │         ║  ║  ║  ║ 
# q_7: ┤ X ├┤3            ├┤3            ├┤3            ├┤3            ├─────────╫──╫──╫──╫─
#      └───┘└─────────────┘└─────────────┘└─────────────┘└─────────────┘         ║  ║  ║  ║ 
# c: 4/══════════════════════════════════════════════════════════════════════════╩══╩══╩══╩═
#                                                                                0  1  2  3 
#
mulmod = mul11m15
mulmod = mul7m15

shor4b =
    multikron(diag(1, 16), t(Conj(QFT(16)))) %*%
    multikron(revctrl(Reduce("%*%", rep(list(mulmod), 8))), diag(1, 8)) %*%
    swap(2, 3, 8) %*% multikron(revctrl(Reduce("%*%", rep(list(mulmod), 4))), diag(1, 8)) %*% swap(2, 3, 8) %*%
    swap(1, 3, 8) %*% multikron(revctrl(Reduce("%*%", rep(list(mulmod), 2))), diag(1, 8)) %*% swap(1, 3, 8) %*%
    swap(0, 3, 8) %*% multikron(revctrl(Reduce("%*%", rep(list(mulmod), 1))), diag(1, 8)) %*% swap(0, 3, 8) %*%
    multikron(multikron(X, III), multikron(rep(list(H), 4))) 

sh4state = shor4b %*% as.qubit(0, 8)

# plot(Mod(sh4state), type='l', col='blue')
which(Mod(sh4state) > 1e-4)

meas4b = c()
for (val in 0:15)
{
    mix = 3
    ires = 0
    for (mix in 0:15)
    {
        bval = c(as.bits(mix, 4), as.bits(val, 4))
        rval = sum(2^(7:0) * bval)
        ires = ires + Mod(t(Conj(as.qubit(rval, 8))) %*% sh4state)^2
    }
    meas4b = c(meas4b, ires)
}

# plot(Mod(meas4b), col='blue')
which(Mod(meas4b) > 1e-4) - 1
diff(which(Mod(meas4b) > 1e-4))

