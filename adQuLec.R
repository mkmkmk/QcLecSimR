# ----------------------
# ----------     AD S. Aaronson, Introduction to Quantum Information Science Lecture Notes
# ----------     https://www.scottaaronson.com/qclec.pdf
# ----------------------
# 
# This script is designed to run line by line --> RStudio and Ctrl + Enter
#
#
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

# ----- AD lecture 4 

# eq. (4.2) Hadamard
matrix(c(1,1,1,-1), nrow=2) / sqrt(2)
H
t(t(1:0))  # |0〉
q0

t(t(0:1))  # |1〉
q1

H %*% q0
round(H %*% (H %*% q0), 10)

H %*% q1
round(H %*% (H %*% q1), 10)
H %*% (H %*% (H %*% q1))


# --------- ad bomb
N = 10
step1_dud = kronecker(ROT(pi/2/N) %*% q0, q0)
step1_bomb = CNOT %*% step1_dud
step1_bomb
step1_dud
all(step1_dud == CNOT %*% CNOT %*% step1_dud) # $$$

# after step 1
#
# |C〉 - control bit, |P〉 - probe bit
#
#     | before meas          | after meas, if P=|0〉 eq. (5.1)
#     |----------------------|---------------------------------
#     | |CP〉| bomb  | dud   |   | bomb         | dud           
#     |------|-------|-------|---|--------------|-------------
#     | |00〉| 0.988 | 0.988 | * | 0.988/0.988  | 0.988/1     
#     | |01〉| 0.000 | 0.000 |   | 0            | 0         
#     | |10〉| 0.000 | 0.156 | * | 0.000/0.988  | 0.156/1 
#     | |11〉| 0.156 | 0.000 |   | 0            | 0         

round(postmeas(step1_bomb, q00, q10), 3)
round(postmeas(step1_dud, q00, q10), 3)
round(postmeas(step1_bomb, list(q00, q10)), 3)

# almost a bomb simulation
bomb = CNOT %*% kronecker(ROT(pi/2/N), diag(1,2))
dud = kronecker(ROT(pi/2/N), diag(1,2))
qcirc = dud
qcirc = bomb
{
    ctrl = q0
    probe = q0
    state = kronecker(ctrl, probe)
    for(i in 1:N)
    {
        state = qcirc %*% state
        # we assume that the bomb did not explode:
        state = postmeas(state, q00, q10)
    }
    round(state, 3)
}


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


# ----- SWAP 
#
#     00 01 10 11    
# 00   1  0  0  0 
# 01   0  0  1  0    
# 10   0  1  0  0
# 11   0  0  0  1
matrix(c(1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1), 4)
SWAP %*% q00 == q00
SWAP %*% q10 == q01 
SWAP %*% q01 == q10
SWAP %*% q11 == q11


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

# ----- AD lecture 7

# --------- ad Bloch sphere
if (FALSE)
{
    plot(exp(1i*pi/2), xlim = c(-1,1), ylim = c(-1,1))
    
    th = 0:99/100*pi
    plot(th, cos(th/2), type = 'l', col = "blue")
    lines(th, sin(th/2), type = 'l', col = "red")
}

th  = c(0, pi, pi/2, pi/2, pi/2, pi/2)
phi = c(0, 0,  0,    pi,   pi/2, -pi/2) 

bloch = toStateTb(th, phi)
bloch
all(Mod(bloch[,1] - q0) < 1e-10)
all(Mod(bloch[,2] - q1) < 1e-10)
all(Mod(bloch[,3] - qp) < 1e-10)
all(Mod(bloch[,4] - qm) < 1e-10)

toBloch(q0)
toBloch(q1)
toBloch(qp)
toBloch(qm)
toBloch((q1 - q0) / sqrt(2))

matrix(c(1, 0, 0, -1), 2, 2)
Z 
matrix(c(0, 1i, -1i, 0), 2, 2)
Y 
X 
NOT
S 
matrix(c(1, 0, 0, 1i), 2, 2)
Tg 
matrix(c(1, 0, 0, exp(1i*pi/4)), 2, 2)

all(Mod(S - (Tg%*%Tg)) < 1e-9)

toBloch(Z %*% qp)
toBloch(qm)
all(Mod(toState(toBloch(Z %*% qp) * pi) - qm) < 1e-9)

toBloch(X %*% q0) == toBloch(q1)

toBloch(S %*% qp)
toBloch(S %*% qm)
S %*% S %*% qm == qp

toState(toBloch(S %*% qm))

# p.7.2 cloning
kronecker(.33*q0+.44*q1, q0)

# ----- scratchpad 
CNOT %*% kronecker(qp, q0)
toBloch(X %*% q0) == toBloch(q1)
toState(toBloch(qp)*pi)
CNOT %*% kronecker(H, I) 
CNOT %*% HI %*% q00
kronecker(H, I) %*% kronecker(q0, q0)
CNOT %*% kronecker(H %*% q0, q0)
CNOT %*% kronecker(qp, q0)
(q00 + q11) / sqrt(2)

CNOT %*% q00 == q00
CNOT %*% q01 == q01
CNOT %*% q10 == q11
CNOT %*% q11 == q10

CNOT %*% kronecker(q0, qm) == kronecker(q0, qm)
CNOT %*% kronecker(q0, qp) == kronecker(q0, qp)
CNOT %*% kronecker(q1, qm) == kronecker(q1, qm)  # FALSE
CNOT %*% kronecker(q1, qp) == kronecker(q1, qp)  

n = 4
multikron(rep(list(H), n)) %*% multikron(rep(list(q0), n))

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
    ctrl = q0
    state = kronecker(ctrl, alfa)
    for(i in 1:N)
    {
        state = postmeas(qcirc %*% state, meas)
        # state = meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]
        round(state, 3)
    }
    round(state, 3)
}
sum(Mod(state)^2)


# -----------------
# -- bomb - swaped bits - ok
#  |αC〉
#  |00〉
#  |01〉
#  |10〉
#  |11〉
# error: swCNOT = kronecker (q0 %*% t(q0), NOT) + kronecker (q1 %*% t(q1),  diag(1,2))
swCNOT = SWAP %*% CNOT %*% SWAP
all(swCNOT == kronecker (diag(1,2), q0 %*% t(q0)) + kronecker (NOT, q1 %*% t(q1)) )
qcirc = swCNOT %*% kronecker(diag(1,2), ROT(pi/2/N))
# --
alfa = q0
meas = list(q00, q01)
meas = list(kronecker(alfa, q0), kronecker(alfa, q1))
# --
alfa = q1
meas = list(q10, q11)
# --
{
    ctrl = q0
    state = kronecker(alfa, ctrl)
    for(i in 1:N)
        state = postmeas(qcirc %*% state, meas)
    state
    all(state == kronecker(alfa, q0))
}
sum(Mod(state)^2)

    
#------------------------
# Nimish Mishra, Understanding the basics of measurements in Quantum Computation
# https://towardsdatascience.com/understanding-basics-of-measurements-in-quantum-computation-4c885879eba0
M0 = q0 %*% t(Conj(q0))

M00 = q00 %*% t(Conj(q00))
M01 = q01 %*% t(Conj(q01))
M10 = q10 %*% t(Conj(q10))
M11 = q11 %*% t(Conj(q11))

state = step1_bomb 
state = step1_dud 
meas = M00 + M10
meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]

postmeas2 = function(state, measBasis)
{
    meas = 0
    for (it in measBasis)
        meas = meas + it %*% t(Conj(it))
    meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]
}

round(postmeas2(step1_bomb, list(q00, q10)), 3)
round(postmeas2(step1_dud, list(q00, q10)), 3)

qmm %*% t(Conj(qmm)) + qpm %*% t(Conj(qpm)) + qmp %*% t(Conj(qmp)) + qpp %*% t(Conj(qpp)) 
state = CNOT %*% kronecker(ROT(pi/2/N) %*% q0, qm)
postmeas(state , list(qmm, qpm))
meas = qmm %*% t(Conj(qmm)) +  qpm %*% t(Conj(qpm))
meas %*% state / sqrt(Conj(t(state)) %*% meas %*% state)[1,1]

postmeas(state , list(qmm, qpm))


# ad flip does NOThing to |+〉
NOT %*% qp == qp
CNOT %*% qpp == qpp 

# CNOT*CNOT == I, no tak --> crypto: c ⊕ k = m ⊕ k ⊕ k = m
all(CNOT %*% CNOT == diag(1,4)) 


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

XI %*% ZI %*% bell # error


# eq. (9.2) - error?
kronecker(I, H)
IH
errBobT = matrix(c(1,1,0,0, 0,0,-1,1, 0,0,1,-1, 1,-1,0,0), 4) / sqrt(2)
bobT = SWAP %*% HI %*% CNOT %*% SWAP 

bobT == errBobT

sqrt(2) * errBobT
sqrt(2) * bobT
sqrt(2) * IH %*% SWAP %*% CNOT %*% SWAP  
sqrt(2) * IH %*% swCNOT

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


# ----- AD lecture 10

# quantum teleportation

# eq. 10.2
alfa = .666
beta = (1-alfa^2)^.5
eq102 = kronecker(alfa*q0 + beta*q1, bell) 
eq102 * sqrt(2) # ok

# eq. (10.3)
eq103 = kronecker(CNOT, I) %*% eq102 # ok
eq103 * sqrt(2) # ok

# eq. (10.4)
eq104 =  kronecker(HI, I) %*% eq103 
2 * eq104 # ok!

# Table 10.1
# if Alice sees:
ifAliceSees = q00
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm

ifAliceSees = q01
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm
postm = kronecker(II, X) %*% postm
postm
all(Mod(postm - multikron(ifAliceSees, alfa*q0 + beta*q1)) < 1e-9)


ifAliceSees = q10
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm
kronecker(II, Z) %*% postm

ifAliceSees = q11
postm = postmeas(eq104, list(kronecker(ifAliceSees, q0), kronecker(ifAliceSees, q1)))
postm
kronecker(II, Z) %*% kronecker(II, X) %*% postm
# ok!!!!

# -------
# teleportation of half of another Bell's pair
# bit 0 is an additional bit, bits 1 2 3 are in Fig. 10.1
state = multikron(I, H, I, I) %*% multikron(I, CNOT, I) %*% kronecker(bell, bell) 
ifAliceSees = q00
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
all(postm == multikron(SWAP, SWAP) %*% multikron(q0, bell, q0))
all(multikron(I, SWAP, I) %*% multikron(SWAP, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q01
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, X) %*% postm
postm
all(postm == multikron(SWAP, SWAP) %*% multikron(q0, bell, q1))
all(multikron(I, SWAP, I) %*% multikron(SWAP, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q10
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, Z) %*% postm
postm
all(postm == multikron(SWAP, SWAP) %*% multikron(q1, bell, q0))
all(multikron(I, SWAP, I) %*% multikron(SWAP, I, I) %*% postm == multikron(ifAliceSees, bell))
## --
ifAliceSees = q11
postm = postmeas(state, multikron(q0, ifAliceSees, q0), multikron(q0, ifAliceSees, q1), multikron(q1, ifAliceSees, q0), multikron(q1, ifAliceSees, q1))
postm
postm = multikron(I,I,I, Z) %*% multikron(I,I,I, X) %*% postm
postm
all(postm == multikron(SWAP, SWAP) %*% multikron(q1, bell, q1))
all(multikron(I, SWAP, I) %*% multikron(SWAP, I, I) %*% postm == multikron(ifAliceSees, bell))

#   ---x--- a
#      |
# a ---x---
#                  
# b ---x---
#      |
#   ---x--- b 

# a ---x-------
#      |          
#   ---x---x---
#          |
#   -------x--- a
#              
# b ----------- b

# -----------------
# ad GHZ
GHZ

ifISee = q0
postmeas(GHZ, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))

ifISee = q1
postmeas(GHZ, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))

someEn = (multikron(q1, q0, q0) + multikron(q0, q1, q0) + multikron(q0, q0, q1)) / sqrt(3)
ifISee = q0
postmeas(someEn, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))
ifISee = q1
postmeas(someEn, multikron(q0, q0, ifISee), multikron(q0, q1, ifISee), multikron(q1, q0, ifISee), multikron(q1, q1, ifISee))



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



# ----- AD lecture 17

#
#---- (17.3) 

x = q0;  f = 0
kronecker(x, NOT %*% as.qubit(f)) 
x = q0;  f = 1
kronecker(x, NOT %*% as.qubit(f)) 
x = q1;  f = 0
kronecker(x, NOT %*% as.qubit(f)) 
x = q1;  f = 1
kronecker(x, NOT %*% as.qubit(f)) 

# (17.1) -> (17.2)  &  (17.1) -> (17.3)
n = 2
for(ix in 1:2^n - 1)
{
    row = as.bits(ix, 2)
    f = row[1]
    x = as.qubit(row[2])
    ok = all(requal(
        2^-.5 * (kronecker(x, as.qubit((0 + f)%%2)) - kronecker(x, as.qubit((1 + f)%%2))),
        (-1)^f * kronecker(x, qm)
    ))
    stopifnot(ok)
}    
# if the second register is placed in the |+〉state, then nothing happens
n = 2
for(ix in 1:2^n - 1)
{
    row = as.bits(ix, 2)
    f = row[1]
    x = as.qubit(row[2])
    ok = all(requal(
        2^-.5 * (kronecker(x, as.qubit((0 + f)%%2)) + kronecker(x, as.qubit((1 + f)%%2))),
        kronecker(x, qp)
        ))
    stopifnot(ok)
}    

x = q1
x = q0
kronecker(x, qm)
kronecker(x, (q0 - q1) / 2^.5)
kronecker(x, q0)/2^.5 - kronecker(x, q1)/2^.5




kronecker(x, qm)


kronecker(x, NOT %*% H %*% as.qubit(f)) 

x = q0;  f = 1
kronecker(H %*% x, NOT %*% H %*% as.qubit(f)) 
kronecker(H %*% x, NOT %*% H %*% as.qubit(f)) 
x = q1;  f = 1
kronecker(H %*% x, NOT %*% H %*% as.qubit(f)) 

x = q0; f = 0


# Deutch
#---- (17.5) 
f = c(0, 1)
(-1)^(1-f)
(-1)^(1) * (-1)^(-f)
-(-1)^(-f)
-(-1)^f
(-1)^(1-f) == -(-1)^f


# ---- (17.9) Deutsch-Josza
multikron(H %*% q0, H %*% q0, H %*% q0, H %*% q0)
multikron(H %*% q0, H %*% q0, H %*% q0, H %*% q1)
multikron(H %*% q0, H %*% q0, H %*% q1, H %*% q0)
multikron(H %*% q0, H %*% q0, H %*% q1, H %*% q1)
multikron(H %*% q0, H %*% q1, H %*% q0, H %*% q0)
multikron(H %*% q0, H %*% q1, H %*% q0, H %*% q1)
multikron(H %*% q0, H %*% q1, H %*% q1, H %*% q0)
multikron(H %*% q0, H %*% q1, H %*% q1, H %*% q1)
multikron(H %*% q1, H %*% q0, H %*% q0, H %*% q0)
multikron(H %*% q1, H %*% q0, H %*% q0, H %*% q1)
multikron(H %*% q1, H %*% q0, H %*% q1, H %*% q0)
multikron(H %*% q1, H %*% q0, H %*% q1, H %*% q1)
multikron(H %*% q1, H %*% q1, H %*% q0, H %*% q0)
multikron(H %*% q1, H %*% q1, H %*% q0, H %*% q1)
multikron(H %*% q1, H %*% q1, H %*% q1, H %*% q0)
multikron(H %*% q1, H %*% q1, H %*% q1, H %*% q1)

# left and right sides of eq. (17.9)
n = 10
n = 8
n = 6
n = 4
for (ix in 1:2^n - 1)
{
    x = as.bits(ix, n)
    right = c()
    for(iy in 1:2^n - 1)
    {
        y = as.bits(iy, n)
        right = c(right, (-1)^sum(x*y))
    }
    right = t(t(right)) / sqrt(2^n)
    
    left = list()
    for(i in 1:length(x))
        left[[i]] = if (x[i] == 0) H %*% q0 else H %*% q1
    
    left = multikron(left)
    stopifnot(max(Mod(left - right)) < 1e-12)
}


# ----- AD lecture 18

# Bernstein-Vazirani, (18.4)
n = 8
n = 4
#allbits = expand.grid(rep(list(0:1), n))
for(si in 1:2^n - 1)
{
    # s = unlist(allbits[1 + si,], use.names = F)
    s = as.bits(si, n)
    ket_s = as.qubit(si, n)
    inner = c()
    for(xrow in 1:2^n - 1)
    {
        x = as.bits(xrow, n)
        inner = c(inner, (-1)^sum(s*x))
    }
    inner = t(t(inner)) / sqrt(2^n)
    Hn = multikron(rep(list(H), n))
    stopifnot(max(Mod(Hn %*% inner - ket_s)) < 1e-12)
}


# ad Birthday paradox 
# if there were N days in the year, then you’d need about √N 
# people in the room for a likely birthday collision
# N - number of days, n - number of people
N = 365
n = 23
# ≈≈
(1 - 1/N)^(n*(n-1)/2)
.5
# ≈≈
(1-1/N)^(n^2/2)
.5
# ≈≈
(1-1/N)^(n^2)
(.5)^2
# ≈≈
log((1-1/N)^(n^2))
-2*log(2)
# ≈≈
n^2*log(1-1/N)
-2*log(2)
# ≈≈
n^2
-2*log(2)/log(1-1/N)
# ≈≈
n^2
2*log(2)/log(N/(N-1))
# ≈≈
# N --> N+1
n^2
2*log(2)/log((N+1)/N)
# ≈≈
n^2
2*log(2)/log(1+1/N)
# ≈≈
# when z is small: log(1 + z) ≈ z
# (1/N) is small => log(1+1/N) ≈ 1/N
n^2
2*log(2)*N
# ≈≈
n
sqrt(2*log(2)) * sqrt(N)
# ≈≈
n
sqrt(N)

# ---------------
# ad https://qiskit.org/textbook/ch-algorithms/deutsch-jozsa.html#3.-Creating-Quantum-balanceds--
# Deutsch-Jozsa
# "One of the ways we can guarantee our circuit is balanced is by performing a CNOT for each qubit..."
#

# oracle bit - leftmost bit (bit0, bit ids: 3210)
CX03l = multikron(I, I, SW) %*% multikron(I, SW, I) %*% multikron(SW %*% CX %*% SW, I, I) %*% multikron(I, SW, I) %*% multikron(I, I, SW) 
CX13l = multikron(I, SW, I) %*% multikron(SW %*% CX %*% SW, I, I) %*% multikron(I, SW, I)
CX23l = multikron(SW %*% CX %*% SW, I, I)
balanced_left = CX23l %*% CX13l %*% CX03l

# oracle bit - rightmost bit (bit3, bit ids: 3210)
CX03r = multikron(SW, I, I) %*% multikron(I, SW, I) %*% multikron(I, I, CX) %*% multikron(I, SW, I) %*% multikron(SW, I, I) 
CX13r = multikron(I, SW, I) %*% multikron(I, I, CX) %*% multikron(I, SW, I)
CX23r = multikron(I, I, CX)
balanced_right = CX23r %*% CX13r %*% CX03r 

stopifnot(all(balanced_left == REV4 %*% balanced_right %*% REV4))

all(CX03r == selSW(1, 3, 4) %*% multikron(I, I, CX) %*% selSW(1, 3, 4))
all(CX03l == selSW(2, 0, 4) %*% multikron(SW %*% CX %*% SW, I, I) %*% selSW(2, 0, 4))

all(multikron(SW %*% CX %*% SW, I, I) == multikron(SW, I, I) %*% multikron(CX, I, I) %*% multikron(SW, I, I) )


balanced_left %*% as.qubit(0, 4)
balanced_left %*% as.qubit(3, 4)
balanced_left %*% as.qubit(5, 4)
balanced_left %*% as.qubit(6, 4)

balanced_left %*% as.qubit(1, 4)
balanced_left %*% as.qubit(2, 4)
balanced_left %*% as.qubit(4, 4)
balanced_left %*% as.qubit(7, 4)

wrap = multikron(I, X, I, I)
wrap = multikron(I, X, I, I)
wrap = multikron(X, X, X, X)

# left vs right
if (F)
{
    # right
    balanced = balanced_right
    balancedwrp = wrap %*% balanced %*% wrap
    const0 = multikron(I, I, I, I)
    const1 = multikron(X, X, X, I)
    
    full_dj = multikron(H, H, H, I) %*% balanced %*% multikron(H, H, H, H %*% X) %*% as.qubit(0, 4)
    full_dj = multikron(H, H, H, I) %*% balancedwrp %*% multikron(H, H, H, H %*% X) %*% as.qubit(0, 4)
    full_dj = multikron(H, H, H, I) %*% const0 %*% multikron(H, H, H, H %*% X) %*% as.qubit(0, 4)
    full_dj = multikron(H, H, H, I) %*% const1 %*% multikron(H, H, H, H %*% X) %*% as.qubit(0, 4)

    # meas 
    # P( |000〉)
    Mod(t(Conj(as.qubit(0, 4))) %*% full_dj)^2 + Mod(t(Conj(as.qubit(1, 4))) %*% full_dj)^2
    # P( |111〉)
    Mod(t(Conj(as.qubit(14, 4))) %*% full_dj)^2 + Mod(t(Conj(as.qubit(15, 4))) %*% full_dj)^2
    
} else
{   
    # left
    balanced = balanced_left
    balancedwrp = wrap %*% balanced %*% wrap
    const0 = multikron(I, I, I, I)
    const1 = multikron(I, X, X, X)
    
    full_dj = multikron(I, H, H, H) %*% balanced %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    full_dj = multikron(I, H, H, H) %*% balancedwrp %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    full_dj = multikron(I, H, H, H) %*% const0 %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    full_dj = multikron(I, H, H, H) %*% const1 %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)

    # meas
    # P( |000〉)
    Mod(t(Conj(as.qubit(0, 4))) %*% full_dj)^2 + Mod(t(Conj(as.qubit(8, 4))) %*% full_dj)^2
    # P( |111〉)
    Mod(t(Conj(as.qubit(7, 4))) %*% full_dj)^2 + Mod(t(Conj(as.qubit(15, 4))) %*% full_dj)^2
    
}

all(kronecker(X, I) %*% q10 == q00)
all(kronecker(I, X) %*% q10 == q11)


#---------
# https://qiskit.org/textbook/ch-algorithms/bernstein-vazirani.html
#
# bits: o210
all(Z %*% H %*% q0 == qm)
bv_func = multikron(I, I, I, I) # 000
bv_func = CX23l %*% CX13l %*% CX03l # 111
bv_func = CX23l %*% CX03l # 101
bv_func = CX23l %*% CX13l # 110
bv_func = CX13l %*% CX03l # 011
{
    full_bv = multikron(I, H, H, H) %*% bv_func %*% multikron(H %*% X, H, H, H) %*% as.qubit(0, 4)
    # meas bits 210
    res = matrix(0, 8, 1)
    for (val in 0:7)
        res[1+val, 1] = Mod(t(Conj(as.qubit(val, 4))) %*% full_bv)^2 + Mod(t(Conj(as.qubit(8+val, 4))) %*% full_bv)^2
    round(res, 10)
    as.bits(which.max(res) - 1, 3)
}

# ad Simon
if (FALSE)
{
    n=10
    t=1:2^(n/2)
    plot(t, (1-(t*(t-1)/2/(2^n-1)))^t, type='l', col ='blue')
}



# ----- AD lecture 20
if (FALSE)
{
    Q = 2^6
    plot(exp(2i*pi/Q)^(0:(Q-1)), xlim = c(-1,1), ylim = c(-1,1), col = "blue")
}




SW %*% CX %*% SW 

qft = function(n)
{
    w = exp(2i*pi/n) 
    res = diag(0, n)
    for(i in 0:(n-1))
        for(j in 0:(n-1))
            res[1+i, 1+j] = w^(i*j)
    res / sqrt(n)
}

stopifnot(max(Mod(qft(2) - H)) < 1e-12)
B2 = matrix(c(1,0,0,1i), 2)
qft(4)

 
