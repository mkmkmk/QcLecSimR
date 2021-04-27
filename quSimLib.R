# ----------------------
# QuSimR "library" code
#
#
#   Author:  Mariusz Krej
#
# ----------------------
epsilon = 1e-12    


q0 = t(t(1:0))  # |0〉
q1 = t(t(0:1))  # |1〉

H = matrix(c(1,1,1,-1), nrow=2) / sqrt(2)

qp = H %*% q0  # |+〉
qm = H %*% q1  # |-〉

q00 = kronecker(q0, q0)
q01 = kronecker(q0, q1)
q10 = kronecker(q1, q0)
q11 = kronecker(q1, q1)

qpp = kronecker(qp,qp) #  |++〉
qpm = kronecker(qp,qm) #  |+-〉
qmp = kronecker(qm,qp) #  |-+〉
qmm = kronecker(qm,qm) #  |--〉

I = diag(1, 2)
NOT = 1 - I
CNOT = kronecker(q0 %*% t(q0), I) + kronecker(q1 %*% t(q1), NOT)
CX = CNOT
ROT = function(th) matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2)

NOTI = kronecker(NOT, I)
H2 = kronecker(H, H)
IH = kronecker(I, H)
HI = kronecker(H, I)
SWAP = matrix(c(1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1), 4)
SW = SWAP

Z = matrix(c(1, 0, 0, -1), 2, 2)
Y = matrix(c(0, 1i, -1i, 0), 2, 2)
X = NOT
S = matrix(c(1, 0, 0, 1i), 2, 2)
Tg = matrix(c(1, 0, 0, exp(.25i*pi)), 2, 2)

XI = kronecker(NOT, I) 
ZI = kronecker(Z, I)
II = kronecker(I, I)

bell = CNOT %*% HI %*% q00

requal = function(a, b, dig=10) round(a, dig) == round(b, dig)

toBloch = function (qbit)
{
    stopifnot(sum(Mod(qbit)^2) < 1 + 1e-10)
    # TODO if not, rotate |0〉
    stopifnot(abs(Im(qbit[1])) < 1e-10)
    th = 2 * acos(qbit[1])
    #sin = (1 - qbit[1]^2)^.5
    #stopifnot(abs(sin - sin(th/2)) < 1e-9)
    #phi = if (abs(th) > 1e-9) acos(Re(qbit[2]) / sin) else 0
    phi = Arg(qbit[2])
    round(c(th, phi) / pi, 8) 
}

toStateTb = function(th, phi)
    round(q0 %*% t(cos(th/2)) + q1 %*% t(exp(1i*phi)*sin(th/2)), 10)

toState = function(bloch)
    toStateTb(bloch[1], bloch[2])

# post partial measurement state
# measBasis - measurement basis-state list
postmeas = function(state, ...)
{
    measBasis = list(...)
    if(is.list(measBasis[[1]]))
        measBasis = unlist(measBasis, recursive = FALSE)
    res = 0
    for (basis in measBasis)
        res = res + (Conj(t(basis)) %*% state)[1,1] * basis
    res / sqrt(sum(Mod(res)^2)) # eq. (5.1)
}

# mutliple kronecker
multikron = function(...)
{
    inp = list(...)
    if(is.list(inp[[1]]))
        inp = unlist(inp, recursive = FALSE)
    res = inp[[1]]
    if (length(inp) < 2)
        return(inp[[1]])
    for (ii in 2:length(inp))
        res = kronecker(res, inp[[ii]])
    res
}

# tests
stopifnot(all(multikron(I) == I))
stopifnot(all(multikron(I,H) == IH))
stopifnot(all(multikron(H,I,I) == kronecker(HI, I)))
stopifnot(all(multikron(I, H, I, I) == kronecker(I, kronecker(HI, I))))

GHZ = (multikron(q0, q0, q0) + multikron(q1, q1, q1)) / sqrt(2)

as.qubit = function(inp, nbits = NA)
{
    stopifnot(inp >= 0)
    if (is.na(nbits))
        nbits = max(1, ceiling(log2(1 + inp)))
    #if (inp == 2^nbits) 
    #    nbits = nbits + 1
    stopifnot(inp < 2^nbits)
    ket = matrix(0, 2^nbits, 1)
    ket[1 + inp] = 1
    ket
}

# tests
stopifnot(as.qubit(0) == q0)
stopifnot(as.qubit(1) == q1)
stopifnot(as.qubit(3) == q11)
stopifnot(as.qubit(4) == kronecker(q1, q00))
stopifnot(as.qubit(0, 2) == q00)
stopifnot(as.qubit(1, 2) == q01)
stopifnot(as.qubit(127)[1 + 127] == 1)

as.bits = function(inp, nbits = NA)
{
    stopifnot(inp >= 0)
    if (is.na(nbits))
        nbits = max(1, ceiling(log2(1 + inp)))
    #if (inp == 2^nbits) inp = inp + 1
    stopifnot(inp < 2^nbits)
    rev(as.integer(intToBits(inp)[1:nbits]))
}

# tests
stopifnot(as.bits(0) == c(0))
stopifnot(as.bits(1) == c(1))
stopifnot(as.bits(2) == c(1, 0))
stopifnot(as.bits(3) == c(1, 1))
stopifnot(as.bits(4) == c(1, 0, 0))
stopifnot(as.bits(7) == c(1, 1, 1))
stopifnot(as.bits(8) == c(1, 0, 0, 0))
stopifnot(as.bits(15) == c(1, 1, 1, 1))
stopifnot(as.bits(16) == c(1, 0, 0, 0, 0))
stopifnot(length(as.bits(2^16 - 1)) == 16)
stopifnot(length(as.bits(2^16)) == 17)

REV4 = multikron(I, SW, I) %*% multikron(SW, SW) %*% multikron(I, SW, I) %*% multikron(SW, SW)
REV3 = multikron(I, SW) %*% multikron(SW, I) %*% multikron(I, SW)
all(REV3 == multikron(SW, I) %*% multikron(I, SW) %*% multikron(SW, I))
stopifnot(REV4 %*% REV4 == multikron(rep(list(I), 4)))
stopifnot(REV3 %*% REV3 == multikron(rep(list(I), 3)))

# reverse bits
revBits = function(n)
{
    res = diag(0, 2^n)
    for(x in 0:(2^n-1))
        res[, 1+x] = as.qubit(sum(as.bits(x,n) * 2^(0:(n-1))), n)
    res
}

# tests
stopifnot(revBits(2) == SW)
stopifnot(revBits(3) == REV3)
stopifnot(revBits(4) == REV4)

# swap selected bits 
# bit indexes: ..543210
swap = function(b1, b2, n)
{
    stopifnot(b1 < n && b2 < n && b1 >= 0 && b2 >= 0)
    b1 = n - 1 - b1
    b2 = n - 1 - b2
    res = diag(0, 2^n)
    for(x in 0:(2^n-1))
    {
        bts = as.bits(x, n)
        mem = bts[1+b1]
        bts[1+b1] = bts[1+b2]
        bts[1+b2] = mem
        res[, 1+x] = as.qubit(sum(2^((n-1):0) * bts), n)
    }
    res
}

# tests
stopifnot(swap(1, 0, 2) == SW)
stopifnot(swap(1, 0, 3) == kronecker(I, SW))
stopifnot(swap(1, 2, 4) == multikron(I, SW, I))
stopifnot(swap(3, 2, 4) == multikron(SW, I, I))
stopifnot(swap(1, 2, 4) == swap(2, 1, 4))


ctr_gate = function(gate)
{
    n = nrow(gate)
    stopifnot(abs(2^round(log2(n)) - n) < epsilon)
    I_n = diag(1, n)
    zeros = diag(0, n)
    rbind(cbind(I_n, zeros), cbind(zeros, gate))
}
stopifnot(Mod(ctr_gate(gate = NOT) - CNOT) < epsilon)
