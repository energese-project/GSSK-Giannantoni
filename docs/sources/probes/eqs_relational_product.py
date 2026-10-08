# Does Giannantoni 2023 EQS (7.1-7.4) follow from the literal relational
# product table (5.1.3-5.1.5) applied to a De Moivre root of unity?
import math, random
# basis index: 0=i (real unit, +1), 1=j, 2=k ; table[a][b] = (coef, basis)
T = {(0,0):(1,0),(0,1):(1,1),(0,2):(1,2),
     (1,0):(1,1),(1,1):(-1,0),(1,2):(1,2),
     (2,0):(1,2),(2,1):(1,2),(2,2):(-1,0)}
def prod(x, y):
    out = [0.0,0.0,0.0]
    for a in range(3):
        for b in range(3):
            c, e = T[(a,b)]; out[e] += c*x[a]*y[b]
    return out
random.seed(1)
worst = 0.0
for _ in range(1000):
    psi = random.uniform(-3,3); S,F,Th = (random.uniform(-2,2) for _ in range(3))
    r = math.sqrt(2)*psi
    B = math.cos(r); C = math.sin(r)/math.sqrt(2)          # Eq. 7.4
    root = [B, C, C]                                        # De Moivre: cos r + (j+k) sin r / sqrt2
    got = prod(root, [S,F,Th])
    want = [B*S - C*(F+Th), B*F + C*S, B*Th + C*S + C*(F+Th)]  # brackets of 7.1.1, 7.2, 7.3
    worst = max(worst, max(abs(g-w) for g,w in zip(got,want)))
print("max |table-product - EQS bracket| over 1000 draws:", worst)
# Same check with an associative alternative (jk+kj=0, jk dropped): does it match 7.3?
T2 = dict(T); T2[(1,2)] = (0,2); T2[(2,1)] = (0,2)
T, Tsave = T2, T
got = prod([B,C,C],[S,F,Th]); print("anticommuting variant k-part:", got[2], " EQS 7.3 wants:", B*Th + C*S + C*(F+Th))
# associativity witness for the literal table
T = Tsave; j=[0,1,0]; k=[0,0,1]
print("(j*j)*k =", prod(prod(j,j),k), "  j*(j*k) =", prod(j,prod(j,k)))
