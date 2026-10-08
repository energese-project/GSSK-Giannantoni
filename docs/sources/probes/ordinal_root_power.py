# Is the EQS root r_l = B + C j + C k (Giannantoni 2023 Eq 7.4, De Moivre form)
# an (N-1)-th root of unity under the literal relational product (Eq 5.1.3-5.1.5),
# with powers taken left to right?  Angle chosen so sqrt(2)*psi_l = 2*pi*l/(N-1).
import math
T = {(0,0):(1,0),(0,1):(1,1),(0,2):(1,2),(1,0):(1,1),(1,1):(-1,0),(1,2):(1,2),
     (2,0):(1,2),(2,1):(1,2),(2,2):(-1,0)}
def prod(x, y):
    out = [0.0]*3
    for a in range(3):
        for b in range(3):
            c, e = T[(a,b)]; out[e] += c*x[a]*y[b]
    return out
for N in (3, 4, 5, 7):
    for l in range(1, N-1):
        r_ang = 2*math.pi*l/(N-1)
        r = [math.cos(r_ang), math.sin(r_ang)/math.sqrt(2), math.sin(r_ang)/math.sqrt(2)]
        p = r[:]
        for _ in range(N-2): p = prod(p, r)          # ((r r) r) ... left to right
        q = r[:]
        for _ in range(N-2): q = prod(r, q)          # r (r (r ...)) right to left
        print(f"N={N} l={l}  left: {[round(v,6) for v in p]}  right: {[round(v,6) for v in q]}")
