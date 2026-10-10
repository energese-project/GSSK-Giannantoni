# Numerics for the First Fundamental Equation's closed form (srs FR-MOP-001, numerics N1).
#   alpha(t) = [ ((a+bt)^q - a^q) / (b (p+k)) ]^k ,  q = (p+k)/k
# Naive evaluation cancels as b t / a -> 0. The unified form
#   S = (a^{p/k} t / k) * E(x),  x = b t / a,  E(x) = expm1(q log1p x) / (q x),  E(0) = 1
# is accurate there and continuous through b = 0. Reference: mpmath-free exact series check.
import math
from fractions import Fraction

def naive(a, b, p, k, t):
    q = (p + k) / k
    return (((a + b*t)**q - a**q) / (b*(p + k)))**k

def unified(a, b, p, k, t):
    q = (p + k) / k
    x = b*t/a
    E = 1.0 if x == 0 else math.expm1(q*math.log1p(x)) / (q*x)
    return (a**(p/k) * t / k * E)**k

def reference(a, b, p, k, t, terms=40):
    # (1+x)^q - 1 = sum_{n>=1} binom(q,n) x^n, exact rational arithmetic for rational inputs
    a, b, p, k, t = map(Fraction, (a, b, p, k, t))
    q = (p + k) / k
    x = b*t/a
    s, c = Fraction(0), Fraction(1)
    for n in range(1, terms):
        c = c * (q - n + 1) / n
        s += c * x**n
    S = float(a)**float(q) * float(s) / float(b*(p + k))
    return S**float(k)

for t in (1e-2, 1e-5, 1e-8, 1e-11):
    a, b, p, k = 1.0, 0.25, 1.0, 2.0
    r = reference(a, b, p, k, t)
    print(f"t={t:.0e}  naive rel err {abs(naive(a,b,p,k,t)-r)/r:.1e}   unified rel err {abs(unified(a,b,p,k,t)-r)/r:.1e}")
print("continuity at b -> 0 (a=1,p=1,k=2,t=1): ",
      [f"{unified(1.0, b, 1.0, 2.0, 1.0):.15f}" for b in (1e-6, 1e-12, 0.0)],
      " exact b=0:", f"{(1.0**0.5*1/2)**2:.15f}")
