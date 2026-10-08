# Giannantoni & Zoli 2009, Eq. 10: incipient Taylor series truncated at n,
# f*(t0+D) = f0 * sum_{k=0..n} (a*D)^k / k!,  a = f'(t0)/f(t0)
from math import factorial
def inc(f0, f1, D, n=2):
    a = f1 / f0
    return f0 * sum((a*D)**k / factorial(k) for k in range(n+1))
print("Eq16 max scen   paper 16.4   got", inc(0.4, 0.32, 10))
print("Eq17 min scen   paper 3.01   got", inc(0.4, 0.11, 10))
print("Eq21 sea L1     paper 178.0  got", inc(18.0, 6.0, 10))
for t0, paper in ((2.0, 154.3), (1.8, 172.06)):
    print(f"Eq19 tau0={t0}  paper {paper}  got", inc(6.0*t0, 6.0, 10-t0))
for t0 in (2.0, 1.8):
    print(f"Eq22 tau0={t0} paper 15-17   got", inc(0.6*t0, 0.6, 10-t0))
# derived percentages as printed
print("max: net", inc(.4,.32,10)-.4, "vs 6.4 ->", (inc(.4,.32,10)-.4)/6.4-1)
print("min: net", inc(.4,.11,10)-.4, "(paper says 1.91 net, 73%)", "3.01-1.1 =", inc(.4,.11,10)-1.1)
