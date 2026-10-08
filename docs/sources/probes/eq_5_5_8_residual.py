# Residual of (d~/dt)^k f = (f'/f)^k * f against beta = (a+b t)^p  (Eq. 5.5.2)
import sys; a,b,p,k = map(float, sys.argv[1:5]) if len(sys.argv) == 5 else (1.0, 0.25, 1.0, 2.0)  # usage: a b p k
beta = lambda t: (a+b*t)**p
def mine(t):   # derived from 5.5.6: f = ( (1/k) * int_0^t beta^(1/k) )^k
    I = ((a+b*t)**((p+k)/k) - a**((p+k)/k)) * k/(b*(p+k))
    return (I/k)**k
def p557(t):   # as printed 5.5.7: f = (1/k) * { int_0^t beta^(1/k) }^k
    I = ((a+b*t)**((p+k)/k) - a**((p+k)/k)) * k/(b*(p+k))
    return (1/k)*I**k
def p558(t):   # as printed 5.5.8 last line
    return (1/(k*b)) * ((k/(p+k))*(a+b*t)**(p/k+1))**k
def resid(f,t,h=1e-6):
    d=(f(t+h)-f(t-h))/(2*h); v=f(t)
    return (d/v)**k*v - beta(t)
for name,f in [("derived",mine),("printed 5.5.7",p557),("printed 5.5.8",p558)]:
    print(name, [round(resid(f,t),6) for t in (0.5,1.0,2.0,4.0)])
