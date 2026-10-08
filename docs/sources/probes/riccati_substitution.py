# Giannantoni 2006 Eq 3.16-3.18: which substitution maps the Riccati
#   f' + Q f + R f^2 = P   onto   R y'' - (R' - Q R) y' - P R^2 y = 0 ?
# Constant Q=1, R=1, P=2 so (3.18) is y'' + y' - 2 y = 0.
import math
h = 1e-4
y  = lambda t: math.exp(t) + math.exp(-2*t)             # solves (3.18)
d  = lambda g, t: (g(t+h) - g(t-h)) / (2*h)
d2 = lambda g, t: (g(t+h) - 2*g(t) + g(t-h)) / h**2
f_std = lambda t: d(y, t) / y(t)                        # standard: f = y'/(R y)
ric = lambda f, t: d(f, t) + f(t) + f(t)**2 - 2
print("standard f = y'/(R y): Riccati residual",
      max(abs(ric(f_std, t)) for t in (0.1, 0.5, 1.0, 2.0)))
# printed 3.17: y = f'/(f R). Take the Riccati solution f_std and build y from it.
y_pr = lambda t: d(f_std, t) / f_std(t)
lin = lambda g, t: d2(g, t) + d(g, t) - 2*g(t)
print("printed y = f'/(f R): (3.18) residual",
      max(abs(lin(y_pr, t)) for t in (0.1, 0.5, 1.0, 2.0)))
