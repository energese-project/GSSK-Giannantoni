# Giannantoni 2006 Eq 3.3-3.7: which reading of the incipient derivative makes the
# printed solutions solve  f~'' + a1 f~' + a0 f = 0 ?  Constant coefficients.
#   pointwise : (d~/dt)^n f = (f'/f)^n f            [02 Eq 14.9.5], [06 Eq 3.2]
#   termwise  : linear extension over exponential terms, (d~/dt)^n e^{ut} = u^n e^{ut}
#   classical : the ordinary derivative
import math
t = 0.8
def residuals(terms, a1, a0):
    # terms: list of (c, u) for f = sum c e^{u t}
    f   = sum(c*math.exp(u*t) for c,u in terms)
    f1  = sum(c*u*math.exp(u*t) for c,u in terms)          # ordinary f'
    f2  = sum(c*u*u*math.exp(u*t) for c,u in terms)        # ordinary f''
    classical = f2 + a1*f1 + a0*f
    termwise  = sum(c*(u*u + a1*u + a0)*math.exp(u*t) for c,u in terms)
    b = f1/f
    pointwise = (b*b + a1*b + a0)*f
    return classical, termwise, pointwise
# Distinct roots 1 and -2: a1 = 1, a0 = -2.  Eq 3.6 superposition.
print("Eq 3.6, roots {1,-2}, f = e^t + e^-2t:  classical %.3g  termwise %.3g  pointwise %.3g"
      % residuals([(1,1.0),(1,-2.0)], 1.0, -2.0))
# Double root alpha = 0.7: a1 = -2 alpha, a0 = alpha^2.
al = 0.7
print("double root, f = e^{at}:                 classical %.3g  termwise %.3g  pointwise %.3g"
      % residuals([(1,al)], -2*al, al*al))
# Eq 3.7 with constant alpha: y2 = int_0^t e^{alpha (t - xi)} dxi = (e^{alpha t} - 1)/alpha
print("Eq 3.7 y2 = (e^{at}-1)/a:                classical %.3g  termwise %.3g  pointwise %.3g"
      % residuals([(1/al,al),(-1/al,0.0)], -2*al, al*al))
# Classical second solution t e^{at}, for comparison (not a sum of exponentials, so no termwise form)
f = t*math.exp(al*t); f1 = (1+al*t)*math.exp(al*t); f2 = (2*al + al*al*t)*math.exp(al*t)
b = f1/f
print("classical t e^{at}:                      classical %.3g  pointwise %.3g"
      % (f2 - 2*al*f1 + al*al*f, (b*b - 2*al*b + al*al)*f))
