/* idc.h — single-variable Incipient Differential Calculus.
 *
 * docs/requirements/icd.md IF-API-002; srs.md §2.1. Each function implements
 * one requirement, named beside it, from the source equation it cites. The
 * derivative of a single function given pointwise is the incipient one,
 *
 *     (d~/dt)^n f = (f'/f)^n f            [02 Eq 14.9.5], [23 Eq 5.5.2]
 *
 * and of a superposition of exponential terms it is taken termwise (srs §1,
 * PLAN R16).
 *
 * Status of record: this header is the incipient calculus. The kernel's
 * "incipient" method is not (PLAN.md §2 E2, ADR 0022).
 */

#ifndef GIA_IDC_H
#define GIA_IDC_H

#include <complex.h>

#include "gia_status.h"

/* FR-IDC-002. (f'/f)^n * f for f != 0 and integer n >= 0, in C. The power is
 * taken by repeated multiplication, never cpow, so its sign and branch are
 * exact (numerics.md §1). f == 0 is GIA_E_DOMAIN; n < 0 or out == NULL is
 * GIA_E_ARG; a non-finite result is GIA_E_RANGE. */
gia_status gia_idc_of(double complex f, double complex df, int n,
                      double complex *out, const char **why);

/* A caller-owned coefficient function a(t). It must be finite on [0, t_max];
 * a non-finite value is GIA_E_DOMAIN. */
typedef double complex (*gia_cfn)(double t, void *ctx);

/* FR-IDC-006 — the second-order incipient LDE with variable coefficients,
 *
 *     f~'' + a1(t) f~' + a0(t) f = 0,   f(0) = f0,   f~'(0) = f1,
 *
 * solved as f = sum_i c_i exp(int_0^t r_i), with r_i(t) the roots of
 * r^2 + a1 r + a0 = 0 labelled by continuity in t ([06 Eq 3.3-3.6],
 * [09 Eq 3, 7], [10 Eq 8.1, 10.1]; numerics.md N3). The derivative is
 * termwise (PLAN R16). Where the discriminant vanishes on the whole interval
 * the solution is the one family c e^{int r}: solved when f1 = r(0) f0,
 * refused (GIA_E_DOMAIN) otherwise, and [06 Eq 3.7] is not used (PLAN X12).
 * Roots that coincide at t = 0 but not everywhere leave the constants
 * undetermined (GIA_E_DOMAIN); roots that collide at an isolated time
 * t* > 0 are GIA_E_CONVERGENCE, with t* written to *t_fail when t_fail is
 * not NULL (the one output written on failure: it is the diagnosis).
 *
 * The solution holds a1, a0 and ctx without copying; ctx must outlive it.
 * Free it with gia_lde2_free. */
typedef struct gia_lde2_sol gia_lde2_sol;

gia_status gia_lde2_solve(gia_cfn a1, gia_cfn a0, void *ctx,
                          double complex f0, double complex f1, double t_max,
                          gia_lde2_sol **sol, double *t_fail, const char **why);

/* f(t), and the traditional residual f'' + a1 f' + a0 f of the same terms
 * (either pointer may be NULL). t must lie in [0, t_max], else GIA_E_ARG.
 * Every exponential is guarded (numerics.md N6): GIA_E_RANGE, never inf. */
gia_status gia_lde2_eval(const gia_lde2_sol *sol, double t,
                         double complex *f, double complex *trad_residual,
                         const char **why);

/* The terms at t: constants c[i], roots r[i](t), and E[i] = exp(int_0^t r_i).
 * f = sum c E, f~' = sum c r E. Exposed so that a test can form the defining
 * equation's residual itself (vv-plan.md §2, "residual first"). A one-family
 * solution reports one term, with c[1] = 0 and r[1] = r[0], E[1] = E[0]. */
gia_status gia_lde2_terms(const gia_lde2_sol *sol, double t,
                          double complex c[2], double complex r[2],
                          double complex E[2], const char **why);

void gia_lde2_free(gia_lde2_sol *sol);

/* FR-IDC-007 — the binary function [02 Eq 14.7.1-14.7.7], [06 Eq 3.10-3.15]:
 *
 *     f' + A f^(1/2) + B f = 0,    two branches sigma = 1, 2,
 *
 * solved per PLAN R9 / numerics.md N4: u1, u2 are the roots of
 * u^2 + A u + B = 0 (one characteristic serves both branches; the printed
 * "distinct exponents per branch" is erratum X7), each branch is
 * f_s = c_s1 e^{u1^2 t} + c_s2 e^{u2^2 t} with half-derivative
 * sum_i c_si u_i e^{u_i^2 t}, and the constants come from
 * [[1, 1], [u1, u2]] c_s = (f_s(0), f_s^(1/2)(0)). u1 = u2 is defined by no
 * source and is refused (GIA_E_UNSUPPORTED). */
typedef struct { double complex u[2]; double complex c[2][2]; } gia_binary_sol;

gia_status gia_binary_solve(double complex A, double complex B,
                            const double complex f0[2], const double complex fhalf0[2],
                            gia_binary_sol *sol, const char **why);

/* f_s(t) and f_s^(1/2)(t) for both branches (fhalf may be NULL). Returns a
 * status rather than nothing so that an overflowing exponential is
 * GIA_E_RANGE, never inf (NFR-NUM-003). */
gia_status gia_binary_eval(const gia_binary_sol *sol, double t,
                           double complex f[2], double complex fhalf[2],
                           const char **why);

#endif /* GIA_IDC_H */
