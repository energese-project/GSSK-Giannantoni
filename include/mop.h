/* mop.h — the Maximum Ordinality Principle: the First Fundamental Equation.
 *
 * docs/requirements/icd.md IF-API-004; srs.md §2.3; ADR 0019. For each couple
 * (i, j) of components the First Equation [23 Eq 5.4.2]
 *
 *     (d~/dt)^k alpha_ij = beta_ij,      (d~/dt)^k f = (f'/f)^k f  [23 Eq 5.5.2]
 *
 * is solved with alpha(0) = 0 by the solution derived from [23 Eq 5.5.6],
 *
 *     alpha(t) = { (1/k) int_0^t beta^{1/k} }^k,
 *
 * never by the printed [23 Eq 5.5.7-5.5.8] (PLAN X1). Couples are independent;
 * their collection, with a zero diagonal, is the Matrioska [23 Eq 5.6.1].
 *
 * This unit does not include gssk.h (NFR-SEP-001) and does not call the
 * harmony constructor (FR-HAR-003).
 */

#ifndef GIA_MOP_H
#define GIA_MOP_H

#include <complex.h>

#include "gia_status.h"
#include "relational.h"

/* The cardinality k as a reduced fraction num/den, den >= 1. */
typedef struct { int num, den; } gia_rational;

typedef enum { GIA_BETA_NONE, GIA_BETA_AFFINE, GIA_BETA_SAMPLES } gia_beta_kind;

/* A couple's boundary condition. AFFINE: beta(t) = (a + b t)^p, principal
 * branch. SAMPLES: the piecewise-linear interpolant of (t[i], v[i]),
 * i = 0..n-1, t strictly increasing from t[0] = 0; the arrays are the
 * caller's. NONE: the couple is unrelated. */
typedef struct {
    gia_beta_kind          kind;
    double complex         a, b;
    double                 p;
    int                    n;
    const double          *t;
    const double complex  *v;
} gia_beta;

/* FR-MOP-001, FR-MOP-002, FR-MOP-004 — one couple at time t >= 0.
 *
 * AFFINE uses numerics.md N1, S(t) = (a^{p/k} t / k) E(x), x = b t / a,
 * E(x) = expm1(q log1p(x)) / (q x), q = (p + k)/k, alpha = S^k: no
 * cancellation near t = 0 and continuous through b = 0. SAMPLES uses N2:
 * adaptive Gauss-Kronrod 7-15 per segment on |beta|^{1/k} e^{i theta/k}, theta
 * continued from the principal argument at t = 0.
 *
 * Domain (FR-MOP-002): integer k >= 1 with real or complex a, b; a non-integer
 * k only with real a > 0, b >= 0 (and, for samples, a real positive S). Refused
 * with GIA_E_DOMAIN: any other k; a + b t = 0 somewhere on [0, t]; p = -k with
 * b != 0; a sample or segment through beta = 0; t past the last sample.
 * Overflow is GIA_E_RANGE (N6); an unmet quadrature tolerance
 * GIA_E_CONVERGENCE (N2). */
gia_status gia_mop_couple(const gia_beta *beta, gia_rational k, double t,
                          double complex *alpha, const char **why);

/* FR-MOP-003 — the Matrioska: N x N, row-major. a[i*N + j] is alpha_ij(t)
 * where related[i*N + j] is 1; unrelated couples (GIA_BETA_NONE) have
 * related 0 and an unspecified a. The diagonal is related 0 and a = 0. */
typedef struct { int N; double complex *a; unsigned char *related; } gia_matrioska;

/* Solves every off-diagonal couple of beta (N*N, row-major, diagonal ignored).
 * Allocates out->a and out->related; free with gia_matrioska_free. On any
 * couple's failure nothing is allocated and its status is returned. */
gia_status gia_mop_solve(int N, const gia_beta *beta, gia_rational k, double t,
                         gia_matrioska *out, const char **why);

void gia_matrioska_free(gia_matrioska *m);

/* NFR-NUM-005 — the adaptive Gauss-Kronrod 7-15 integrator N2 uses, exposed
 * so that its error estimate and its refusal can be tested directly. Integrates
 * g over [lo, hi] to |err| <= tol (1 + |I|), bisecting at most `max_depth`
 * levels deep; returns the integral and the summed error estimate, or
 * GIA_E_CONVERGENCE when the tolerance is not met within the depth. */
typedef double complex (*gia_quad_fn)(double t, void *ctx);
gia_status gia_quad_gk15(gia_quad_fn g, void *ctx, double lo, double hi, double tol,
                         int max_depth, double complex *integral, double *err,
                         const char **why);

/* FR-MOP-007 — a relational-valued couple, beta = beta_i i (+) beta_j j (+)
 * beta_k k with each component a real boundary condition. For k = 1 the
 * incipient derivative of order one is the ordinary one, (alpha'/alpha) alpha
 * = alpha', so the equation separates and is solved componentwise. For k > 1
 * the sources define neither division nor non-integer powers in the
 * relational algebra (PLAN §6), and this refuses: GIA_E_UNSUPPORTED. */
gia_status gia_mop_couple_rel(const gia_beta beta[3], gia_rational k, double t,
                              rel_t *alpha, const char **why);

#endif /* GIA_MOP_H */
