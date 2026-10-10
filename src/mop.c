/* mop.c — the First Fundamental Equation of the MOP (include/mop.h).
 *
 * docs/requirements/srs.md §2.3, numerics.md N1, N2, N6; ADR 0019. Outputs
 * are written only on GIA_OK and `why` only on failure (IF-API-001); no
 * file-scope mutable state (NFR-REE-001); signed and complex values are never
 * clamped (NFR-NUM-006). This unit calls no harmony constructor (FR-HAR-003).
 */

#include "mop.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define EXP_MAX 709.78                  /* log(DBL_MAX), numerics N6 */
#define QUAD_TOL 1e-10                  /* N2: 1e-10 (1 + |I|) */
#define QUAD_DEPTH 50                   /* N2: at most 50 bisections */

static gia_status fail(gia_status s, const char *reason, const char **why) {
    if (why) *why = reason;
    return s;
}

static int finite_c(double complex z) {
    return isfinite(creal(z)) && isfinite(cimag(z));
}

/* ------------------------------------------------------------------ *
 * Complex log1p / expm1 (numerics N1)
 * ------------------------------------------------------------------ */

/* Kahan: log(1 + z) = z if 1 + z == 1, else log(1+z) z / ((1+z) - 1). */
static double complex clog1p_(double complex z) {
    double complex u = 1.0 + z;
    if (u == 1.0) return z;
    return clog(u) * z / (u - 1.0);
}

/* expm1(x + iy) = (expm1(x) cos y - 2 sin^2(y/2)) + i e^x sin y. */
static double complex cexpm1_(double complex z) {
    double x = creal(z), y = cimag(z), s = sin(y / 2.0);
    return (expm1(x) * cos(y) - 2.0 * s * s) + I * exp(x) * sin(y);
}

/* z^n, integer n >= 1, by binary exponentiation (numerics §1). */
static double complex ipow_c(double complex z, int n) {
    double complex r = 1.0;
    while (n > 0) {
        if (n & 1) r *= z;
        n >>= 1;
        if (n) z *= z;
    }
    return r;
}

/* alpha = S^k: integer k by multiplication, rational k only for real S > 0.
 * Overflow (Re(k log S) > log DBL_MAX) is GIA_E_RANGE. */
static gia_status power_k(double complex S, gia_rational k, double complex *out,
                          const char **why) {
    double kr = (double)k.num / (double)k.den;
    if (S == 0.0) { *out = 0.0; return GIA_OK; }
    if (kr * log(cabs(S)) > EXP_MAX)
        return fail(GIA_E_RANGE, "gia_mop_couple: alpha = S^k overflows (numerics N1, N6)", why);
    if (k.den == 1) {
        *out = ipow_c(S, k.num);
    } else {
        if (cimag(S) != 0.0 || creal(S) <= 0.0)
            return fail(GIA_E_DOMAIN, "gia_mop_couple: a non-integer k needs a real positive "
                        "S = (1/k) int beta^{1/k} (FR-MOP-002)", why);
        *out = exp(kr * log(creal(S)));
    }
    if (!finite_c(*out))
        return fail(GIA_E_RANGE, "gia_mop_couple: alpha is not finite (NFR-NUM-003)", why);
    return GIA_OK;
}

/* ------------------------------------------------------------------ *
 * N1 — affine-power beta
 * ------------------------------------------------------------------ */

static gia_status couple_affine(const gia_beta *be, gia_rational k, double t,
                                double complex *alpha, const char **why) {
    double         kr = (double)k.num / (double)k.den, p = be->p;
    double         q = (p + kr) / kr;
    double complex a = be->a, b = be->b, x, E, S;
    int            real = cimag(a) == 0.0 && cimag(b) == 0.0;

    if (!finite_c(a) || !finite_c(b) || !isfinite(p))
        return fail(GIA_E_DOMAIN, "gia_mop_couple: a, b and p must be finite", why);
    if (k.den != 1 && !(real && creal(a) > 0.0 && creal(b) >= 0.0))
        return fail(GIA_E_DOMAIN, "gia_mop_couple: a non-integer k is defined only for real "
                    "a > 0, b >= 0 (FR-MOP-002, PLAN R1)", why);
    if (a == 0.0)
        return fail(GIA_E_DOMAIN, "gia_mop_couple: a + b t vanishes at t = 0 (FR-MOP-002)", why);
    if (b != 0.0) {
        /* a + b s = 0 at s = -a/b; refused when that s is real and in [0, t]. */
        double complex s = -a / b;
        if (fabs(cimag(s)) <= 1e-15 * fmax(1.0, cabs(s)) && creal(s) >= 0.0 && creal(s) <= t)
            return fail(GIA_E_DOMAIN, "gia_mop_couple: a + b t = 0 on [0, t], where beta^{1/k} "
                        "is not defined (FR-MOP-002)", why);
        if (p + kr == 0.0)
            return fail(GIA_E_DOMAIN, "gia_mop_couple: p = -k with b != 0 makes b (p + k) = 0 "
                        "(FR-MOP-002); the integral is logarithmic, not a power", why);
    }
    if (t == 0.0) { *alpha = 0.0; return GIA_OK; }  /* [23 Eq 5.5.5]: lower limit 0 */

    /* E(x) = expm1(q log1p x) / (q x), E(0) = 1 (N1). */
    x = b * t / a;
    if (x == 0.0) {
        E = 1.0;
    } else if (real) {
        double xr = creal(x);
        E = expm1(q * log1p(xr)) / (q * xr);
    } else {
        E = cexpm1_(q * clog1p_(x)) / (q * x);
    }
    /* S = a^{p/k} t / k E(x), a^{p/k} principal. */
    /* Real a with a positive base or an integer exponent: the real pow, so
     * a real negative beta keeps an exactly real alpha. */
    if (real && (creal(a) > 0.0 || p / kr == floor(p / kr)))
        S = pow(creal(a), p / kr) * t / kr * E;
    else
        S = cexp((p / kr) * clog(a)) * t / kr * E;
    if (!finite_c(S))
        return fail(GIA_E_RANGE, "gia_mop_couple: S is not finite (NFR-NUM-003)", why);
    return power_k(S, k, alpha, why);
}

/* ------------------------------------------------------------------ *
 * N2 — adaptive Gauss-Kronrod 7-15
 * ------------------------------------------------------------------ */

static const double GK_X[8] = {
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.000000000000000000000000000000000 };
static const double GK_WK[8] = {
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714 };
static const double GK_WG[4] = {      /* Gauss 7 weights, at GK_X[1], [3], [5], [7] */
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327 };

static void gk15(gia_quad_fn g, void *ctx, double lo, double hi,
                 double complex *I_, double *err) {
    double         c = (lo + hi) / 2.0, h = (hi - lo) / 2.0;
    double complex k = 0.0, gs = 0.0, f0 = g(c, ctx);
    int            j;
    k  = GK_WK[7] * f0;
    gs = GK_WG[3] * f0;
    for (j = 0; j < 7; j++) {
        double complex fp = g(c + h * GK_X[j], ctx), fm = g(c - h * GK_X[j], ctx);
        k += GK_WK[j] * (fp + fm);
        if (j % 2 == 1) gs += GK_WG[j / 2] * (fp + fm);
    }
    *I_  = k * h;
    *err = cabs((k - gs) * h);
}

/* Recursive bisection: returns 0 when the depth is exhausted unmet. */
static int gk_adapt(gia_quad_fn g, void *ctx, double lo, double hi, double tol_abs,
                    int depth, double complex *I_, double *err) {
    double complex I1, I2;
    double         e1, e2, mid;
    gk15(g, ctx, lo, hi, I_, err);
    if (!finite_c(*I_) || !isfinite(*err)) return 0;
    if (*err <= tol_abs) return 1;
    if (depth <= 0) return 0;
    mid = (lo + hi) / 2.0;
    if (!gk_adapt(g, ctx, lo, mid, tol_abs / 2.0, depth - 1, &I1, &e1)) return 0;
    if (!gk_adapt(g, ctx, mid, hi, tol_abs / 2.0, depth - 1, &I2, &e2)) return 0;
    *I_ = I1 + I2; *err = e1 + e2;
    return 1;
}

gia_status gia_quad_gk15(gia_quad_fn g, void *ctx, double lo, double hi, double tol,
                         int max_depth, double complex *integral, double *err,
                         const char **why) {
    double complex I0, I_;
    double         e0, e;
    if (!g || !integral || !err)
        return fail(GIA_E_ARG, "gia_quad_gk15: NULL argument", why);
    if (!isfinite(lo) || !isfinite(hi) || !(tol > 0.0) || max_depth < 0)
        return fail(GIA_E_ARG, "gia_quad_gk15: bad interval, tolerance or depth", why);
    if (lo == hi) { *integral = 0.0; *err = 0.0; return GIA_OK; }
    /* The tolerance is relative to the integral's size: estimate it first. */
    gk15(g, ctx, lo, hi, &I0, &e0);
    if (!gk_adapt(g, ctx, lo, hi, tol * (1.0 + cabs(I0)), max_depth, &I_, &e))
        return fail(GIA_E_CONVERGENCE, "gia_quad_gk15: tolerance not met within the "
                    "subdivision limit (numerics N2, NFR-NUM-005)", why);
    *integral = I_; *err = e;
    return GIA_OK;
}

/* ------------------------------------------------------------------ *
 * N2 — sampled beta
 * ------------------------------------------------------------------ */

typedef struct {
    double complex b0, b1;      /* segment end values */
    double         t0, t1;
    double         theta0;      /* continued argument at t0 */
    double         invk;        /* 1/k */
} seg_ctx;

/* g(t) = |beta|^{1/k} e^{i theta / k}, theta = theta0 + arg(beta(t)/beta(t0)). */
static double complex seg_g(double t, void *c) {
    const seg_ctx *s = (const seg_ctx *)c;
    double         w = (t - s->t0) / (s->t1 - s->t0);
    double complex bt = s->b0 + w * (s->b1 - s->b0);
    double         th = s->theta0 + carg(bt / s->b0);
    return exp(s->invk * log(cabs(bt))) * cexp(I * th * s->invk);
}

/* Does the segment from u to v pass through 0? */
static int segment_hits_zero(double complex u, double complex v) {
    double complex d = v - u;
    double         dd = creal(d) * creal(d) + cimag(d) * cimag(d), s, cr;
    if (u == 0.0 || v == 0.0) return 1;
    if (dd == 0.0) return 0;
    s  = -(creal(u) * creal(d) + cimag(u) * cimag(d)) / dd;      /* closest point */
    cr = creal(u) * cimag(d) - cimag(u) * creal(d);              /* collinear with 0? */
    return s >= 0.0 && s <= 1.0 && fabs(cr) <= 1e-15 * (cabs(u) + cabs(v)) * cabs(d);
}

static gia_status couple_samples(const gia_beta *be, gia_rational k, double t,
                                 double complex *alpha, const char **why) {
    double         kr = (double)k.num / (double)k.den, theta;
    double complex sum = 0.0;
    int            i;

    if (be->n < 2 || !be->t || !be->v)
        return fail(GIA_E_ARG, "gia_mop_couple: samples need n >= 2 and both arrays", why);
    if (be->t[0] != 0.0)
        return fail(GIA_E_DOMAIN, "gia_mop_couple: samples must start at t = 0", why);
    for (i = 0; i < be->n; i++) {
        if (!isfinite(be->t[i]) || !finite_c(be->v[i]))
            return fail(GIA_E_DOMAIN, "gia_mop_couple: a sample is not finite", why);
        if (i > 0 && !(be->t[i] > be->t[i - 1]))
            return fail(GIA_E_DOMAIN, "gia_mop_couple: sample times must increase", why);
        if (i > 0 && segment_hits_zero(be->v[i - 1], be->v[i]))
            return fail(GIA_E_DOMAIN, "gia_mop_couple: beta passes through 0 between samples, "
                        "where beta^{1/k} has no continued branch (numerics N2)", why);
    }
    if (t > be->t[be->n - 1])
        return fail(GIA_E_DOMAIN, "gia_mop_couple: t is past the last sample (FR-MOP-004)", why);
    if (t == 0.0) { *alpha = 0.0; return GIA_OK; }

    theta = carg(be->v[0]);                     /* principal at t = 0 */
    for (i = 0; i + 1 < be->n && be->t[i] < t; i++) {
        seg_ctx        sc;
        double complex Iseg;
        double         err, hi = be->t[i + 1] < t ? be->t[i + 1] : t;
        gia_status     st;
        sc.b0 = be->v[i]; sc.b1 = be->v[i + 1];
        sc.t0 = be->t[i]; sc.t1 = be->t[i + 1];
        sc.theta0 = theta; sc.invk = 1.0 / kr;
        st = gia_quad_gk15(seg_g, &sc, be->t[i], hi, QUAD_TOL, QUAD_DEPTH, &Iseg, &err, why);
        if (st != GIA_OK) return st;
        sum += Iseg;
        theta += carg(be->v[i + 1] / be->v[i]);   /* unwrap to the next sample */
    }
    return power_k(sum / kr, k, alpha, why);
}

/* ------------------------------------------------------------------ *
 * Public
 * ------------------------------------------------------------------ */

gia_status gia_mop_couple(const gia_beta *beta, gia_rational k, double t,
                          double complex *alpha, const char **why) {
    double complex al;
    gia_status     st;
    if (!beta || !alpha) return fail(GIA_E_ARG, "gia_mop_couple: NULL argument", why);
    if (k.den <= 0) return fail(GIA_E_ARG, "gia_mop_couple: k's denominator must be >= 1", why);
    if (k.num <= 0)
        return fail(GIA_E_DOMAIN, "gia_mop_couple: the cardinality k must be > 0 (FR-MOP-002)", why);
    if (!(t >= 0.0) || !isfinite(t)) return fail(GIA_E_ARG, "gia_mop_couple: t must be >= 0", why);
    switch (beta->kind) {
    case GIA_BETA_AFFINE:  st = couple_affine(beta, k, t, &al, why);  break;
    case GIA_BETA_SAMPLES: st = couple_samples(beta, k, t, &al, why); break;
    default: return fail(GIA_E_ARG, "gia_mop_couple: the couple is unrelated (GIA_BETA_NONE)", why);
    }
    if (st == GIA_OK) *alpha = al;
    return st;
}

gia_status gia_mop_solve(int N, const gia_beta *beta, gia_rational k, double t,
                         gia_matrioska *out, const char **why) {
    double complex *a;
    unsigned char  *rel;
    int             i, j;
    gia_status      st;

    if (N < 1 || !beta || !out) return fail(GIA_E_ARG, "gia_mop_solve: bad argument", why);
    a   = (double complex *)calloc((size_t)N * (size_t)N, sizeof(double complex));
    rel = (unsigned char *)calloc((size_t)N * (size_t)N, 1);
    if (!a || !rel) {
        free(a); free(rel);
        return fail(GIA_E_NOMEM, "gia_mop_solve: out of memory", why);
    }
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
            const gia_beta *b = &beta[(size_t)i * (size_t)N + (size_t)j];
            if (i == j || b->kind == GIA_BETA_NONE) continue;
            st = gia_mop_couple(b, k, t, &a[(size_t)i * (size_t)N + (size_t)j], why);
            if (st != GIA_OK) { free(a); free(rel); return st; }
            rel[(size_t)i * (size_t)N + (size_t)j] = 1;
        }
    out->N = N; out->a = a; out->related = rel;
    return GIA_OK;
}

void gia_matrioska_free(gia_matrioska *m) {
    if (!m) return;
    free(m->a); free(m->related);
    m->a = NULL; m->related = NULL; m->N = 0;
}

/* FR-MOP-007. */
gia_status gia_mop_couple_rel(const gia_beta beta[3], gia_rational k, double t,
                              rel_t *alpha, const char **why) {
    double complex c[3];
    int            n;
    gia_status     st;
    if (!beta || !alpha) return fail(GIA_E_ARG, "gia_mop_couple_rel: NULL argument", why);
    if (!(k.num == 1 && k.den == 1))
        return fail(GIA_E_UNSUPPORTED, "a relational-valued First Equation is defined only for "
                    "k = 1: the sources define no division or non-integer power in the "
                    "relational algebra (FR-MOP-007, PLAN §6)", why);
    for (n = 0; n < 3; n++) {
        if ((st = gia_mop_couple(&beta[n], k, t, &c[n], why)) != GIA_OK) return st;
        if (cimag(c[n]) != 0.0)
            return fail(GIA_E_DOMAIN, "gia_mop_couple_rel: each component's boundary condition "
                        "must be real", why);
    }
    alpha->i = creal(c[0]); alpha->j = creal(c[1]); alpha->k = creal(c[2]);
    return GIA_OK;
}

/* FR-MOP-006 — [23 Eq 7.1-7.5]; numerics N8. */
gia_status gia_eqs(const gia_eqs_params *p, rel_t ref, int l, double out[3],
                   const char **why) {
    const double pi = 3.14159265358979323846;
    double       r2psi, B, C, E[3], expo, res[3];
    rel_t        root, br;
    int          i;

    if (!p || !out) return fail(GIA_E_ARG, "gia_eqs: NULL argument", why);
    if (p->N < 3 || l < 1 || l > p->N - 1)
        return fail(GIA_E_ARG, "gia_eqs: needs N >= 3 and l = 1..N-1", why);
    if (p->eps[1] != p->eps[2])
        return fail(GIA_E_DOMAIN, "[23 Eq 7.4] writes one angle sqrt2 psi for both the j and the "
                    "k spinor, which holds only for eps_2 = eps_3 (PLAN X11)", why);
    for (i = 0; i < 3; i++)
        if (!isfinite(p->psi1[i]) || !isfinite(p->eps[i]))
            return fail(GIA_E_DOMAIN, "gia_eqs: parameters must be finite", why);
    if (!isfinite(p->psi2) || !isfinite(p->A) || !isfinite(ref.i) || !isfinite(ref.j) ||
        !isfinite(ref.k))
        return fail(GIA_E_DOMAIN, "gia_eqs: parameters must be finite", why);

    r2psi = p->psi2 * (p->eps[1] + 2.0 * pi * l) / (double)(p->N - 1);     /* (7.5) */
    B = cos(r2psi);
    C = sin(r2psi) / sqrt(2.0);                                             /* (7.4) */
    for (i = 0; i < 3; i++) E[i] = (p->eps[i] + 4.0 * pi * l) / (double)(p->N - 1);
    root.i = B; root.j = C; root.k = C;
    br = rel_mul(root, ref);            /* the three brackets, PLAN R2 */

    expo = p->psi1[0] * E[0] * br.i;
    if (expo > 709.78)
        return fail(GIA_E_RANGE, "gia_eqs: rho = A e^S overflows (numerics N6)", why);
    res[0] = p->A * exp(expo);                                              /* (7.1) */
    res[1] = p->psi1[1] * E[1] * br.j;                                      /* (7.2) */
    res[2] = p->psi1[2] * E[2] * br.k;                                      /* (7.3) */
    for (i = 0; i < 3; i++)
        if (!isfinite(res[i]))
            return fail(GIA_E_RANGE, "gia_eqs: a coordinate is not finite", why);
    for (i = 0; i < 3; i++) out[i] = res[i];
    return GIA_OK;
}

/* FR-MOP-005. */
gia_status gia_mop_second(double complex alpha12_0, double c1, double c2, int N, double t,
                          gia_second *out, gia_matrioska *r, const char **why) {
    const double    two_pi = 6.283185307179586476925286766559;
    double complex  w, A, eB11, *a = NULL;
    unsigned char  *rel = NULL;
    int             j;

    if (!out || N < 2 || !finite_c(alpha12_0) || !isfinite(c1) || !isfinite(c2) ||
        !(t >= 0.0) || !isfinite(t))
        return fail(GIA_E_ARG, "gia_mop_second: N >= 2, t >= 0 and finite arguments required", why);
    /* c1 + c2 s is affine in s, so it is positive on [0, t] iff it is at both ends. */
    if (!(c1 > 0.0) || !(c1 + c2 * t > 0.0))
        return fail(GIA_E_DOMAIN, "gia_mop_second: c1 + c2 s <= 0 on [0, t]; ln is undefined "
                    "(FR-MOP-005)", why);
    w    = cexp(I * (two_pi / (double)(N - 1)));
    A    = alpha12_0 * w + log(c1 + c2 * t);
    eB11 = 1.0 + (cexp(2.0 * A) - 1.0) / 2.0;
    if (!finite_c(A) || !finite_c(eB11))
        return fail(GIA_E_RANGE, "gia_mop_second: e^{2A} overflows (N6)", why);
    if (r) {
        a   = (double complex *)calloc((size_t)N * (size_t)N, sizeof(double complex));
        rel = (unsigned char *)calloc((size_t)N * (size_t)N, 1);
        if (!a || !rel) {
            free(a); free(rel);
            return fail(GIA_E_NOMEM, "gia_mop_second: out of memory", why);
        }
        /* Row 1, couples 1j for j = 2..N: 0-based index j, exponent j - 1. */
        for (j = 1; j < N; j++) {
            a[j]   = eB11 * cexp(I * (two_pi * (double)(j - 1) / (double)(N - 1)));
            rel[j] = 1;
        }
        r->N = N; r->a = a; r->related = rel;
    }
    out->A       = A;
    out->B[0][0] = A;  out->B[0][1] = -A;
    out->B[1][0] = -A; out->B[1][1] = A;
    return GIA_OK;
}

/* ------------------------------------------------------------------ *
 * The harmony detector (FR-HAR-001, FR-HAR-002)
 * ------------------------------------------------------------------ */

#define GIA_TWO_PI   6.283185307179586476925286766559
#define HARMONY_TOL  1e-9                  /* vv-plan.md §6 */
#define HARMONY_DELTA 1e-3

const char *gia_verdict_str(gia_verdict v) {
    switch (v) {
    case GIA_H_IMPOSED:     return "imposed";
    case GIA_H_TRANSPORTED: return "transported";
    case GIA_H_PRESENT:     return "present";
    case GIA_H_ABSENT:      return "absent";
    default:                return "unknown verdict";
    }
}

/* FR-HAR-001. */
gia_status gia_harmony_residual(const gia_matrioska *alpha, double *R_H, const char **why) {
    double complex a12;
    double         worst = 0.0;
    int            N, j;

    if (!alpha || !R_H || !alpha->a || !alpha->related || alpha->N < 3)
        return fail(GIA_E_ARG, "gia_harmony_residual: N >= 3 required ([23 Eq 5.6.5])", why);
    N = alpha->N;
    for (j = 1; j < N; j++)
        if (!alpha->related[j])
            return fail(GIA_E_DOMAIN, "gia_harmony_residual: a couple of row 1 is unrelated, so "
                        "[23 Eq 5.6.5] has no ratio for it", why);
    a12 = alpha->a[1];
    if (!finite_c(a12) || !(cabs(a12) >= DBL_MIN))
        return fail(GIA_E_DOMAIN, "gia_harmony_residual: alpha_12 = 0, so the ratios of "
                    "[23 Eq 5.6.5] are undefined (FR-HAR-001)", why);
    for (j = 2; j < N; j++) {
        const double complex ratio = alpha->a[j] / a12;
        const double complex root  = cexp(I * (GIA_TWO_PI * (double)(j - 1) / (double)(N - 1)));
        double               d;
        if (!finite_c(ratio))
            return fail(GIA_E_RANGE, "gia_harmony_residual: a ratio alpha_1j/alpha_12 is not "
                        "finite (N6)", why);
        d = cabs(ratio - root);
        if (d > worst) worst = d;
    }
    *R_H = worst;
    return GIA_OK;
}

gia_status gia_harmony_observed(const gia_matrioska *alpha, gia_verdict *v, const char **why) {
    double     R;
    gia_status st;
    if (!v) return fail(GIA_E_ARG, "gia_harmony_observed: NULL verdict", why);
    if ((st = gia_harmony_residual(alpha, &R, why)) != GIA_OK) return st;
    *v = R <= HARMONY_TOL ? GIA_H_PRESENT : GIA_H_ABSENT;
    return GIA_OK;
}

/* vv-plan.md §6's 64-bit LCG, mapped to [-1, 1]. */
static double harmony_lcg(uint64_t *x) {
    *x = 6364136223846793005ULL * *x + 1442695040888963407ULL;
    return 2.0 * ((double)(*x >> 11) * (1.0 / 9007199254740992.0)) - 1.0;
}

/* Input m (0: harmonic; 1..8: perturbed) and its own Matrioska at t = 1. */
static void harmony_input(int N, int m, gia_beta *in, gia_matrioska *at1) {
    uint64_t x = (uint64_t)m;
    int      j;
    memset(in, 0, (size_t)N * (size_t)N * sizeof(gia_beta));
    for (j = 1; j < N; j++) {
        const double complex w = cexp(I * (GIA_TWO_PI * (double)(j - 1) / (double)(N - 1)));
        double complex       f = 1.0;
        if (m > 0) {
            const double xr = harmony_lcg(&x), xi = harmony_lcg(&x);
            f = 1.0 + HARMONY_DELTA * (xr + xi * I);
        }
        in[j].kind = GIA_BETA_AFFINE;
        in[j].a    = w * f;
        in[j].b    = 0.25 * w * f;
        in[j].p    = 1.0;
        at1->a[j]       = in[j].a + in[j].b;      /* (a + b t)^1 at t = 1 */
        at1->related[j] = 1;
    }
}

/* FR-HAR-002, by vv-plan.md §6. */
gia_status gia_harmony_verdict(gia_construction c, void *ctx, int N, gia_verdict *v,
                               const char **why) {
    gia_beta      *in;
    gia_matrioska  at1;
    int            m, imposed = 1, transported = 1, some_in = 0;
    gia_status     st = GIA_OK;

    if (!c || !v || N < 3)
        return fail(GIA_E_ARG, "gia_harmony_verdict: a construction, and N >= 3", why);
    in          = (gia_beta *)malloc((size_t)N * (size_t)N * sizeof(gia_beta));
    at1.N       = N;
    at1.a       = (double complex *)calloc((size_t)N * (size_t)N, sizeof(double complex));
    at1.related = (unsigned char *)calloc((size_t)N * (size_t)N, 1);
    if (!in || !at1.a || !at1.related) {
        free(in); gia_matrioska_free(&at1);
        return fail(GIA_E_NOMEM, "gia_harmony_verdict: out of memory", why);
    }
    for (m = 0; m <= 8 && st == GIA_OK; m++) {
        gia_matrioska out;
        double        r_in = 0.0, r_out = 0.0;
        int           has_out;
        harmony_input(N, m, in, &at1);
        if ((st = gia_harmony_residual(&at1, &r_in, why)) != GIA_OK) break;
        if (r_in > HARMONY_TOL) some_in = 1;
        memset(&out, 0, sizeof out);
        if ((st = c(in, N, &out, ctx)) != GIA_OK) break;
        has_out = out.N == N && gia_harmony_residual(&out, &r_out, NULL) == GIA_OK;
        gia_matrioska_free(&out);
        if (!has_out) { imposed = transported = 0; continue; }
        if (!(r_out <= HARMONY_TOL)) imposed = 0;
        if (!(fabs(r_out - r_in) <= HARMONY_TOL * (r_in > 1.0 ? r_in : 1.0))) transported = 0;
    }
    free(in);
    gia_matrioska_free(&at1);
    if (st != GIA_OK) return st;
    *v = imposed ? GIA_H_IMPOSED : (transported && some_in) ? GIA_H_TRANSPORTED : GIA_H_ABSENT;
    return GIA_OK;
}

/* beta_12's value at t, from a related affine-power couple. */
static int beta12_at(const gia_beta *in, double t, double complex *v) {
    const gia_beta *b = &in[1];
    if (b->kind != GIA_BETA_AFFINE) return 0;
    *v = cpow(b->a + b->b * t, b->p);
    return finite_c(*v);
}

gia_status gia_construct_first(const gia_beta *in, int N, gia_matrioska *out, void *ctx) {
    static const gia_rational k1 = {1, 1};
    if (!in || !out || N < 3) return GIA_E_ARG;
    return gia_mop_solve(N, in, ctx ? *(const gia_rational *)ctx : k1, 1.0, out, NULL);
}

gia_status gia_construct_second(const gia_beta *in, int N, gia_matrioska *out, void *ctx) {
    gia_second     sec;
    double complex a0;
    (void)ctx;
    if (!in || !out || N < 3 || !beta12_at(in, 0.0, &a0)) return GIA_E_ARG;
    return gia_mop_second(a0, 1.0, 0.5, N, 1.0, &sec, out, NULL);
}

gia_status gia_construct_eqs(const gia_beta *in, int N, gia_matrioska *out, void *ctx) {
    gia_eqs_params p;
    const rel_t    ref = {1.0, 0.0, 0.0};
    double complex b12, *a;
    unsigned char *rel;
    gia_status     st;
    int            l;
    (void)ctx;
    if (!in || !out || N < 3 || !beta12_at(in, 1.0, &b12)) return GIA_E_ARG;
    p.psi1[0] = p.psi1[1] = p.psi1[2] = 1.0;
    p.psi2    = 1.0;
    p.eps[0]  = p.eps[1] = p.eps[2] = -GIA_TWO_PI;
    p.A       = 1.0;
    p.N       = N;
    a   = (double complex *)calloc((size_t)N * (size_t)N, sizeof(double complex));
    rel = (unsigned char *)calloc((size_t)N * (size_t)N, 1);
    if (!a || !rel) { free(a); free(rel); return GIA_E_NOMEM; }
    for (l = 1; l < N; l++) {
        /* With ref = {1, 0, 0} the brackets of [23 Eq 7.1-7.3] are B and C, so
         * rho = e^{E_l1 B} and phi = E_l2 C (psi1 = 1, A = 1). */
        const double E = (p.eps[0] + 2.0 * GIA_TWO_PI * (double)l) / (double)(N - 1);
        double       c[3], B, C;
        if ((st = gia_eqs(&p, ref, l, c, NULL)) != GIA_OK) { free(a); free(rel); return st; }
        B = log(c[0]) / E;
        C = c[1] / E;
        a[l]   = b12 * (B + I * sqrt(2.0) * C);
        rel[l] = 1;
    }
    out->N = N; out->a = a; out->related = rel;
    return GIA_OK;
}
