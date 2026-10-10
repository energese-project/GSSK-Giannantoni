/* idc.c — single-variable Incipient Differential Calculus (include/idc.h).
 *
 * Every function here implements one requirement of docs/requirements/srs.md
 * §2.1 by the method numerics.md names for it, and returns a gia_status
 * (IF-API-001): outputs are written only on GIA_OK, `why` only on failure,
 * and nothing prints, exits or aborts. No file-scope mutable state
 * (NFR-REE-001), and no clamping of signed or complex values (NFR-NUM-006).
 */

#include "idc.h"

#include <math.h>
#include <stdlib.h>

/* gia_status_str lives with the first unit that returns a gia_status; mop.c
 * and relational.c link idc.o for it. */
const char *gia_status_str(gia_status s) {
    switch (s) {
    case GIA_OK:            return "ok";
    case GIA_E_ARG:         return "invalid argument";
    case GIA_E_DOMAIN:      return "outside the function's domain";
    case GIA_E_RANGE:       return "result out of range";
    case GIA_E_CONVERGENCE: return "did not converge";
    case GIA_E_UNSUPPORTED: return "not defined by the available sources";
    case GIA_E_LIMIT:       return "a fixed limit was exceeded";
    case GIA_E_NOMEM:       return "out of memory";
    }
    return "unknown status";
}

static gia_status fail(gia_status s, const char *reason, const char **why) {
    if (why) *why = reason;
    return s;
}

static int finite_c(double complex z) {
    return isfinite(creal(z)) && isfinite(cimag(z));
}

/* z^n for integer n >= 0 by binary exponentiation: multiplication only, so
 * the sign and branch are exact (numerics.md §1, "Integer powers"). */
static double complex ipow_c(double complex z, int n) {
    double complex r = 1.0;
    while (n > 0) {
        if (n & 1) r *= z;
        n >>= 1;
        if (n) z *= z;
    }
    return r;
}

/* FR-IDC-002 — [02 Eq 14.9.5], [23 Eq 5.5.2]: (d~/dt)^n f = (f'/f)^n f. */
gia_status gia_idc_of(double complex f, double complex df, int n,
                      double complex *out, const char **why) {
    double complex r;
    if (!out)  return fail(GIA_E_ARG, "gia_idc_of: out is NULL", why);
    if (n < 0) return fail(GIA_E_ARG, "gia_idc_of: the cardinality n must be >= 0", why);
    if (!finite_c(f) || !finite_c(df))
        return fail(GIA_E_DOMAIN, "gia_idc_of: f and f' must be finite", why);
    if (f == 0.0)
        return fail(GIA_E_DOMAIN,
                    "gia_idc_of: (f'/f)^n f needs f != 0 ([02 Eq 14.9.5]); "
                    "the incipient derivative is undefined where f vanishes", why);
    if (n == 0) { *out = f; return GIA_OK; }
    r = ipow_c(df / f, n) * f;
    if (!finite_c(r))
        return fail(GIA_E_RANGE, "gia_idc_of: (f'/f)^n f overflows (NFR-NUM-003)", why);
    *out = r;
    return GIA_OK;
}

/* ================================================================== *
 * FR-IDC-006 — the second-order incipient LDE (numerics.md N3)
 * ================================================================== */

#define N3_TOL   1e-12      /* discriminant, step agreement, IC singularity */
#define EXP_MAX  709.78     /* log(DBL_MAX): exp beyond this overflows (N6) */
#define N3_MAX_NODES 4000000

/* 5-point Gauss-Legendre on [-1, 1]. */
static const double GL5_X[5] = {
    -0.906179845938663992797626878299, -0.538469310105683091036314420700, 0.0,
     0.538469310105683091036314420700,  0.906179845938663992797626878299 };
static const double GL5_W[5] = {
     0.236926885056189087514264040720,  0.478628670499366468041291514836,
     0.568888888888888888888888888889,
     0.478628670499366468041291514836,  0.236926885056189087514264040720 };

struct gia_lde2_sol {
    gia_cfn         a1, a0;
    void           *ctx;
    double          t_max;
    int             one_family;     /* double root everywhere: c e^{int r}  */
    double complex  c[2];
    int             n, cap;         /* accepted step ends, t[0] = 0         */
    double         *t;
    double complex *r;              /* labelled roots at t[k]: r[2k], r[2k+1] */
    double complex *Ic;              /* int_0^{t[k]} r_i: Ic[2k], Ic[2k+1]       */
};

/* The coefficients and roots at t. The stable form of N3:
 * s = csqrt(a1^2 - 4 a0), q = a1 + sigma s with |q| largest, r1 = -q/2,
 * r2 = a0/r1. `dbl` reports a double root by N3's discriminant test. */
static gia_status roots_at(const gia_lde2_sol *sl, double t, double complex r[2],
                           double complex *a1o, double complex *a0o, int *dbl) {
    double complex a1 = sl->a1(t, sl->ctx), a0 = sl->a0(t, sl->ctx), disc, sq, q;
    if (!finite_c(a1) || !finite_c(a0)) return GIA_E_DOMAIN;
    disc = a1 * a1 - 4.0 * a0;
    if (dbl) *dbl = cabs(disc) <= N3_TOL * (cabs(a1) * cabs(a1) + cabs(a0));
    if (sl->one_family) {
        /* The double root itself, exactly -a1/2, not two rounded copies. */
        r[0] = r[1] = -a1 / 2.0;
    } else {
        sq = csqrt(disc);
        q  = cabs(a1 + sq) >= cabs(a1 - sq) ? a1 + sq : a1 - sq;
        if (q == 0.0) { r[0] = r[1] = 0.0; }
        else          { r[0] = -q / 2.0; r[1] = a0 / r[0]; }
    }
    if (a1o) *a1o = a1;
    if (a0o) *a0o = a0;
    return GIA_OK;
}

/* Order q so that q[i] continues prev[i]: the assignment with the smaller
 * total move (N3, "labelling by continuity"). */
static void label(const double complex prev[2], double complex q[2]) {
    if (cabs(q[0] - prev[1]) + cabs(q[1] - prev[0]) <
        cabs(q[0] - prev[0]) + cabs(q[1] - prev[1])) {
        double complex tmp = q[0]; q[0] = q[1]; q[1] = tmp;
    }
}

/* GL5 over [ta, tb], labelling each node from the previous one, starting from
 * the labelled roots `prev` at ta. Writes the integrals, the labelled roots at
 * tb, the largest move of either root away from `prev`, and whether every node
 * was a double root. */
static gia_status gl5_step(const gia_lde2_sol *sl, double ta, double tb,
                           const double complex prev[2], double complex Ic[2],
                           double complex end[2], double *move, int *all_dbl) {
    double         half = (tb - ta) / 2.0, mid = (ta + tb) / 2.0;
    double complex last[2], q[2];
    int            k, dbl;
    gia_status     st;

    Ic[0] = Ic[1] = 0.0;
    last[0] = prev[0]; last[1] = prev[1];
    *move = 0.0; *all_dbl = 1;
    for (k = 0; k < 5; k++) {
        st = roots_at(sl, mid + half * GL5_X[k], q, NULL, NULL, &dbl);
        if (st != GIA_OK) return st;
        label(last, q);
        if (!dbl) *all_dbl = 0;
        Ic[0] += GL5_W[k] * q[0];
        Ic[1] += GL5_W[k] * q[1];
        if (cabs(q[0] - prev[0]) > *move) *move = cabs(q[0] - prev[0]);
        if (cabs(q[1] - prev[1]) > *move) *move = cabs(q[1] - prev[1]);
        last[0] = q[0]; last[1] = q[1];
    }
    Ic[0] *= half; Ic[1] *= half;
    st = roots_at(sl, tb, q, NULL, NULL, &dbl);
    if (st != GIA_OK) return st;
    label(last, q);
    if (!dbl) *all_dbl = 0;
    if (cabs(q[0] - prev[0]) > *move) *move = cabs(q[0] - prev[0]);
    if (cabs(q[1] - prev[1]) > *move) *move = cabs(q[1] - prev[1]);
    end[0] = q[0]; end[1] = q[1];
    return GIA_OK;
}

static int push_node(gia_lde2_sol *sl, double t, const double complex r[2],
                     const double complex Ic[2]) {
    if (sl->n == sl->cap) {
        int cap = sl->cap ? 2 * sl->cap : 64;
        double *nt; double complex *nr, *nI;
        if (cap > N3_MAX_NODES) return 0;
        nt = (double *)realloc(sl->t, (size_t)cap * sizeof(double));
        if (!nt) return 0;
        sl->t = nt;
        nr = (double complex *)realloc(sl->r, (size_t)cap * 2 * sizeof(double complex));
        if (!nr) return 0;
        sl->r = nr;
        nI = (double complex *)realloc(sl->Ic, (size_t)cap * 2 * sizeof(double complex));
        if (!nI) return 0;
        sl->Ic = nI;
        sl->cap = cap;
    }
    sl->t[sl->n] = t;
    sl->r[2 * sl->n] = r[0]; sl->r[2 * sl->n + 1] = r[1];
    sl->Ic[2 * sl->n] = Ic[0]; sl->Ic[2 * sl->n + 1] = Ic[1];
    sl->n++;
    return 1;
}

void gia_lde2_free(gia_lde2_sol *sol) {
    if (!sol) return;
    free(sol->t); free(sol->r); free(sol->Ic);
    free(sol);
}

/* The adaptive sweep of N3 from 0 to t_max. */
static gia_status sweep(gia_lde2_sol *sl, double *t_fail, const char **why) {
    const double   hmin = N3_TOL * (sl->t_max > 1.0 ? sl->t_max : 1.0);
    double         ta = 0.0, h = sl->t_max;
    double complex cur[2], cum[2] = { 0.0, 0.0 };
    gia_status     st;

    cur[0] = sl->r[0]; cur[1] = sl->r[1];
    while (ta < sl->t_max) {
        double         tb = ta + h, tm, move, move2, sep;
        double complex If[2], Ia[2], Ib[2], ef[2], ea[2], eb[2];
        int            dblf, dbla, dblb, agree, i;

        if (tb > sl->t_max || sl->t_max - tb < hmin) tb = sl->t_max;
        tm = (ta + tb) / 2.0;
        if ((st = gl5_step(sl, ta, tb, cur, If, ef, &move, &dblf)) != GIA_OK ||
            (st = gl5_step(sl, ta, tm, cur, Ia, ea, &move2, &dbla)) != GIA_OK)
            return fail(st, "gia_lde2_solve: a coefficient is not finite on [0, t_max]", why);
        if (move2 > move) move = move2;
        if ((st = gl5_step(sl, tm, tb, ea, Ib, eb, &move2, &dblb)) != GIA_OK)
            return fail(st, "gia_lde2_solve: a coefficient is not finite on [0, t_max]", why);

        if (sl->one_family && !(dblf && dbla && dblb))
            return fail(GIA_E_DOMAIN,
                        "gia_lde2_solve: the roots coincide at t = 0 but separate after it, "
                        "so f(0) and f~'(0) cannot fix the two constants (PLAN R13, numerics N3)",
                        why);

        agree = 1;
        for (i = 0; i < 2; i++)
            if (cabs(If[i] - (Ia[i] + Ib[i])) > N3_TOL * (1.0 + cabs(Ia[i] + Ib[i]))) agree = 0;
        sep = cabs(cur[0] - cur[1]);
        if (agree && (sl->one_family || move <= 0.5 * sep)) {
            double complex Inew[2];
            Inew[0] = cum[0] + Ia[0] + Ib[0];
            Inew[1] = cum[1] + Ia[1] + Ib[1];
            if (!push_node(sl, tb, eb, Inew))
                return fail(GIA_E_LIMIT, "gia_lde2_solve: the sweep needed more steps than it may store", why);
            cum[0] = Inew[0]; cum[1] = Inew[1];
            cur[0] = eb[0];   cur[1] = eb[1];
            ta = tb;
            h *= 2.0;
            continue;
        }
        h /= 2.0;
        if (h < hmin) {
            if (t_fail) *t_fail = ta;
            return fail(GIA_E_CONVERGENCE,
                        "gia_lde2_solve: the roots of r^2 + a1 r + a0 collide at an isolated "
                        "time (reported in t_fail), where they cannot be labelled by "
                        "continuity (numerics N3)", why);
        }
    }
    return GIA_OK;
}

gia_status gia_lde2_solve(gia_cfn a1, gia_cfn a0, void *ctx,
                          double complex f0, double complex f1, double t_max,
                          gia_lde2_sol **sol, double *t_fail, const char **why) {
    gia_lde2_sol  *sl;
    double complex r0[2], zero[2] = { 0.0, 0.0 };
    int            dbl0;
    gia_status     st;

    if (!a1 || !a0 || !sol)
        return fail(GIA_E_ARG, "gia_lde2_solve: a1, a0 and sol must not be NULL", why);
    if (!(t_max > 0.0) || !isfinite(t_max))
        return fail(GIA_E_ARG, "gia_lde2_solve: t_max must be finite and > 0", why);
    if (!finite_c(f0) || !finite_c(f1))
        return fail(GIA_E_DOMAIN, "gia_lde2_solve: initial conditions must be finite", why);

    sl = (gia_lde2_sol *)calloc(1, sizeof(*sl));
    if (!sl) return fail(GIA_E_NOMEM, "gia_lde2_solve: out of memory", why);
    sl->a1 = a1; sl->a0 = a0; sl->ctx = ctx; sl->t_max = t_max;

    if ((st = roots_at(sl, 0.0, r0, NULL, NULL, &dbl0)) != GIA_OK) {
        gia_lde2_free(sl);
        return fail(st, "gia_lde2_solve: a coefficient is not finite at t = 0", why);
    }
    if (dbl0) {
        /* A double root at 0. If it stays double the solution is one family
         * (FR-IDC-006); the sweep refuses if it does not. */
        sl->one_family = 1;
        (void)roots_at(sl, 0.0, r0, NULL, NULL, NULL);
        if (cabs(f1 - r0[0] * f0) >
            N3_TOL * fmax(1.0, fmax(cabs(f1), cabs(r0[0] * f0)))) {
            gia_lde2_free(sl);
            return fail(GIA_E_DOMAIN,
                        "gia_lde2_solve: a double root admits only f = c e^{int r}, so "
                        "f~'(0) must equal r(0) f(0); [06 Eq 3.7]'s second solution solves "
                        "the equation under no reading (PLAN X12)", why);
        }
        sl->c[0] = f0; sl->c[1] = 0.0;
    } else {
        double complex det = r0[1] - r0[0];
        sl->c[0] = (f0 * r0[1] - f1) / det;
        sl->c[1] = (f1 - f0 * r0[0]) / det;
    }
    if (!push_node(sl, 0.0, r0, zero) || (st = sweep(sl, t_fail, why)) != GIA_OK) {
        if (sl->n == 0) st = fail(GIA_E_NOMEM, "gia_lde2_solve: out of memory", why);
        gia_lde2_free(sl);
        return st;
    }
    *sol = sl;
    return GIA_OK;
}

/* Labelled roots and integrals at t, by GL5 from the last stored node. */
static gia_status state_at(const gia_lde2_sol *sl, double t, double complex r[2],
                           double complex Ic[2]) {
    int            lo = 0, hi = sl->n - 1, dbl;
    double         move;
    double complex dI[2];
    gia_status     st;

    while (lo < hi) {                       /* largest k with t[k] <= t */
        int mid = (lo + hi + 1) / 2;
        if (sl->t[mid] <= t) lo = mid; else hi = mid - 1;
    }
    if (sl->t[lo] == t) {
        r[0] = sl->r[2 * lo]; r[1] = sl->r[2 * lo + 1];
        Ic[0] = sl->Ic[2 * lo]; Ic[1] = sl->Ic[2 * lo + 1];
        return GIA_OK;
    }
    st = gl5_step(sl, sl->t[lo], t, &sl->r[2 * lo], dI, r, &move, &dbl);
    if (st != GIA_OK) return st;
    Ic[0] = sl->Ic[2 * lo] + dI[0];
    Ic[1] = sl->Ic[2 * lo + 1] + dI[1];
    return GIA_OK;
}

gia_status gia_lde2_terms(const gia_lde2_sol *sol, double t,
                          double complex c[2], double complex r[2],
                          double complex E[2], const char **why) {
    double complex rr[2], Ic[2];
    int            i;
    if (!sol || !c || !r || !E)
        return fail(GIA_E_ARG, "gia_lde2_terms: NULL argument", why);
    if (!(t >= 0.0 && t <= sol->t_max))
        return fail(GIA_E_ARG, "gia_lde2_terms: t outside [0, t_max]", why);
    if (state_at(sol, t, rr, Ic) != GIA_OK)
        return fail(GIA_E_DOMAIN, "gia_lde2_terms: a coefficient is not finite", why);
    for (i = 0; i < 2; i++)
        if (creal(Ic[i]) > EXP_MAX)
            return fail(GIA_E_RANGE, "gia_lde2_terms: exp(int r) overflows (numerics N6)", why);
    for (i = 0; i < 2; i++) { c[i] = sol->c[i]; r[i] = rr[i]; E[i] = cexp(Ic[i]); }
    return GIA_OK;
}

/* r_i'(t) by central difference of the labelled roots, one-sided and second
 * order at either end of [0, t_max] (vv-plan.md §3: h = 1e-5 max(1, t)). */
static gia_status root_slopes(const gia_lde2_sol *sl, double t, const double complex r[2],
                              double complex dr[2]) {
    double         h = 1e-5 * (t > 1.0 ? t : 1.0);
    double complex p[2], q[2];
    gia_status     st;
    int            i;
    if (t - h >= 0.0 && t + h <= sl->t_max) {
        if ((st = roots_at(sl, t + h, p, NULL, NULL, NULL)) != GIA_OK) return st;
        if ((st = roots_at(sl, t - h, q, NULL, NULL, NULL)) != GIA_OK) return st;
        label(r, p); label(r, q);
        for (i = 0; i < 2; i++) dr[i] = (p[i] - q[i]) / (2.0 * h);
    } else {
        double s = (t - h < 0.0) ? h : -h;     /* step into the interval */
        if ((st = roots_at(sl, t + s, p, NULL, NULL, NULL)) != GIA_OK) return st;
        if ((st = roots_at(sl, t + 2.0 * s, q, NULL, NULL, NULL)) != GIA_OK) return st;
        label(r, p); label(p, q);
        for (i = 0; i < 2; i++) dr[i] = (-3.0 * r[i] + 4.0 * p[i] - q[i]) / (2.0 * s);
    }
    return GIA_OK;
}

gia_status gia_lde2_eval(const gia_lde2_sol *sol, double t,
                         double complex *f, double complex *trad_residual,
                         const char **why) {
    double complex c[2], r[2], E[2], dr[2], a1, a0, fv, res = 0.0;
    gia_status     st;
    int            i, nt;

    if (!sol) return fail(GIA_E_ARG, "gia_lde2_eval: sol is NULL", why);
    if ((st = gia_lde2_terms(sol, t, c, r, E, why)) != GIA_OK) return st;
    nt = sol->one_family ? 1 : 2;
    fv = 0.0;
    for (i = 0; i < nt; i++) fv += c[i] * E[i];
    if (trad_residual) {
        double complex unused[2];
        /* r is the terms' labelled roots; the slopes are labelled from them. */
        if (roots_at(sol, t, unused, &a1, &a0, NULL) != GIA_OK ||
            root_slopes(sol, t, r, dr) != GIA_OK)
            return fail(GIA_E_DOMAIN, "gia_lde2_eval: a coefficient is not finite", why);
        /* f'' + a1 f' + a0 f, termwise and traditional: d/dt (r E) = (r' + r^2) E. */
        for (i = 0; i < nt; i++)
            res += c[i] * (dr[i] + r[i] * r[i] + a1 * r[i] + a0) * E[i];
    }
    if (!finite_c(fv) || (trad_residual && !finite_c(res)))
        return fail(GIA_E_RANGE, "gia_lde2_eval: result not finite (NFR-NUM-003)", why);
    if (f) *f = fv;
    if (trad_residual) *trad_residual = res;
    return GIA_OK;
}
