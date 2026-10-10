/* test_mop.c — the Giannantoni kernel's verification and validation suite.
 *
 * docs/requirements/vv-plan.md §7 is the catalogue; every test here carries a
 * `Verifies:` tag naming its requirements and catalogue ID, and opens with the
 * source equation it rests on. The oracle policy (vv-plan.md §2, ADR 0018) is
 * strict: an oracle is a source equation's residual, a hand-derived closed
 * form, or a number printed in a source -- never output the code produced.
 *
 * Tolerances are vv-plan.md §3's and are never widened without a PLAN.md §5
 * erratum. A "must fail" assertion clears its tolerance by at least 100x.
 */

#define _POSIX_C_SOURCE 200809L

#include "engine.h"
#include "gia_status.h"
#include "idc.h"
#include "mop.h"
#include "relational.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* vv-plan.md §3 */
#define TOL_CLOSED   1e-12   /* closed form, well-conditioned: relative      */
#define TOL_RESIDUAL 1e-6    /* defining-equation residual: relative         */
#define TOL_QUAD     1e-9    /* quadrature: relative                         */

static int failures = 0;

static void ok(const char *what, int cond) {
    printf("  %-64s %s\n", what, cond ? "PASS" : "FAIL");
    if (!cond) failures++;
}

/* |got - want| <= tol * max(1, |want|) */
static int near_c(double complex got, double complex want, double tol) {
    double scale = cabs(want) > 1.0 ? cabs(want) : 1.0;
    return isfinite(creal(got)) && isfinite(cimag(got)) &&
           cabs(got - want) <= tol * scale;
}

static void close_c(const char *what, double complex got, double complex want,
                    double tol) {
    int c = near_c(got, want, tol);
    printf("  %-64s %s", what, c ? "PASS" : "FAIL");
    if (!c)
        printf("  (got %.15g%+.15gi want %.15g%+.15gi)", creal(got), cimag(got),
               creal(want), cimag(want));
    printf("\n");
    if (!c) failures++;
}

/* Sentinels: a value no function could produce, to show outputs untouched. */
#define SENTINEL_C (12345.0 - 678.0 * I)
static const char *const SENTINEL_WHY = "sentinel: why not written";

/* ------------------------------------------------------------------ *
 * IF-API-001 — status codes and the `why` out-parameter
 * ------------------------------------------------------------------ */

/* Source: docs/requirements/icd.md IF-API-001 */
/* Verifies: IF-API-001 (T-API-01) */
static void test_status_contract(void) {
    gia_status   s;
    int          all_named = 1, distinct = 1;
    const char  *why;
    double complex out;

    printf("\n[T-API-01] status codes and the why out-parameter\n");
    for (s = GIA_OK; s <= GIA_E_NOMEM; s++) {
        gia_status t;
        const char *n = gia_status_str(s);
        if (!n || !*n || strcmp(n, "unknown status") == 0) all_named = 0;
        for (t = GIA_OK; t < s; t++)
            if (n && gia_status_str(t) && strcmp(n, gia_status_str(t)) == 0) distinct = 0;
    }
    ok("every status has its own non-empty string", all_named && distinct);
    ok("an out-of-range status reads \"unknown status\"",
       strcmp(gia_status_str((gia_status)99), "unknown status") == 0);

    /* Failure: why set, output untouched. */
    why = SENTINEL_WHY; out = SENTINEL_C;
    s = gia_idc_of(0.0, 1.0, 2, &out, &why);
    ok("failure sets why to a reason", s != GIA_OK && why != SENTINEL_WHY && why && *why);
    ok("failure leaves the output untouched", out == SENTINEL_C);
    ok("failure with why == NULL is safe", gia_idc_of(0.0, 1.0, 2, &out, NULL) == GIA_E_DOMAIN);

    /* Success: why untouched. */
    why = SENTINEL_WHY;
    s = gia_idc_of(2.0, 1.0, 1, &out, &why);
    ok("success leaves why untouched", s == GIA_OK && why == SENTINEL_WHY);
}

/* ------------------------------------------------------------------ *
 * 7.1 Incipient calculus
 * ------------------------------------------------------------------ */

/* Source: [02 Eq 14.9.5], [23 Eq 5.5.2] — (d~/dt)^n f = (f'/f)^n f.
 *
 * Oracle: f = 1 + t^2 by hand, so f' = 2t and the incipient derivative is
 * (2t/(1+t^2))^n (1+t^2), written here as pow() of the closed-form ratio --
 * an independent route from the function's repeated multiplication. Then a
 * complex f, whose oracle is the polar form, and the refusal at f = 0. */
/* Verifies: FR-IDC-002 (T-IDC-03) */
static void test_idc_general_f(void) {
    static const double ts[] = { -2.5, -1.0, -0.3, 0.0, 0.4, 1.0, 3.0, 17.0 };
    size_t i;
    int    n, all = 1, refused, ratio_matters = 1, shown = 0;
    char   label[96];
    double complex out;
    const char *why = NULL;

    printf("\n[T-IDC-03] incipient derivative of a general function\n");
    for (i = 0; i < sizeof ts / sizeof ts[0]; i++) {
        double t = ts[i], f = 1.0 + t * t, df = 2.0 * t;
        for (n = 0; n <= 8; n++) {
            double want = pow(df / f, (double)n) * f;
            out = SENTINEL_C;
            if (gia_idc_of(f, df, n, &out, &why) != GIA_OK ||
                !near_c(out, want, TOL_CLOSED)) {
                if (shown++ < 3) {   /* the first few say where; the summary counts */
                    snprintf(label, sizeof label, "  f = 1+t^2, t = %g, n = %d", t, n);
                    close_c(label, out, want, TOL_CLOSED);
                    failures--;
                }
                all = 0;
            }
        }
    }
    ok("f = 1+t^2: (2t/(1+t^2))^n (1+t^2), 8 times x n = 0..8, 1e-12", all);

    /* Mutation "drop /f": with f != 1 the result must differ from df^n f. */
    for (n = 2; n <= 8; n++) {
        double t = 3.0, f = 1.0 + t * t, df = 2.0 * t;
        (void)gia_idc_of(f, df, n, &out, &why);
        if (near_c(out, pow(df, (double)n) * f, 100 * TOL_CLOSED)) ratio_matters = 0;
    }
    ok("the ratio f'/f, not f', is raised (\"drop /f\" would fail)", ratio_matters);

    /* Complex f: (f'/f)^n f in polar form, |r|^n e^{i n arg r} * f. */
    {
        double complex f = 1.5 - 0.75 * I, df = -0.2 + 2.0 * I, r = df / f;
        all = 1;
        for (n = 0; n <= 8; n++) {
            double complex want = pow(cabs(r), n) * cexp(I * n * carg(r)) * f;
            if (gia_idc_of(f, df, n, &out, &why) != GIA_OK || !near_c(out, want, TOL_CLOSED))
                all = 0;
        }
        ok("complex f: matches |f'/f|^n e^{i n arg(f'/f)} f, n = 0..8", all);
    }

    /* n = 0 returns f itself, whatever f' is. */
    ok("n = 0 is f", gia_idc_of(-4.0 + 2.0 * I, 1e300, 0, &out, &why) == GIA_OK &&
                     out == -4.0 + 2.0 * I);

    /* A negative f keeps its sign: no clamping (NFR-NUM-006). */
    ok("negative f, odd n: sign survives",
       gia_idc_of(-2.0, 1.0, 3, &out, &why) == GIA_OK &&
       near_c(out, pow(-0.5, 3) * -2.0, TOL_CLOSED));

    /* Refusals. */
    why = NULL;
    refused = gia_idc_of(0.0, 1.0, 2, &out, &why) == GIA_E_DOMAIN;
    ok("f = 0 is refused with GIA_E_DOMAIN", refused);
    ok("  and the reason cites [02 Eq 14.9.5]", why && strstr(why, "14.9.5") != NULL);
    ok("n < 0 is GIA_E_ARG", gia_idc_of(1.0, 1.0, -1, &out, &why) == GIA_E_ARG);
    ok("out == NULL is GIA_E_ARG", gia_idc_of(1.0, 1.0, 1, NULL, &why) == GIA_E_ARG);
    ok("overflow is GIA_E_RANGE, never inf",
       gia_idc_of(1e-300, 1e300, 3, &out, &why) == GIA_E_RANGE);
    ok("non-finite input is GIA_E_DOMAIN",
       gia_idc_of(NAN, 1.0, 1, &out, &why) == GIA_E_DOMAIN &&
       gia_idc_of(1.0, INFINITY, 1, &out, &why) == GIA_E_DOMAIN);
}

/* Coefficient functions for the LDE tests. ctx is unused except where a
 * test counts calls. */
static double complex zero_fn(double t, void *ctx)    { (void)t; (void)ctx; return 0.0; }
static double complex neg_sq1p(double t, void *ctx)   { (void)ctx; return -(1.0 + t) * (1.0 + t); }
static double complex neg2t(double t, void *ctx)      { (void)ctx; return -2.0 * t; }
static double complex tsq(double t, void *ctx)        { (void)ctx; return t * t; }
static double complex neg_sq_half(double t, void *ctx){ (void)ctx; return -(t - 0.5) * (t - 0.5); }
static double complex onept(double t, void *ctx)      { (void)ctx; return 1.0 + t; }
static double complex neg_tsq(double t, void *ctx)    { (void)ctx; return -t * t; }
static double complex neg_sq1psin(double t, void *ctx){ (void)ctx; return -(1.0 + sin(t)) * (1.0 + sin(t)); }
static double complex neg_sq1pit(double t, void *ctx) { (void)ctx; return -(1.0 + I * t) * (1.0 + I * t); }

/* Source: [06 Eq 3.3-3.6], [09 Eq 3, 7], [10 Eq 8.1, 10.1]; PLAN R13, R16,
 * X12; numerics.md N3.
 *
 * (i) a1 = 0, a0 = -(1+t)^2: the incipient characteristic r^2 - (1+t)^2 = 0
 *     has roots r = +-(1+t), so by hand
 *         f = c+ e^{phi} + c- e^{-phi},   phi = t + t^2/2,
 *     with c+ + c- = f0 and c+ - c- = f1 (since r(0) = +-1). The termwise
 *     incipient residual is sum c (r^2 + a1 r + a0) E = 0, and the traditional
 *     residual is sum c r' E = c+ e^{phi} - c- e^{-phi}, because r' = +-1.
 *     (The catalogue writes this case as a0 = -t^2, roots +-t; those roots
 *     coincide at t = 0, where N3 itself refuses the initial conditions, so
 *     the shifted coefficient keeps the test's intent with ICs that
 *     determine the constants. vv-plan.md T-IDC-04 is revised to match.)
 * (ii) a1 = -2t, a0 = t^2: discriminant 4t^2 - 4t^2 = 0 everywhere, a double
 *     root r = t. f1 = r(0) f0 = 0 is solved as f0 e^{t^2/2}; f1 != 0 is
 *     refused, and [06 Eq 3.7]'s second solution is not used (X12).
 * (iii) a1 = 0, a0 = -(t - 1/2)^2: roots +-(t - 1/2), distinct at 0, collide
 *     at t = 1/2: refused, naming the time. */
/* Verifies: FR-IDC-006, NFR-NUM-001 (T-IDC-04) */
static void test_lde2(void) {
    static const double ts[] = { 0.0, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0 };
    gia_lde2_sol  *sol = NULL;
    const char    *why = NULL;
    double         t_fail = -1.0;
    size_t         i;
    int            all, res_ok, trad_ok, ic_ok, roots_ok;
    double complex c[2], r[2], E[2], f, trad;
    const double complex f0 = 2.0, f1 = 0.5;
    const double complex cp = (f0 + f1) / 2.0, cm = (f0 - f1) / 2.0;

    printf("\n[T-IDC-04] second-order incipient LDE, variable coefficients\n");

    /* (i) */
    ok("(i) a1 = 0, a0 = -(1+t)^2 solves",
       gia_lde2_solve(zero_fn, neg_sq1p, NULL, f0, f1, 2.0, &sol, &t_fail, &why) == GIA_OK && sol);
    if (sol) {
        all = res_ok = trad_ok = roots_ok = 1;
        for (i = 0; i < sizeof ts / sizeof ts[0]; i++) {
            double t = ts[i], phi = t + t * t / 2.0;
            double complex want = cp * exp(phi) + cm * exp(-phi);
            double complex a1 = 0.0, a0 = -(1.0 + t) * (1.0 + t), resid = 0.0;
            int k;
            if (gia_lde2_eval(sol, t, &f, &trad, &why) != GIA_OK ||
                gia_lde2_terms(sol, t, c, r, E, &why) != GIA_OK) { all = 0; continue; }
            if (!near_c(f, want, TOL_CLOSED)) all = 0;
            for (k = 0; k < 2; k++) {
                /* each term is one of the two hand-derived ones */
                double sgn = creal(r[k]) > 0 ? 1.0 : -1.0;
                if (!near_c(r[k], sgn * (1.0 + t), TOL_CLOSED) ||
                    !near_c(E[k], exp(sgn * phi), TOL_CLOSED) ||
                    !near_c(c[k], sgn > 0 ? cp : cm, TOL_CLOSED)) roots_ok = 0;
                resid += c[k] * (r[k] * r[k] + a1 * r[k] + a0) * E[k];
            }
            if (!(cabs(resid) <= TOL_CLOSED * cabs(want))) res_ok = 0;
            if (!near_c(trad, cp * exp(phi) - cm * exp(-phi), TOL_RESIDUAL)) trad_ok = 0;
        }
        ok("(i) f = c+ e^{t+t^2/2} + c- e^{-(t+t^2/2)} by hand, 7 times, 1e-12", all);
        ok("(i) roots +-(1+t), E = e^{+-phi}, constants from the ICs", roots_ok);
        ok("(i) termwise incipient residual sum c(r^2+a1 r+a0)E = 0", res_ok);
        ok("(i) traditional residual = c+ e^{phi} - c- e^{-phi} (not 0)", trad_ok);
        ic_ok = gia_lde2_terms(sol, 0.0, c, r, E, &why) == GIA_OK &&
                near_c(c[0] * E[0] + c[1] * E[1], f0, TOL_CLOSED) &&
                near_c(c[0] * r[0] * E[0] + c[1] * r[1] * E[1], f1, TOL_CLOSED);
        ok("(i) f(0) = f0 and f~'(0) = sum c r(0) = f1", ic_ok);
        ok("eval outside [0, t_max] is GIA_E_ARG",
           gia_lde2_eval(sol, 2.5, &f, NULL, &why) == GIA_E_ARG &&
           gia_lde2_eval(sol, -0.1, &f, NULL, &why) == GIA_E_ARG);
        gia_lde2_free(sol); sol = NULL;
    }

    /* (i), with coefficients the quadrature cannot integrate exactly:
     * a0 = -(1 + sin t)^2, roots +-(1 + sin t), psi = t + 1 - cos t; and a
     * complex one, a0 = -(1 + i t)^2, roots +-(1 + i t), psi = t + i t^2/2. */
    {
        int k2;
        for (k2 = 0; k2 < 2; k2++) {
            gia_cfn a0f = k2 == 0 ? neg_sq1psin : neg_sq1pit;
            all = 1; res_ok = 1;
            if (gia_lde2_solve(zero_fn, a0f, NULL, f0, f1, 6.0, &sol, NULL, &why) != GIA_OK) {
                ok(k2 == 0 ? "(i') a0 = -(1+sin t)^2 solves" : "(i'') a0 = -(1+it)^2 solves", 0);
                continue;
            }
            for (i = 0; i <= 24; i++) {
                double t = 0.25 * (double)i;
                double complex psi = k2 == 0 ? t + 1.0 - cos(t) : t + I * t * t / 2.0;
                double complex want = cp * cexp(psi) + cm * cexp(-psi), resid = 0.0;
                double complex a0v = a0f(t, NULL);
                int k;
                if (gia_lde2_eval(sol, t, &f, NULL, &why) != GIA_OK ||
                    gia_lde2_terms(sol, t, c, r, E, &why) != GIA_OK ||
                    !near_c(f, want, TOL_CLOSED)) all = 0;
                for (k = 0; k < 2; k++) resid += c[k] * (r[k] * r[k] + a0v) * E[k];
                if (!(cabs(resid) <= TOL_CLOSED * cabs(want))) res_ok = 0;
            }
            ok(k2 == 0 ? "(i') a0 = -(1+sin t)^2: f by hand on [0, 6], 1e-12"
                       : "(i'') complex a0 = -(1+it)^2: f by hand on [0, 6], 1e-12", all && res_ok);
            gia_lde2_free(sol); sol = NULL;
        }
    }

    /* (ii) */
    ok("(ii) double root, consistent ICs (f1 = r(0) f0 = 0) solves",
       gia_lde2_solve(neg2t, tsq, NULL, 3.0, 0.0, 2.0, &sol, NULL, &why) == GIA_OK && sol);
    if (sol) {
        all = 1;
        for (i = 0; i < sizeof ts / sizeof ts[0]; i++) {
            double t = ts[i];
            if (gia_lde2_eval(sol, t, &f, &trad, &why) != GIA_OK ||
                !near_c(f, 3.0 * exp(t * t / 2.0), TOL_CLOSED)) all = 0;
            /* traditional residual of f0 e^{t^2/2} under a1 = -2t, a0 = t^2:
             * (1 + t^2) - 2t^2 + t^2 = 1, times f */
            if (!near_c(trad, 3.0 * exp(t * t / 2.0), TOL_RESIDUAL)) all = 0;
        }
        ok("(ii) f = f0 e^{t^2/2}; traditional residual = f, by hand", all);
        gia_lde2_free(sol); sol = NULL;
    }
    why = NULL;
    ok("(ii) double root, inconsistent ICs refused (GIA_E_DOMAIN)",
       gia_lde2_solve(neg2t, tsq, NULL, 3.0, 1.0, 2.0, &sol, NULL, &why) == GIA_E_DOMAIN && !sol);
    ok("  and the reason cites X12, not [06 Eq 3.7]'s second solution",
       why && strstr(why, "X12") != NULL);

    /* (iii) */
    why = NULL; t_fail = -1.0;
    ok("(iii) roots colliding at t = 1/2: GIA_E_CONVERGENCE",
       gia_lde2_solve(zero_fn, neg_sq_half, NULL, 1.0, 0.2, 1.0, &sol, &t_fail, &why)
           == GIA_E_CONVERGENCE && !sol);
    ok("  with the collision time reported (|t* - 1/2| < 1e-6)", fabs(t_fail - 0.5) < 1e-6);

    /* Roots that coincide at t = 0 only (a0 = -t^2, the catalogue's original
     * case): [[1, 1], [r1(0), r2(0)]] is singular, so the ICs cannot fix c. */
    why = NULL;
    ok("roots coinciding at t = 0 only (a0 = -t^2): refused, GIA_E_DOMAIN",
       gia_lde2_solve(zero_fn, neg_tsq, NULL, 1.0, 0.0, 1.0, &sol, NULL, &why) == GIA_E_DOMAIN &&
       !sol && why != NULL);
    ok("t_max <= 0 or NULL coefficient is GIA_E_ARG",
       gia_lde2_solve(zero_fn, neg_sq1p, NULL, 1.0, 0.0, 0.0, &sol, NULL, &why) == GIA_E_ARG &&
       gia_lde2_solve(NULL, neg_sq1p, NULL, 1.0, 0.0, 1.0, &sol, NULL, &why) == GIA_E_ARG);
}

/* Source: [10 App. Eq 28-34] (catalogue reading; [10] is not in docs/).
 *
 * u^2 + psi_f u = 0 with psi_f = 1 + t: roots 0 and -(1+t), so
 * g = C1 + C2 e^{-(t + t^2/2)}. With f = e^{int psi_f} = e^{t + t^2/2}, the
 * product F = f g~' = C2 (-(1+t)) is not constant although g satisfies the
 * incipient condition: the paper's "24.1 does not imply 24". */
/* Verifies: FR-IDC-006, BR-001, BR-008 (T-IDC-05, VAL-05) */
static void test_lde2_zero_root(void) {
    gia_lde2_sol  *sol = NULL;
    const char    *why = NULL;
    const double   C1 = 1.0, C2 = 0.5;
    double complex c[2], r[2], E[2], g, F0 = 0.0, F1 = 0.0;
    int            all = 1, res_ok = 1;
    double         t;

    printf("\n[T-IDC-05 / VAL-05] [10 App. Eq 28-34]: a zero root\n");
    ok("g~'' + (1+t) g~' = 0 solves",
       gia_lde2_solve(onept, zero_fn, NULL, C1 + C2, -C2, 2.0, &sol, NULL, &why) == GIA_OK && sol);
    if (!sol) return;
    for (t = 0.0; t <= 2.0 + 1e-12; t += 0.25) {
        double phi = t + t * t / 2.0;
        double complex gp = 0.0, resid = 0.0, F;
        int k;
        if (gia_lde2_eval(sol, t, &g, NULL, &why) != GIA_OK ||
            gia_lde2_terms(sol, t, c, r, E, &why) != GIA_OK) { all = 0; continue; }
        if (!near_c(g, C1 + C2 * exp(-phi), TOL_CLOSED)) all = 0;
        for (k = 0; k < 2; k++) {
            gp    += c[k] * r[k] * E[k];
            resid += c[k] * (r[k] * r[k] + (1.0 + t) * r[k]) * E[k];
        }
        if (!(cabs(resid) <= TOL_CLOSED)) res_ok = 0;
        F = exp(phi) * gp;
        if (!near_c(F, -C2 * (1.0 + t), TOL_CLOSED)) all = 0;
        if (t == 0.0) F0 = F;
        if (t == 1.0) F1 = F;
    }
    ok("g = C1 + C2 e^{-(t+t^2/2)}; F = f g~' = -C2 (1+t), by hand", all);
    ok("the incipient condition holds (termwise residual 0)", res_ok);
    ok("F is not constant: |F(1) - F(0)| > 100x tolerance",
       cabs(F1 - F0) > 100.0 * TOL_CLOSED);
    gia_lde2_free(sol);
}

/* Source: [06b] "linearly dependent on initial conditions".
 * solve(ic1 + ic2) = solve(ic1) + solve(ic2), pointwise. */
/* Verifies: FR-IDC-012 (T-IDC-10) */
static void test_lde2_linear_in_ics(void) {
    gia_lde2_sol  *s1 = NULL, *s2 = NULL, *s12 = NULL;
    const char    *why = NULL;
    int            all = 1;
    double         t;

    printf("\n[T-IDC-10] linearity in the initial conditions: LDE (N3)\n");
    if (gia_lde2_solve(zero_fn, neg_sq1p, NULL, 2.0, 0.5, 2.0, &s1, NULL, &why) != GIA_OK ||
        gia_lde2_solve(zero_fn, neg_sq1p, NULL, -1.0 + 0.5 * I, 3.0, 2.0, &s2, NULL, &why) != GIA_OK ||
        gia_lde2_solve(zero_fn, neg_sq1p, NULL, 1.0 + 0.5 * I, 3.5, 2.0, &s12, NULL, &why) != GIA_OK)
        all = 0;
    for (t = 0.0; all && t <= 2.0 + 1e-12; t += 0.125) {
        double complex a, b, ab;
        (void)gia_lde2_eval(s1, t, &a, NULL, &why);
        (void)gia_lde2_eval(s2, t, &b, NULL, &why);
        (void)gia_lde2_eval(s12, t, &ab, NULL, &why);
        if (!near_c(ab, a + b, TOL_CLOSED)) all = 0;
    }
    ok("LDE: solve(ic1 + ic2) = solve(ic1) + solve(ic2), 17 times", all);
    gia_lde2_free(s1); gia_lde2_free(s2); gia_lde2_free(s12);
}

/* Source: [02 Eq 14.7.1-14.7.7], [06 Eq 3.10-3.15]; PLAN R9, X7; numerics N4.
 *
 * A = 1, B = -2: u^2 + u - 2 = (u - 1)(u + 2), so u = 1, -2 and the
 * exponents are u^2 = 1, 4. By hand, for ICs (F0, H0):
 *     c2 = (F0 - H0)/3, c1 = F0 - c2,
 *     f = c1 e^t + c2 e^{4t},   f^(1/2) = c1 e^t - 2 c2 e^{4t}.
 * The defining residual f' + A f^(1/2) + B f is formed from eval's own
 * outputs with f' by central difference (residual first, vv-plan §2). */
/* Verifies: FR-IDC-007 (T-IDC-06) */
static void test_binary(void) {
    const double complex F0[2] = { 2.0, -1.0 + 0.5 * I }, H0[2] = { 0.5, 3.0 };
    gia_binary_sol sol;
    const char    *why = NULL;
    double complex f[2], fh[2], fp[2], fm[2];
    int            all = 1, res_ok = 1, ic_ok, s, k;
    double         t;

    printf("\n[T-IDC-06] the binary function f' + A f^(1/2) + B f = 0\n");
    ok("A = 1, B = -2 solves", gia_binary_solve(1.0, -2.0, F0, H0, &sol, &why) == GIA_OK);
    for (k = 0; k <= 8; k++) {
        t = 0.125 * (double)k;
        if (gia_binary_eval(&sol, t, f, fh, &why) != GIA_OK) { all = res_ok = 0; continue; }
        for (s = 0; s < 2; s++) {
            double complex c2 = (F0[s] - H0[s]) / 3.0, c1 = F0[s] - c2;
            if (!near_c(f[s], c1 * exp(t) + c2 * exp(4.0 * t), TOL_CLOSED) ||
                !near_c(fh[s], c1 * exp(t) - 2.0 * c2 * exp(4.0 * t), TOL_CLOSED)) all = 0;
        }
        if (t > 0.0) {
            double h = 1e-5 * (t > 1.0 ? t : 1.0);
            (void)gia_binary_eval(&sol, t + h, fp, NULL, &why);
            (void)gia_binary_eval(&sol, t - h, fm, NULL, &why);
            for (s = 0; s < 2; s++) {
                double complex d = (fp[s] - fm[s]) / (2.0 * h);
                double complex res = d + 1.0 * fh[s] - 2.0 * f[s];
                if (!(cabs(res) <= TOL_RESIDUAL * cabs(d))) res_ok = 0;
            }
        }
    }
    ok("both branches: f and f^(1/2) by hand, 9 times, 1e-12", all);
    ok("residual f' + A f^(1/2) + B f <= 1e-6 (f' by central difference)", res_ok);
    ic_ok = gia_binary_eval(&sol, 0.0, f, fh, &why) == GIA_OK;
    for (s = 0; s < 2; s++)
        ic_ok = ic_ok && near_c(f[s], F0[s], TOL_CLOSED) && near_c(fh[s], H0[s], TOL_CLOSED);
    ok("the four initial conditions are reproduced", ic_ok);
    ok("one characteristic for both branches: u = {1, -2} (X7)",
       (near_c(sol.u[0], 1.0, TOL_CLOSED) && near_c(sol.u[1], -2.0, TOL_CLOSED)) ||
       (near_c(sol.u[1], 1.0, TOL_CLOSED) && near_c(sol.u[0], -2.0, TOL_CLOSED)));

    /* Complex coefficients: A = -(2+i), B = 1+i has roots 1 and 1+i. */
    res_ok = gia_binary_solve(-(2.0 + I), 1.0 + I, F0, H0, &sol, &why) == GIA_OK;
    for (k = 1; res_ok && k <= 8; k++) {
        double h = 1e-5;
        t = 0.125 * (double)k;
        (void)gia_binary_eval(&sol, t, f, fh, &why);
        (void)gia_binary_eval(&sol, t + h, fp, NULL, &why);
        (void)gia_binary_eval(&sol, t - h, fm, NULL, &why);
        for (s = 0; s < 2; s++) {
            double complex d = (fp[s] - fm[s]) / (2.0 * h);
            if (!(cabs(d - (2.0 + I) * fh[s] + (1.0 + I) * f[s]) <= TOL_RESIDUAL * cabs(d)))
                res_ok = 0;
        }
    }
    ok("complex A = -(2+i), B = 1+i: residual <= 1e-6 on both branches", res_ok);

    why = NULL;
    ok("u1 = u2 (A = 2, B = 1) is refused: GIA_E_UNSUPPORTED",
       gia_binary_solve(2.0, 1.0, F0, H0, &sol, &why) == GIA_E_UNSUPPORTED && why);
    ok("an overflowing exponential is GIA_E_RANGE",
       gia_binary_solve(1.0, -2.0, F0, H0, &sol, &why) == GIA_OK &&
       gia_binary_eval(&sol, 200.0, f, fh, &why) == GIA_E_RANGE);
}

/* Source: [06b] "linearly dependent on initial conditions" (the binary half). */
/* Verifies: FR-IDC-012 (T-IDC-10) */
static void test_binary_linear_in_ics(void) {
    const double complex a0[2] = { 2.0, -1.0 }, ah[2] = { 0.5, 3.0 };
    const double complex b0[2] = { -0.25 + I, 4.0 }, bh[2] = { 1.0, -2.0 * I };
    double complex s0[2], sh[2], fa[2], fb[2], fs[2], ha[2], hb[2], hs[2];
    gia_binary_sol sa, sb, ss;
    const char    *why = NULL;
    int            all, k, s;

    printf("\n[T-IDC-10] linearity in the initial conditions: binary (N4)\n");
    for (s = 0; s < 2; s++) { s0[s] = a0[s] + b0[s]; sh[s] = ah[s] + bh[s]; }
    all = gia_binary_solve(1.0, -2.0, a0, ah, &sa, &why) == GIA_OK &&
          gia_binary_solve(1.0, -2.0, b0, bh, &sb, &why) == GIA_OK &&
          gia_binary_solve(1.0, -2.0, s0, sh, &ss, &why) == GIA_OK;
    for (k = 0; all && k <= 16; k++) {
        double t = 0.0625 * (double)k;
        (void)gia_binary_eval(&sa, t, fa, ha, &why);
        (void)gia_binary_eval(&sb, t, fb, hb, &why);
        (void)gia_binary_eval(&ss, t, fs, hs, &why);
        for (s = 0; s < 2; s++)
            if (!near_c(fs[s], fa[s] + fb[s], TOL_CLOSED) ||
                !near_c(hs[s], ha[s] + hb[s], TOL_CLOSED)) all = 0;
    }
    ok("binary: solve(ic1 + ic2) = solve(ic1) + solve(ic2), 17 times", all);
}

/* Riccati coefficients. ctx points at {Q, R, P} for the constant cases. */
static double complex rc_Q(double t, void *c)  { (void)t; return ((const double *)c)[0]; }
static double complex rc_R(double t, void *c)  { (void)t; return ((const double *)c)[1]; }
static double complex rc_dR(double t, void *c) { (void)t; (void)c; return 0.0; }
static double complex rc_P(double t, void *c)  { (void)t; return ((const double *)c)[2]; }
static double complex rv_Q(double t, void *c)  { (void)c; return t; }
static double complex rv_R(double t, void *c)  { (void)c; return 1.0 + t; }
static double complex rv_dR(double t, void *c) { (void)c; (void)t; return 1.0; }
static double complex rv_P(double t, void *c)  { (void)c; (void)t; return 1.0; }

/* By hand, for constant Q, R, P: the LDE y'' + Q y' - P R y = 0 has roots
 * r1, r2 of r^2 + Q r - P R = 0; y = c1 e^{r1 t} + c2 e^{r2 t} with
 * c1 + c2 = 1 and c1 r1 + c2 r2 = R f0; f = y'/(R y). */
static double ric_closed(const double *qrp, double f0, double t) {
    double Q = qrp[0], R = qrp[1], P = qrp[2];
    double d = sqrt(Q * Q + 4.0 * P * R), r1 = (-Q + d) / 2.0, r2 = (-Q - d) / 2.0;
    double c2 = (r1 - R * f0) / (r1 - r2), c1 = 1.0 - c2;
    double y = c1 * exp(r1 * t) + c2 * exp(r2 * t);
    double dy = c1 * r1 * exp(r1 * t) + c2 * r2 * exp(r2 * t);
    return dy / (R * y);
}

/* The printed [06 Eq 3.17] substitution, y = f'/(f R), applied to the true f,
 * fed to [06 Eq 3.18]: R y'' - (R' - Q R) y' - P R^2 y, by central
 * differences (h = 1e-4, as in probes/riccati_substitution.py). */
static double printed_317_residual(const double *qrp, double f0, double t) {
    const double h = 1e-4, Q = qrp[0], R = qrp[1], P = qrp[2];
    double y[3];
    int k;
    for (k = 0; k < 3; k++) {
        double tt = t + (double)(k - 1) * h, f = ric_closed(qrp, f0, tt);
        y[k] = (P - Q * f - R * f * f) / (f * R);      /* f' from the ODE */
    }
    return fabs(R * (y[2] - 2.0 * y[1] + y[0]) / (h * h) + Q * R * (y[2] - y[0]) / (2.0 * h)
                - P * R * R * y[1]);
}

/* Source: [06 Eq 3.16-3.18]; PLAN R11, X2; numerics.md N5.
 *
 * Constant coefficients, two cases: Q = R = 1, P = 2 (the probe's), and
 * Q = 0.5, R = 2, P = 1.5 (R != 1, so dropping R from f = y'/(R y) shows).
 * Oracle: the closed form above, and the Riccati residual of eval's f with
 * f' by central difference. X2: the printed substitution leaves an Eq 3.18
 * residual > 100x tolerance. Variable coefficients Q = t, R = 1 + t, P = 1:
 * the residual eval reports equals the central-difference residual of its
 * own f (it is reported, not assumed zero). */
/* Verifies: FR-IDC-008 (T-IDC-07) */
static void test_riccati(void) {
    static const double cases[2][3] = { { 1.0, 1.0, 2.0 }, { 0.5, 2.0, 1.5 } };
    const double f0 = 0.5;
    gia_lde2_sol  *sol = NULL;
    const char    *why = NULL;
    double complex f, res, fp, fm;
    int            k, i, all, res_ok, rep_ok;
    double         worst_printed = 0.0;

    printf("\n[T-IDC-07] Riccati by linearisation\n");
    for (k = 0; k < 2; k++) {
        double *qrp = (double *)cases[k];
        char    label[96];
        all = res_ok = rep_ok = 1;
        if (gia_riccati_solve(rc_Q, rc_R, rc_dR, rc_P, qrp, f0, 2.0, &sol, NULL, &why) != GIA_OK) {
            ok("constant Q, R, P solves", 0);
            for (i = 1; i <= 15; i++)
                if (printed_317_residual(qrp, f0, 0.125 * (double)i) > worst_printed)
                    worst_printed = printed_317_residual(qrp, f0, 0.125 * (double)i);
            continue;
        }
        for (i = 1; i <= 15; i++) {
            double t = 0.125 * (double)i, h = 1e-5;
            if (gia_riccati_eval(sol, t, &f, &res, &why) != GIA_OK ||
                gia_riccati_eval(sol, t + h, &fp, NULL, &why) != GIA_OK ||
                gia_riccati_eval(sol, t - h, &fm, NULL, &why) != GIA_OK) { all = res_ok = 0; continue; }
            if (!near_c(f, ric_closed(qrp, f0, t), TOL_CLOSED)) all = 0;
            {
                double complex d = (fp - fm) / (2.0 * h);
                double complex r = d + qrp[0] * f + qrp[1] * f * f - qrp[2];
                if (!(cabs(r) <= TOL_RESIDUAL * (1.0 + cabs(d)))) res_ok = 0;
            }
            if (!(cabs(res) <= TOL_RESIDUAL)) rep_ok = 0;
        }
        /* X2 is a fact about the printed formula, independent of the solver. */
        for (i = 1; i <= 15; i++)
            if (printed_317_residual(qrp, f0, 0.125 * (double)i) > worst_printed)
                worst_printed = printed_317_residual(qrp, f0, 0.125 * (double)i);
        snprintf(label, sizeof label, "Q=%g R=%g P=%g: f = y'/(R y) by hand, 15 times, 1e-12",
                 qrp[0], qrp[1], qrp[2]);
        ok(label, all);
        snprintf(label, sizeof label, "Q=%g R=%g P=%g: Riccati residual of f <= 1e-6", qrp[0],
                 qrp[1], qrp[2]);
        ok(label, res_ok);
        ok("  and the residual eval reports is <= 1e-6", rep_ok);
        ok("  f(0) = f0", gia_riccati_eval(sol, 0.0, &f, NULL, &why) == GIA_OK &&
                          near_c(f, f0, TOL_CLOSED));
        gia_lde2_free(sol); sol = NULL;
    }
    printf("    (printed [06 Eq 3.17]: worst Eq 3.18 residual %.4g)\n", worst_printed);
    ok("X2: the printed y = f'/(f R) leaves an Eq 3.18 residual > 100x tol",
       worst_printed > 100.0 * TOL_RESIDUAL);

    /* Variable coefficients. */
    rep_ok = gia_riccati_solve(rv_Q, rv_R, rv_dR, rv_P, NULL, 0.25, 1.0, &sol, NULL, &why) == GIA_OK;
    for (i = 1; rep_ok && i <= 7; i++) {
        double t = 0.125 * (double)i, h = 1e-5;
        double complex d;
        if (gia_riccati_eval(sol, t, &f, &res, &why) != GIA_OK ||
            gia_riccati_eval(sol, t + h, &fp, NULL, &why) != GIA_OK ||
            gia_riccati_eval(sol, t - h, &fm, NULL, &why) != GIA_OK) { rep_ok = 0; break; }
        d = (fp - fm) / (2.0 * h);
        if (!near_c(res, d + t * f + (1.0 + t) * f * f - 1.0, TOL_RESIDUAL)) rep_ok = 0;
    }
    ok("Q = t, R = 1+t, P = 1: reported residual = residual of its own f", rep_ok);
    gia_lde2_free(sol); sol = NULL;

    ok("R(0) = 0 is refused", gia_riccati_solve(rc_Q, rc_dR, rc_dR, rc_P, (void *)cases[0], 1.0,
                                                1.0, &sol, NULL, &why) == GIA_E_DOMAIN && !sol);
}

/* Source: PLAN §6; [06 Eq 3.22] (X3), [06 Eq 3.25-3.27], [06 §4 (i)]. */
/* Verifies: FR-IDC-013 (T-IDC-13) */
static void test_idc_refusals(void) {
    static const struct { const char *name, *cite; } r[] = {
        { "riccati_duet", "3.22" }, { "abel_net", "3.25" }, { "solution_drift", "FR-IDC-011" } };
    const char *why;
    size_t      i;
    int         all = 1;

    printf("\n[T-IDC-13] refusals name their source\n");
    for (i = 0; i < sizeof r / sizeof r[0]; i++) {
        why = NULL;
        if (gia_idc_refuse(r[i].name, &why) != GIA_E_UNSUPPORTED || !why ||
            !strstr(why, r[i].cite)) {
            printf("    %s: %s\n", r[i].name, why ? why : "(no reason)");
            all = 0;
        }
    }
    ok("riccati_duet, abel_net, solution_drift: GIA_E_UNSUPPORTED + source", all);
    ok("an unknown feature is GIA_E_ARG", gia_idc_refuse("teleportation", &why) == GIA_E_ARG &&
                                          gia_idc_refuse(NULL, &why) == GIA_E_ARG);
}

/* Source: [09 Eq 10], [10 Eq 13]; PLAN R10; numerics.md N9.
 * Oracle: the direct sum f0 sum (a dt)^k/k!, evaluated independently of the
 * function's Horner scheme, for n = 0..20 over a spread of a dt. */
/* Verifies: FR-IDC-010 (T-IDC-09) */
static void test_idc_taylor(void) {
    static const double adt[] = { -3.0, -0.5, 0.0, 0.1, 0.8, 2.0, 7.5 };
    const char *why = NULL;
    double      out = 0.0;
    size_t      i;
    int         n, all = 1;

    printf("\n[T-IDC-09] incipient Taylor projection\n");
    for (i = 0; i < sizeof adt / sizeof adt[0]; i++)
        for (n = 0; n <= 20; n++) {
            double f0 = -1.7, dt = 2.0, a = adt[i] / dt, sum = 0.0, term = 1.0;
            int k;
            for (k = 0; k <= n; k++) { sum += term; term *= adt[i] / (double)(k + 1); }
            if (gia_idc_taylor(f0, a * f0, dt, n, &out, &why) != GIA_OK ||
                !near_c(out, f0 * sum, TOL_CLOSED)) all = 0;
        }
    ok("Horner = direct sum, 7 values of a dt x n = 0..20, 1e-12", all);
    ok("f(t0) = 0 is refused (GIA_E_DOMAIN)", gia_idc_taylor(0.0, 1.0, 1.0, 2, &out, &why) == GIA_E_DOMAIN);
    ok("n > 170 is GIA_E_LIMIT; n < 0 is GIA_E_ARG",
       gia_idc_taylor(1.0, 1.0, 1.0, 171, &out, &why) == GIA_E_LIMIT &&
       gia_idc_taylor(1.0, 1.0, 1.0, -1, &out, &why) == GIA_E_ARG);
    ok("overflow is GIA_E_RANGE", gia_idc_taylor(1e300, 1e300, 1e6, 2, &out, &why) == GIA_E_RANGE);
}

/* Agreement with a number printed to `digits` decimal places: within one
 * unit of the last printed digit (the source truncates as often as it
 * rounds: 172.0667 is printed 172.06). */
static int printed_as(double got, double printed, int digits) {
    return fabs(got - printed) < pow(10.0, -digits);
}

/* Source: [09 Eq 16-22] (Giannantoni & Zoli 2009), with n = 2 (PLAN R10);
 * errata X4, X5. Inputs are those of probes/zoli_2009_reproduction.py. */
/* Verifies: BR-001, BR-002, BR-008, FR-IDC-010 (VAL-01) */
static void test_val_zoli_2009(void) {
    const char *why = NULL;
    double      v16 = NAN, v17 = NAN, v21 = NAN, v19a = NAN, v19b = NAN, v22a = NAN, v22b = NAN;

    printf("\n[VAL-01] [09] Giannantoni & Zoli 2009, n = 2\n");
    (void)gia_idc_taylor(0.4, 0.32, 10.0, 2, &v16, &why);
    (void)gia_idc_taylor(0.4, 0.11, 10.0, 2, &v17, &why);
    (void)gia_idc_taylor(18.0, 6.0, 10.0, 2, &v21, &why);
    (void)gia_idc_taylor(6.0 * 1.8, 6.0, 10.0 - 1.8, 2, &v19a, &why);
    (void)gia_idc_taylor(6.0 * 2.0, 6.0, 10.0 - 2.0, 2, &v19b, &why);
    (void)gia_idc_taylor(0.6 * 2.0, 0.6, 10.0 - 2.0, 2, &v22a, &why);
    (void)gia_idc_taylor(0.6 * 1.8, 0.6, 10.0 - 1.8, 2, &v22b, &why);
    printf("    Eq16 %.10g  Eq17 %.10g  Eq21 %.10g  Eq19(1.8) %.10g  Eq19(2) %.10g\n"
           "    Eq22(2) %.10g  Eq22(1.8) %.10g\n", v16, v17, v21, v19a, v19b, v22a, v22b);
    ok("[09 Eq 16] max scenario: 16.4", printed_as(v16, 16.4, 1));
    ok("[09 Eq 17] min scenario: 3.01 (3.0125)", printed_as(v17, 3.01, 2) &&
                                               near_c(v17, 3.0125, TOL_CLOSED));
    ok("[09 Eq 21] sea level L1: 178.0", printed_as(v21, 178.0, 1));
    ok("[09 Eq 19] tau0 = 1.8: 172.06 (172.0667)", printed_as(v19a, 172.06, 2) &&
                                                  near_c(v19a, 172.0 + 1.0 / 15.0, TOL_CLOSED));
    ok("[09 Eq 22] 15-17 cm: 15.6 and 17.2067, each 15..17 to the printed digit",
       near_c(v22a, 15.6, TOL_CLOSED) && near_c(v22b, 17.2 + 1.0 / 150.0, TOL_CLOSED) &&
       floor(v22a + 0.5) >= 15.0 && floor(v22a + 0.5) <= 17.0 &&
       floor(v22b + 0.5) >= 15.0 && floor(v22b + 0.5) <= 17.0);
    /* X4: the same formula at tau0 = 2 gives 156.0, not the printed 154.3. */
    ok("X4: [09 Eq 19] at tau0 = 2 is 156.0", near_c(v19b, 156.0, TOL_CLOSED));
    ok("X4: the printed 154.3 does not follow (|156.0 - 154.3| > 100x tol)",
       near_c(v19b, 156.0, TOL_CLOSED) && fabs(v19b - 154.3) > 100.0 * TOL_CLOSED * 156.0);
    /* X5: the min scenario's net increase is f* - f0 = 2.6125; the printed
     * 1.91 subtracts 1.1 where the max scenario subtracts 0.4. */
    ok("X5: min-scenario net increase is 2.6125, not the printed 1.91",
       near_c(v17 - 0.4, 2.6125, TOL_CLOSED) && fabs((v17 - 0.4) - 1.91) > 0.5);
}

/* The defining residual of [02 Eq 14.10.1] for F = e^{ut} at t, each
 * incipient derivative taken by gia_idc_of on the function's own value and
 * slope: F (d~/dt)^2 (F^2) + A F^2 (d~/dt) F + B F^3, divided by F^3. */
static double complex nl1410_residual(double complex u, double complex A, double complex B,
                                      double t) {
    double complex F = cexp(u * t), dF = u * F, F2 = F * F, dF2 = 2.0 * u * F2;
    double complex d2F2 = 0.0, d1F = 0.0;
    (void)gia_idc_of(F2, dF2, 2, &d2F2, NULL);
    (void)gia_idc_of(F, dF, 1, &d1F, NULL);
    return (F * d2F2 + A * F2 * d1F + B * F * F2) / (F * F2);
}

/* Source: [02 Eq 14.10.1-14.10.2]; PLAN X8.
 * Oracle: the equation's incipient residual for F = e^{ut}, formed with
 * gia_idc_of (FR-IDC-002, already verified) at several t, for each returned
 * root; and the count: a quadratic has two roots, and the [02 Eq 14.10.2]
 * "triplet" has no third. */
/* Verifies: FR-IDC-009 (T-IDC-08) */
static void test_nl1410(void) {
    static const double complex AB[3][2] = { { 3.0, -1.0 }, { -2.0 + I, 0.5 }, { 1.0, 5.0 } };
    double complex u[2];
    const char    *why = NULL;
    int            k, all = 1, distinct = 1, third = 0;

    printf("\n[T-IDC-08] [02 Eq 14.10.1]: 4u^2 + Au + B = 0\n");
    for (k = 0; k < 3; k++) {
        int    i;
        double t;
        if (gia_nl1410_roots(AB[k][0], AB[k][1], u, &why) != GIA_OK) { all = distinct = 0; continue; }
        if (cabs(u[0] - u[1]) < 1e-6) distinct = 0;
        for (i = 0; i < 2; i++)
            for (t = 0.0; t <= 1.0; t += 0.25)
                if (!(cabs(nl1410_residual(u[i], AB[k][0], AB[k][1], t)) <=
                      TOL_CLOSED * (1.0 + cabs(u[i]) * cabs(u[i])))) all = 0;
        /* No third solution: u = 0 and u = -A/4 (the remaining natural
         * candidates for a "triplet") leave a residual. */
        if (cabs(nl1410_residual(0.0, AB[k][0], AB[k][1], 0.5)) <= 1e-6 ||
            cabs(nl1410_residual(-AB[k][0] / 4.0, AB[k][0], AB[k][1], 0.5)) <= 1e-6) third = 1;
    }
    ok("each root: incipient residual of Eq 14.10.1 = 0, 3 (A, B) x 5 t, 1e-12", all);
    ok("two distinct roots for each (A, B)", distinct);
    ok("X8: no third solution (u = 0 and u = -A/4 leave a residual)", !third);
    ok("non-finite A is GIA_E_DOMAIN; NULL u is GIA_E_ARG",
       gia_nl1410_roots(NAN, 1.0, u, &why) == GIA_E_DOMAIN &&
       gia_nl1410_roots(1.0, 1.0, NULL, &why) == GIA_E_ARG);
}

/* ------------------------------------------------------------------ *
 * Network drift (FR-IDC-011, FR-IDC-014)
 * ------------------------------------------------------------------ */

/* Two storages, a -> b linear k = 0.5: a constant flow matrix with two modes. */
static const char *TWO_MODE =
    "{\"system_name\":\"two_mode\",\"nodes\":["
    " {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
    " {\"id\":\"b\",\"type\":\"storage\",\"current_level\":1.0}],"
    "\"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\",\"weight\":0.5}],"
    "\"simulation_params\":{\"t_val\":2.0,\"derivative_order\":2,\"generative_mode\":false}}";

/* A third constant-matrix seed: a chain with a held source and a sink. */
static const char *CHAIN =
    "{\"system_name\":\"chain\",\"nodes\":["
    " {\"id\":\"s\",\"type\":\"source\",\"initial_value\":3.0},"
    " {\"id\":\"x\",\"type\":\"storage\",\"current_level\":2.0},"
    " {\"id\":\"y\",\"type\":\"storage\",\"current_level\":0.5},"
    " {\"id\":\"z\",\"type\":\"sink\"}],"
    "\"edges\":[{\"source\":\"s\",\"target\":\"x\",\"weight\":0.3},"
    " {\"source\":\"x\",\"target\":\"y\",\"weight\":0.7},"
    " {\"source\":\"y\",\"target\":\"z\",\"weight\":0.2}],"
    "\"simulation_params\":{\"t_val\":3.0,\"derivative_order\":2,\"generative_mode\":false}}";

/* The example seed's shape: an interaction module, so the flow matrix depends
 * on the state (gia_flow_matrix_is_constant is false). */
static const char *INTERACTION =
    "{\"system_name\":\"interaction\",\"nodes\":["
    " {\"id\":\"source_1\",\"type\":\"source\",\"initial_value\":2.0},"
    " {\"id\":\"interaction_1\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
    " {\"id\":\"store_1\",\"type\":\"storage\",\"current_level\":10.0},"
    " {\"id\":\"consumer_1\",\"type\":\"storage\",\"current_level\":1.0}],"
    "\"edges\":["
    " {\"source\":\"source_1\",\"target\":\"interaction_1\",\"role\":\"energy\"},"
    " {\"source\":\"store_1\",\"target\":\"interaction_1\",\"role\":\"control\",\"use_ratio\":0},"
    " {\"source\":\"interaction_1\",\"target\":\"store_1\",\"weight\":1.0},"
    " {\"source\":\"store_1\",\"target\":\"consumer_1\",\"weight\":0.5}],"
    "\"simulation_params\":{\"t_val\":1.5,\"derivative_order\":2,\"generative_mode\":false}}";

static int load_seed(const char *json, gia_model *m, cJSON **root) {
    *root = cJSON_Parse(json);
    memset(m, 0, sizeof(*m));
    if (!*root) return 0;
    if (!gia_model_load(m, *root)) { cJSON_Delete(*root); *root = NULL; return 0; }
    return 1;
}

/* Source: [06 §4 (i)]; PLAN E4b, R16.
 * Oracle: exactly 0.0 (not "small"), for every node and every order 1..4, on
 * every constant-matrix seed here; and the refusal on a non-constant one.
 * The negative control shows why E4b had to be guarded: drift taken from
 * phi = ln Q_b on the two-mode network is phi'' != 0, so the mutation
 * "phi = ln Q_i" would report drift where [06 §4 (i)] says there is none. */
/* Verifies: FR-IDC-011, FR-IDC-013 (T-IDC-11, T-IDC-13) */
static void test_solution_drift(void) {
    const char *seeds[2] = { TWO_MODE, CHAIN };
    gia_model   m;
    cJSON      *root;
    const char *why = NULL;
    double      out[8], q[2], qp[2], qm[2];
    int         k, n, i, all = 1;

    printf("\n[T-IDC-11] solution drift on network trajectories\n");
    for (k = 0; k < 2; k++) {
        if (!load_seed(seeds[k], &m, &root)) { all = 0; continue; }
        if (!gia_flow_matrix_is_constant(&m)) all = 0;
        for (n = 1; n <= 4; n++) {
            for (i = 0; i < 8; i++) out[i] = 99.0;
            if (gia_solution_drift(&m, n, out, &why) != GIA_OK) { all = 0; continue; }
            for (i = 0; i < m.n_nodes; i++) if (out[i] != 0.0) all = 0;
        }
        gia_model_free(&m); cJSON_Delete(root);
    }
    ok("constant flow matrix: exactly 0.0, every node, n = 1..4, 2 seeds", all);

    /* The hazard E4b: phi = ln Q_b on the two-mode network. phi'' at t = 1 by
     * central difference of the engine's own Q_b. */
    if (load_seed(TWO_MODE, &m, &root)) {
        const double h = 1e-3;
        double phi2;
        (void)gia_network_state(&m, 1.0, q, NULL);
        (void)gia_network_state(&m, 1.0 + h, qp, NULL);
        (void)gia_network_state(&m, 1.0 - h, qm, NULL);
        phi2 = (log(qp[1]) - 2.0 * log(q[1]) + log(qm[1])) / (h * h);
        printf("    (phi = ln Q_b would give drift phi'' = %.6g at t = 1)\n", phi2);
        ok("negative control: phi = ln Q_i is not affine on a two-mode network",
           fabs(phi2) > 1e-3);
        ok("order outside 1..4 is GIA_E_ARG",
           gia_solution_drift(&m, 0, out, &why) == GIA_E_ARG &&
           gia_solution_drift(&m, 5, out, &why) == GIA_E_ARG);
        gia_model_free(&m); cJSON_Delete(root);
    }

    why = NULL;
    if (load_seed(INTERACTION, &m, &root)) {
        ok("non-constant flow matrix: refused, GIA_E_UNSUPPORTED",
           gia_solution_drift(&m, 2, out, &why) == GIA_E_UNSUPPORTED);
        ok("  and the reason cites FR-IDC-011 / [06 §4]",
           why && strstr(why, "FR-IDC-011") && strstr(why, "06"));
        gia_model_free(&m); cJSON_Delete(root);
    } else {
        ok("interaction seed loads", 0);
    }
}

/* Source: [09 Eq 9-13] at k = 2; numerics.md N7.
 * Two-mode network, by hand (k = 0.5, Q_a(0) = 10, Q_b(0) = 1):
 *   Q_a = 10 e^{-kt}, Q_b = 1 + 10 (1 - e^{-kt}),
 *   Q_a' = -k Q_a, Q_b' = k Q_a,  Q_a'' = k^2 Q_a,  Q_b'' = -k^2 Q_a,
 * so drift_a = 0 (a single exponential projects exactly) and
 * drift_b = (-k^2 Q_a - (k Q_a)^2 / Q_b) dt^2 / 2.
 * Interaction seed: Q' and Q'' against central differences of the engine's
 * own Q (h = 1e-3), and the drift scales as dt^2 (Q'' does not depend on dt). */
/* Verifies: FR-IDC-014 (T-IDC-12) */
static void test_drift_projection(void) {
    gia_model   m;
    cJSON      *root;
    const char *why = NULL;
    double      out[8], out2[8];
    const double k = 0.5, dt = 0.2;
    int         i, all = 1;

    printf("\n[T-IDC-12] output-projection drift\n");
    if (!load_seed(TWO_MODE, &m, &root)) { ok("two-mode seed loads", 0); return; }
    for (i = 0; i <= 8; i++) {
        double t = 0.25 * (double)i, qa = 10.0 * exp(-k * t);
        double qb = 1.0 + 10.0 * (1.0 - exp(-k * t));
        double want_b = (-k * k * qa - (k * qa) * (k * qa) / qb) * dt * dt / 2.0;
        if (gia_drift_projection(&m, t, dt, out, &why) != GIA_OK ||
            !near_c(out[0], 0.0, TOL_CLOSED * k * k * qa * dt * dt) ||
            !near_c(out[1], want_b, TOL_CLOSED)) {
            printf("    t = %g: got %.15g %.15g want 0 %.15g\n", t, out[0], out[1], want_b);
            all = 0;
        }
    }
    ok("constant matrix: drift_a = 0, drift_b by hand, 9 times, 1e-12", all);
    gia_model_free(&m); cJSON_Delete(root);

    if (!load_seed(INTERACTION, &m, &root)) { ok("interaction seed loads", 0); return; }
    all = 1;
    for (i = 1; i <= 5; i++) {
        double t = 0.25 * (double)i, h = 1e-3, q[4], qp[4], qm[4];
        int c;
        if (gia_drift_projection(&m, t, dt, out, &why) != GIA_OK ||
            gia_drift_projection(&m, t, dt / 2.0, out2, &why) != GIA_OK) { all = 0; continue; }
        (void)gia_network_state(&m, t, q, NULL);
        (void)gia_network_state(&m, t + h, qp, NULL);
        (void)gia_network_state(&m, t - h, qm, NULL);
        for (c = 2; c <= 3; c++) {            /* store_1, consumer_1 */
            double d1 = (qp[c] - qm[c]) / (2.0 * h);
            double d2 = (qp[c] - 2.0 * q[c] + qm[c]) / (h * h);
            double want = (d2 - d1 * d1 / q[c]) * dt * dt / 2.0;
            if (!(fabs(out[c] - want) <= 1e-4 * (fabs(want) + dt * dt))) {
                printf("    t = %g, node %d: got %.10g, finite differences %.10g\n", t, c,
                       out[c], want);
                all = 0;
            }
            if (!near_c(out2[c] * 4.0, out[c], 1e-9)) all = 0;
        }
        if (out[0] != 0.0 || out[1] != 0.0) all = 0;     /* source, module: not components */
    }
    ok("interaction: converges; agrees with differences of Q; scales as dt^2", all);
    ok("Q = 0 or a non-component: 0, not a division by zero",
       gia_drift_projection(&m, 0.0, dt, out, &why) == GIA_OK && out[1] == 0.0);
    gia_model_free(&m); cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 7.2 Emergy algebra in IDC form
 * ------------------------------------------------------------------ */

/* A source at transformity 1000 feeding a process that sends its output down
 * n pathways, all replicating (co-production) or all partitioning. */
static char *fan_seed(int n, const char *mode) {
    static char buf[16384];
    int  k, len;
    len = snprintf(buf, sizeof buf,
        "{\"nodes\":[{\"id\":\"sun\",\"type\":\"source\",\"value\":10.0,\"quality_input\":1000.0},"
        "{\"id\":\"proc\",\"type\":\"storage\",\"current_level\":5.0}");
    for (k = 0; k < n; k++)
        len += snprintf(buf + len, sizeof buf - (size_t)len,
                        ",{\"id\":\"out%d\",\"type\":\"storage\",\"current_level\":0.0}", k);
    len += snprintf(buf + len, sizeof buf - (size_t)len,
        "],\"edges\":[{\"source\":\"sun\",\"target\":\"proc\",\"weight\":1.0}");
    for (k = 0; k < n; k++)
        len += snprintf(buf + len, sizeof buf - (size_t)len,
                        ",{\"source\":\"proc\",\"target\":\"out%d\",\"weight\":0.5,"
                        "\"output_mode\":\"%s\"}", k, mode);
    snprintf(buf + len, sizeof buf - (size_t)len,
        "],\"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "\"generative_mode\":false}}");
    return buf;
}

/* Source: [02 Eq 3.6-3.8]; [02 p. 23 rule 2].
 * Oracle: a co-production with n products creates (n - 1) Em(u), where
 * Em(u) is the emergy the process receives; a partition creates none. The
 * process's source term, the model's emergy excess, and each product's
 * emergy (the whole Em(u)) are checked for n = 2, 3, 4. */
/* Verifies: FR-EM-002 (T-EM-03) */
static void test_em_coproduction(void) {
    int n, all = 1, part = 1;

    printf("\n[T-EM-03] co-production creates (n - 1) Em(u)\n");
    for (n = 2; n <= 4; n++) {
        gia_model   m;
        cJSON      *root;
        const char *why = NULL;
        double      phi = -1.0, em[8];
        if (!load_seed(fan_seed(n, "replicate"), &m, &root)) { all = 0; continue; }
        if (!gia_emergy_at(&m, 1.0, em, NULL) ||
            gia_emergy_source_term(&m, 1.0, 1, &phi, &why) != GIA_OK ||
            !near_c(phi, (n - 1) * em[1], TOL_CLOSED) ||
            !near_c(gia_emergy_excess(&m, 1.0), (n - 1) * em[1], TOL_CLOSED) ||
            !near_c(em[2], em[1], TOL_CLOSED)) {
            printf("    n = %d: phi %.10g excess %.10g want %.10g\n", n, phi,
                   gia_emergy_excess(&m, 1.0), (n - 1) * em[1]);
            all = 0;
        }
        gia_model_free(&m); cJSON_Delete(root);
        if (!load_seed(fan_seed(n, "partition"), &m, &root)) { part = 0; continue; }
        if (gia_emergy_source_term(&m, 1.0, 1, &phi, &why) != GIA_OK ||
            !near_c(phi, 0.0, TOL_CLOSED * 1e4)) part = 0;
        gia_model_free(&m); cJSON_Delete(root);
    }
    ok("co-production: Phi = excess = (n-1) Em(u), each product Em(u); n = 2,3,4", all);
    ok("partition: Phi = 0", part);
}

/* The measured gate of ADR 0017: sun (quality 1) as energy, H (quality 1000)
 * as a drawn control (use_ratio 0.01), F = k sun H = 0.3 * 10 * 2 = 6. */
static const char *MEASURED_GATE =
    "{\"nodes\":["
    " {\"id\":\"sun\",\"type\":\"source\",\"value\":10.0,\"quality_input\":1.0},"
    " {\"id\":\"H\",\"type\":\"source\",\"value\":2.0,\"quality_input\":1000.0},"
    " {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.3}},"
    " {\"id\":\"out\",\"type\":\"storage\",\"current_level\":0.0},"
    " {\"id\":\"heat\",\"type\":\"sink\"}],"
    "\"edges\":["
    " {\"source\":\"sun\",\"target\":\"gate\",\"role\":\"energy\"},"
    " {\"source\":\"H\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0.01},"
    " {\"source\":\"gate\",\"target\":\"out\",\"weight\":1.0},"
    " {\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"}],"
    "\"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,\"generative_mode\":false}}";

/* Source: [02 Eq 3.9, 3.12, 3.15]; ADR 0017.
 * Oracle, by hand: Em(u1) = F * 1 = 6 (energy), Em(u2) = s F * 1000 = 60
 * (drawn control), so Em(y) = Em(u1) + Em(u2) = 66 and the interaction's
 * source term Phi(u1, u2) = 0. The mutation "Em(y) = k Em1 Em2" would give
 * 0.3 * 6 * 60 = 108. */
/* Verifies: FR-EM-004 (T-EM-04) */
static void test_em_interaction(void) {
    gia_model   m;
    cJSON      *root;
    const char *why = NULL;
    double      em[5], phi = -1.0;

    printf("\n[T-EM-04] a drawn interaction: Em(y) = Em(u1) + Em(u2), Phi = 0\n");
    if (!load_seed(MEASURED_GATE, &m, &root)) { ok("measured gate loads", 0); return; }
    ok("Em(y) = 6 + 60 = 66", gia_emergy_at(&m, 1.0, em, NULL) && near_c(em[3], 66.0, 1e-9));
    ok("Phi(gate) = 0", gia_emergy_source_term(&m, 1.0, 2, &phi, &why) == GIA_OK &&
                        fabs(phi) <= 1e-9 * 66.0);
    ok("a boundary node is not a process: GIA_E_ARG",
       gia_emergy_source_term(&m, 1.0, 0, &phi, &why) == GIA_E_ARG &&
       gia_emergy_source_term(&m, 1.0, 9, &phi, &why) == GIA_E_ARG);
    gia_model_free(&m); cJSON_Delete(root);
}

/* Source: [02 Eq 3.23-3.26] on the totals of [02 Fig. 3.4] (Brown 1993):
 * S = 10,000, F = 20,000, Z = 30,000, Y = 7,500.
 *   Case A, 1 S + 1 F = 1/2 Z + 1/2 * 4 Y: 30,000 = 15,000 + 15,000.
 *   Case B, S + F + Phi_D + Phi_E = Z + 6 Y with Phi_E = Phi_D / 2:
 *     Phi_D + Phi_E = 75,000 - 30,000 = 45,000, so Phi_D = 30,000, Phi_E = 15,000. */
/* Verifies: FR-EM-007, BR-008 (T-EM-07, VAL-04) */
static void test_em_global_balance(void) {
    const gia_balance_term inA[2]  = { { 10000.0, 1.0 }, { 20000.0, 1.0 } };
    const gia_balance_term outA[2] = { { 30000.0, 0.5 }, { 7500.0, 0.5 * 4.0 } };
    const gia_balance_term outB[2] = { { 30000.0, 1.0 }, { 7500.0, 6.0 } };
    const double           w[2]    = { 1.0, 0.5 };
    double                 res = -1.0, phi[2] = { 0.0, 0.0 };
    const char            *why = NULL;

    printf("\n[T-EM-07 / VAL-04] [02 Eq 3.23-3.26] global balance, Fig. 3.4 totals\n");
    ok("case A balances: residual 0",
       gia_emergy_global_balance(inA, 2, outA, 2, &res, &why) == GIA_OK && res == 0.0);
    ok("case B: Phi_D = 30,000, Phi_E = 15,000",
       gia_emergy_balance_solve(inA, 2, outB, 2, w, 2, phi, &why) == GIA_OK &&
       phi[0] == 30000.0 && phi[1] == 15000.0);
    ok("no source terms, or weights summing to 0: GIA_E_ARG / GIA_E_DOMAIN",
       gia_emergy_balance_solve(inA, 2, outB, 2, w, 0, phi, &why) == GIA_E_ARG &&
       gia_emergy_balance_solve(inA, 2, outB, 2, (const double[]){ 1.0, -1.0 }, 2, phi, &why)
           == GIA_E_DOMAIN);
}

/* Source: [02 p. 23 rules 1-4], [02 Eq 3.8, 3.12, 3.16-3.17] — through the
 * network engine: the four emergy rules in one model. A split shares emergy
 * in proportion to flow (rule 3); a co-production gives each product the whole
 * (rule 2); a drawn interaction sums its inputs (Eq 3.12); co-products
 * reunited at one component count once (rule 4a). */
/* Verifies: BR-003 (VAL-03) */
static void test_val_emergy_rules(void) {
    gia_model   m;
    cJSON      *root;
    double      em[8], q[8];
    int         good = 1;
    /* sun -> p; p splits to a and b (partition); a replicates to c1, c2; both
     * c1 and c2 flow into r (reunion). */
    static const char *seed =
        "{\"nodes\":["
        " {\"id\":\"sun\",\"type\":\"source\",\"value\":10.0,\"quality_input\":100.0},"
        " {\"id\":\"p\",\"type\":\"storage\",\"current_level\":4.0},"
        " {\"id\":\"a\",\"type\":\"storage\",\"current_level\":2.0},"
        " {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0},"
        " {\"id\":\"c1\",\"type\":\"storage\",\"current_level\":1.0},"
        " {\"id\":\"c2\",\"type\":\"storage\",\"current_level\":1.0},"
        " {\"id\":\"r\",\"type\":\"storage\",\"current_level\":0.0}],"
        "\"edges\":["
        " {\"source\":\"sun\",\"target\":\"p\",\"weight\":1.0},"
        " {\"source\":\"p\",\"target\":\"a\",\"weight\":0.3},"
        " {\"source\":\"p\",\"target\":\"b\",\"weight\":0.1},"
        " {\"source\":\"a\",\"target\":\"c1\",\"weight\":0.2,\"output_mode\":\"replicate\"},"
        " {\"source\":\"a\",\"target\":\"c2\",\"weight\":0.2,\"output_mode\":\"replicate\"},"
        " {\"source\":\"c1\",\"target\":\"r\",\"weight\":0.5},"
        " {\"source\":\"c2\",\"target\":\"r\",\"weight\":0.5}],"
        "\"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,\"generative_mode\":false}}";

    printf("\n[VAL-03] [02] the four emergy rules, through the network engine\n");
    if (!load_seed(seed, &m, &root)) { ok("seed loads", 0); return; }
    good = gia_emergy_at(&m, 1.0, em, NULL) && gia_network_state(&m, 1.0, q, NULL);
    /* Rule 1: the sun's emergy reaches p: Em(p) = F(sun->p) * 100 = q_sun * 1 * 100. */
    ok("rule 1: Em(p) = 10 * 1.0 * 100 = 1000", good && near_c(em[1], 1000.0, 1e-9));
    /* Rule 3: a and b split p's emergy 0.3 : 0.1 of p's outflow. */
    ok("rule 3: Em(a) : Em(b) = 0.3 : 0.1, summing to Em(p)",
       good && near_c(em[2], 750.0, 1e-9) && near_c(em[3], 250.0, 1e-9));
    /* Rule 2: each co-product carries the whole of a's emergy. */
    ok("rule 2: Em(c1) = Em(c2) = Em(a)", good && near_c(em[4], em[2], 1e-9) &&
                                          near_c(em[5], em[2], 1e-9));
    /* Rule 4a: reunited at r, the co-products count once -- c1 and c2 each pass
     * all their emergy to r, and r takes the larger, not the sum. */
    ok("rule 4: Em(r) = max, not sum: Em(r) = Em(a)", good && near_c(em[6], em[2], 1e-9));
    gia_model_free(&m); cJSON_Delete(root);
}

/* Source: NFR-LIM-001 (BR-009). 65 inflows into one component, and 65
 * components with a co-production: refused with GIA_E_LIMIT, where the
 * baseline silently dropped the 65th inflow and ignored ancestry past 64. */
/* Verifies: NFR-LIM-001 (T-LIM-01) */
static void test_em_limits(void) {
    static char buf[16384];
    gia_model   m;
    cJSON      *root;
    const char *why = NULL;
    double      phi, em[70];
    int         k, len, refused, tidy = 1;

    printf("\n[T-LIM-01] emergy limits refuse, never truncate\n");
    /* 65 sources, each into `hub`. */
    len = snprintf(buf, sizeof buf, "{\"nodes\":[{\"id\":\"hub\",\"type\":\"storage\","
                   "\"current_level\":0.0}");
    for (k = 0; k < 65; k++)
        len += snprintf(buf + len, sizeof buf - (size_t)len,
                        ",{\"id\":\"s%d\",\"type\":\"source\",\"value\":1.0,"
                        "\"quality_input\":1.0}", k);
    len += snprintf(buf + len, sizeof buf - (size_t)len, "],\"edges\":[");
    for (k = 0; k < 65; k++)
        len += snprintf(buf + len, sizeof buf - (size_t)len,
                        "%s{\"source\":\"s%d\",\"target\":\"hub\",\"weight\":1.0}",
                        k ? "," : "", k);
    snprintf(buf + len, sizeof buf - (size_t)len,
             "],\"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
             "\"generative_mode\":false}}");
    if (!load_seed(buf, &m, &root)) { ok("65-inflow seed loads", 0); return; }
    why = NULL;
    refused = gia_emergy_source_term(&m, 1.0, 0, &phi, &why) == GIA_E_LIMIT &&
              gia_emergy_check_limits(&m, &why) == GIA_E_LIMIT &&
              !gia_emergy_at(&m, 1.0, em, NULL);
    ok("65 inflows into one node: GIA_E_LIMIT, and gia_emergy_at refuses", refused);
    ok("  and the reason names NFR-LIM-001", why && strstr(why, "NFR-LIM-001"));
    gia_model_free(&m); cJSON_Delete(root);

    /* 66 nodes (sun, proc, 64 products) with a co-production. */
    if (!load_seed(fan_seed(64, "replicate"), &m, &root)) { ok("66-node seed loads", 0); return; }
    ok("66 nodes with a co-production: GIA_E_LIMIT",
       gia_emergy_check_limits(&m, &why) == GIA_E_LIMIT);
    gia_model_free(&m); cJSON_Delete(root);

    /* At the limit, it works. */
    if (!load_seed(fan_seed(60, "replicate"), &m, &root)) { ok("62-node seed loads", 0); return; }
    tidy = gia_emergy_check_limits(&m, &why) == GIA_OK &&
           gia_emergy_source_term(&m, 1.0, 1, &phi, &why) == GIA_OK;
    ok("62 nodes, 60 co-products: within the limits, Phi computed", tidy);
    gia_model_free(&m); cJSON_Delete(root);
}

/* Source: [22 Eq 6-8], [06b Eq 6, 10].
 * Oracle: the shapes and entries the equations print. A binary is a column
 * of two equal branches, a duet a row of its two inputs, a duet-binary the
 * specular 2 x 2. The mutation "transpose" turns the binary into a row and
 * the duet into a column, and fails the shape checks. */
/* Verifies: FR-EM-005 (T-EM-05) */
static void test_em_ordinal_forms(void) {
    gia_oform b = gia_oform_binary(7.0), d = gia_oform_duet(3.0, 5.0);
    gia_oform f = gia_oform_duet_binary(2.0, -1.5);

    printf("\n[T-EM-05] ordinal forms of the three processes\n");
    ok("co-production: binary, a 2 x 1 column (Em(u); Em(u))",
       b.kind == GIA_OF_BINARY && b.rows == 2 && b.cols == 1 &&
       b.v[0][0] == 7.0 && b.v[1][0] == 7.0);
    ok("interaction: duet, a 1 x 2 row [Em(u1), Em(u2)]",
       d.kind == GIA_OF_DUET && d.rows == 1 && d.cols == 2 &&
       d.v[0][0] == 3.0 && d.v[0][1] == 5.0);
    ok("feedback: duet-binary [[a1, a2], [a2, a1]]",
       f.kind == GIA_OF_DUET_BINARY && f.rows == 2 && f.cols == 2 &&
       f.v[0][0] == 2.0 && f.v[0][1] == -1.5 && f.v[1][0] == -1.5 && f.v[1][1] == 2.0);
    ok("signs survive (no clamping, NFR-NUM-006)", f.v[0][1] < 0.0);
}

/* Source: [06b Eq 2], [02 Eq 14.11.4-14.11.5]; PLAN R12.
 * (a1; a2) o [b1, b2] = [(a1 b1; a2 b1), (a1 b2; a2 b2)]: unreduced, each
 * entry the pair of factors; reduced, their products. l o l is the pair
 * [l, l], and only its reduction is l^2. */
/* Verifies: FR-EM-006 (T-EM-06) */
static void test_em_circle_product(void) {
    gia_oform  a = gia_oform_binary(0.0), b = gia_oform_duet(5.0, 7.0), l = gia_oform_scalar(3.0);
    gia_oform  r;
    gia_circle c;
    const char *why = NULL;
    int        i, j, pairs = 1, red = 1;

    printf("\n[T-EM-06] the circle product and its cardinal reduction\n");
    a.v[0][0] = 2.0; a.v[1][0] = -3.0;            /* (a1; a2) = (2; -3) */
    ok("(a1; a2) o [b1, b2] is defined", gia_circle_product(&a, &b, &c, &why) == GIA_OK);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            if (c.pair[i][j][0] != a.v[i][0] || c.pair[i][j][1] != b.v[0][j]) pairs = 0;
    ok("unreduced: entry (i, j) keeps the pair (a_i, b_j) [06b Eq 2]",
       c.rows == 2 && c.cols == 2 && pairs);
    r = gia_circle_reduce(&c);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            if (r.v[i][j] != a.v[i][0] * b.v[0][j]) red = 0;
    ok("reduced: (a1 b1; a2 b1), (a1 b2; a2 b2)", r.rows == 2 && r.cols == 2 && red);

    ok("l o l is defined", gia_circle_product(&l, &l, &c, &why) == GIA_OK);
    ok("l o l is stored as the du-et [l, l], not as l^2 [02 Eq 14.11.5]",
       c.rows == 1 && c.cols == 1 && c.pair[0][0][0] == 3.0 && c.pair[0][0][1] == 3.0);
    r = gia_circle_reduce(&c);
    ok("and reduces to l^2 = 9", r.rows == 1 && r.cols == 1 && r.v[0][0] == 9.0);

    why = NULL;
    ok("a row on the left, or a column on the right: GIA_E_ARG",
       gia_circle_product(&b, &a, &c, &why) == GIA_E_ARG && why != NULL);
}

/* ------------------------------------------------------------------ *
 * 7.3 MOP — the First Fundamental Equation
 * ------------------------------------------------------------------ */

static gia_beta affine(double complex a, double complex b, double p) {
    gia_beta be;
    memset(&be, 0, sizeof(be));
    be.kind = GIA_BETA_AFFINE; be.a = a; be.b = b; be.p = p;
    return be;
}

static double complex beta_at(const gia_beta *be, double t) {
    return cpow(be->a + be->b * t, be->p);
}

/* (alpha'/alpha)^k alpha - beta at t, with alpha' by central difference
 * (h = 1e-5 max(1, t)) of gia_mop_couple's own output, relative to |beta|.
 * k integer uses gia_idc_of (FR-IDC-002); k = 1/2 the principal square root. */
static double first_eq_residual(const gia_beta *be, gia_rational k, double t) {
    double         h = 1e-5 * (t > 1.0 ? t : 1.0);
    double complex al, ap, am, inc = 0.0, b = beta_at(be, t);
    if (gia_mop_couple(be, k, t, &al, NULL) != GIA_OK ||
        gia_mop_couple(be, k, t + h, &ap, NULL) != GIA_OK ||
        gia_mop_couple(be, k, t - h, &am, NULL) != GIA_OK) return INFINITY;
    if (k.den == 1) {
        if (gia_idc_of(al, (ap - am) / (2.0 * h), k.num, &inc, NULL) != GIA_OK) return INFINITY;
    } else {
        inc = cpow((ap - am) / (2.0 * h) / al, (double)k.num / (double)k.den) * al;
    }
    return cabs(inc - b) / cabs(b);
}

/* Source: [23 Eq 5.4.2, 5.5.2]; PLAN R1; numerics N1.
 * Oracle: the defining equation's residual, on the grid
 * k in {1, 2, 3} x p in {0, 1/2, 1, 2} x b in {0, 0.25, 1}, with a real
 * (a = 1.3) and complex (a = 1 + 0.5i, b scaled by 1 - 0.4i), at t = 0.5, 1, 2;
 * plus k = 1/2 with real positive a, b. */
/* Verifies: FR-MOP-001, NFR-NUM-001 (T-MOP-01) */
static void test_mop_first_residual(void) {
    static const double ps[] = { 0.0, 0.5, 1.0, 2.0 }, bs[] = { 0.0, 0.25, 1.0 };
    static const double ts[] = { 0.5, 1.0, 2.0 };
    double worst = 0.0, worst_half = 0.0;
    int    k, ip, ib, it, cx;

    printf("\n[T-MOP-01] First Equation: residual of (d~/dt)^k alpha = beta\n");
    for (k = 1; k <= 3; k++)
        for (ip = 0; ip < 4; ip++)
            for (ib = 0; ib < 3; ib++)
                for (cx = 0; cx < 2; cx++) {
                    double complex a = cx ? 1.0 + 0.5 * I : 1.3;
                    double complex b = cx ? bs[ib] * (1.0 - 0.4 * I) : bs[ib];
                    gia_beta be = affine(a, b, ps[ip]);
                    gia_rational kk = { k, 1 };
                    for (it = 0; it < 3; it++) {
                        double r = first_eq_residual(&be, kk, ts[it]);
                        if (!(r <= worst)) worst = r;
                    }
                }
    printf("    worst relative residual over 216 cases: %.3g\n", worst);
    ok("k = 1,2,3 x p = 0,1/2,1,2 x b = 0,.25,1, real and complex: <= 1e-6",
       worst <= TOL_RESIDUAL);
    for (ip = 0; ip < 4; ip++)
        for (ib = 0; ib < 3; ib++) {
            gia_beta be = affine(1.3, bs[ib], ps[ip]);
            gia_rational half = { 1, 2 };
            for (it = 0; it < 3; it++) {
                double r = first_eq_residual(&be, half, ts[it]);
                if (!(r <= worst_half)) worst_half = r;
            }
        }
    ok("k = 1/2, real positive a, b: <= 1e-6", worst_half <= TOL_RESIDUAL);
}

/* The printed forms of [23 Eq 5.5.7] and [23 Eq 5.5.8], as
 * probes/eq_5_5_8_residual.py writes them. */
static double printed_557(double a, double b, double p, double k, double t) {
    double I_ = (pow(a + b * t, (p + k) / k) - pow(a, (p + k) / k)) * k / (b * (p + k));
    return (1.0 / k) * pow(I_, k);
}
static double printed_558(double a, double b, double p, double k, double t) {
    return (1.0 / (k * b)) * pow((k / (p + k)) * pow(a + b * t, p / k + 1.0), k);
}
static double printed_residual(double (*f)(double, double, double, double, double), double t) {
    const double a = 1.0, b = 0.25, p = 1.0, k = 2.0, h = 1e-6;
    double d = (f(a, b, p, k, t + h) - f(a, b, p, k, t - h)) / (2.0 * h), v = f(a, b, p, k, t);
    return pow(d / v, k) * v - pow(a + b * t, p);
}

/* Source: [23 Eq 5.5.7-5.5.8]; PLAN X1; probes/eq_5_5_8_residual.py.
 * The printed solutions leave residuals 1.25 and -0.625 at (1, 0.25, 1, 2),
 * t = 1, where the derived one leaves 0. */
/* Verifies: FR-MOP-001 (T-MOP-02) */
static void test_mop_x1(void) {
    double r7 = printed_residual(printed_557, 1.0), r8 = printed_residual(printed_558, 1.0);
    gia_beta be = affine(1.0, 0.25, 1.0);
    gia_rational k2 = { 2, 1 };

    printf("\n[T-MOP-02] X1: the printed [23 Eq 5.5.7-5.5.8] do not solve Eq 5.5.2\n");
    printf("    printed 5.5.7: %.6g   printed 5.5.8: %.6g   derived: %.3g\n", r7, r8,
           first_eq_residual(&be, k2, 1.0));
    ok("printed 5.5.7 residual = 1.25 (> 100x tol)", fabs(r7 - 1.25) < 1e-4 && fabs(r7) > 100 * TOL_RESIDUAL);
    ok("printed 5.5.8 residual = -0.625 (> 100x tol)", fabs(r8 + 0.625) < 1e-4 && fabs(r8) > 100 * TOL_RESIDUAL);
    ok("the derived solution's residual is <= 1e-6", first_eq_residual(&be, k2, 1.0) <= TOL_RESIDUAL);
}

/* Source: PLAN R1 ([23 Eq 5.5.5]: lower limit 0); numerics N1.
 * alpha(0) = 0 exactly; as b -> 0 the solution tends to the b = 0 one,
 * (beta^{1/k} t / k)^k, continuously and to 1e-12 (no separate naive branch). */
/* Verifies: FR-MOP-001, NFR-NUM-002 (T-MOP-03) */
static void test_mop_b_to_zero(void) {
    double complex al;
    gia_rational   k2 = { 2, 1 };
    gia_beta       be;
    double         b;
    int            all = 1;
    /* a = 1, p = 1, k = 2, t = 1: alpha(b = 0) = (1 * 1 / 2)^2 = 0.25. */
    printf("\n[T-MOP-03] alpha(0) = 0; continuity through b = 0\n");
    be = affine(1.0, 0.25, 1.0);
    ok("alpha(0) = 0 exactly", gia_mop_couple(&be, k2, 0.0, &al, NULL) == GIA_OK && al == 0.0);
    be = affine(1.0, 0.0, 1.0);
    ok("b = 0: alpha = (beta^{1/k} t / k)^k = 0.25",
       gia_mop_couple(&be, k2, 1.0, &al, NULL) == GIA_OK && near_c(al, 0.25, TOL_CLOSED));
    for (b = 1e-6; b >= 1e-14; b /= 100.0) {
        /* exact: [((1+b)^{3/2} - 1)/(3b/... )]^2 -> 0.25 (1 + b/2 + ...)^2; the
         * first-order term is 0.25 * (1 + 3b/4 ... ); within 1e-12 once b <= 1e-12 */
        be = affine(1.0, b, 1.0);
        if (gia_mop_couple(&be, k2, 1.0, &al, NULL) != GIA_OK ||
            !(cabs(al - 0.25) <= 0.25 * (b + 1e-14))) all = 0;
    }
    ok("b = 1e-6 .. 1e-14: alpha -> 0.25 with error O(b), no jump", all);
}

/* Source: [23 Eq 5.6.1]. */
/* Verifies: FR-MOP-003 (T-MOP-04) */
static void test_mop_matrioska(void) {
    gia_beta      be[16];
    gia_matrioska M;
    gia_rational  k1 = { 1, 1 };
    const char   *why = NULL;
    int           i, j, related = 0, diag = 1, indep = 1;

    printf("\n[T-MOP-04] the Matrioska: zero diagonal, N(N-1) independent couples\n");
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            be[i * 4 + j] = affine(1.0 + i, 0.1 * (j + 1), 1.0);
    be[1 * 4 + 3].kind = GIA_BETA_NONE;                 /* one unrelated couple */
    memset(&M, 0, sizeof M);
    ok("solves N = 4", gia_mop_solve(4, be, k1, 1.0, &M, &why) == GIA_OK && M.N == 4);
    if (!M.a || !M.related) {
        ok("diagonal, couples and independence (not reached)", 0);
        return;
    }
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++) {
            double complex want;
            if (i == j) { if (M.related[i * 4 + j] || M.a[i * 4 + j] != 0.0) diag = 0; continue; }
            if (!M.related[i * 4 + j]) continue;
            related++;
            (void)gia_mop_couple(&be[i * 4 + j], k1, 1.0, &want, NULL);
            if (M.a[i * 4 + j] != want) indep = 0;
        }
    ok("diagonal: unrelated, 0", diag);
    ok("N(N-1) - 1 = 11 related couples; the NONE couple is unrelated, not 0",
       related == 11 && !M.related[1 * 4 + 3]);
    ok("each entry is its couple solved alone", indep);
    gia_matrioska_free(&M);
    be[0 * 4 + 1] = affine(1.0, -1.0, 1.0);             /* a + b t = 0 at t = 1 */
    ok("a failing couple: the status, nothing allocated",
       gia_mop_solve(4, be, k1, 2.0, &M, &why) == GIA_E_DOMAIN);
}

/* Source: PLAN R1; numerics N1, §1. */
/* Verifies: FR-MOP-002, NFR-NUM-004 (T-MOP-05) */
static void test_mop_domain(void) {
    double complex al = SENTINEL_C;
    const char    *why;
    gia_beta       be;
    gia_rational   half = { 1, 2 }, k1 = { 1, 1 }, k0 = { 0, 1 }, kneg = { -1, 1 }, kbad = { 1, 0 };
    int            untouched = 1;

    printf("\n[T-MOP-05] the First Equation's domain\n");
    why = NULL; be = affine(1.0 + 0.5 * I, 0.25, 1.0);
    ok("non-integer k with complex a: GIA_E_DOMAIN, with a reason",
       gia_mop_couple(&be, half, 1.0, &al, &why) == GIA_E_DOMAIN && why);
    why = NULL; be = affine(1.0, -0.25, 1.0);
    ok("non-integer k with b < 0: GIA_E_DOMAIN",
       gia_mop_couple(&be, half, 1.0, &al, &why) == GIA_E_DOMAIN);
    why = NULL; be = affine(1.0, -1.0, 1.0);
    ok("a + b t = 0 on [0, t] (t* = 1 < 2): GIA_E_DOMAIN",
       gia_mop_couple(&be, k1, 2.0, &al, &why) == GIA_E_DOMAIN && why);
    ok("  but t = 0.5, before the zero, solves", gia_mop_couple(&be, k1, 0.5, &al, &why) == GIA_OK);
    al = SENTINEL_C; why = NULL; be = affine(1.0, 0.25, -1.0);
    ok("p = -k with b != 0: GIA_E_DOMAIN",
       gia_mop_couple(&be, k1, 1.0, &al, &why) == GIA_E_DOMAIN && why);
    if (al != SENTINEL_C) untouched = 0;
    ok("k = 0, k < 0, den = 0: refused",
       gia_mop_couple(&be, k0, 1.0, &al, &why) == GIA_E_DOMAIN &&
       gia_mop_couple(&be, kneg, 1.0, &al, &why) == GIA_E_DOMAIN &&
       gia_mop_couple(&be, kbad, 1.0, &al, &why) == GIA_E_ARG);
    ok("outputs untouched on refusal", untouched);
    ok("t < 0 is GIA_E_ARG", gia_mop_couple(&be, k1, -1.0, &al, &why) == GIA_E_ARG);
}

/* Source: [23 Eq 5.5.5-5.5.7]; PLAN R1; numerics N2. */
/* Verifies: FR-MOP-004, NFR-NUM-004 (T-MOP-06) */
static void test_mop_samples(void) {
    static const double         ta[3] = { 0.0, 1.0, 2.0 };
    static const double complex va[3] = { 1.0 + 0.5 * I, 1.25 + 0.4 * I, 1.5 + 0.3 * I };
    /* Crosses the negative real axis between t = 1 and 2 (at -1 + 0i). */
    static const double         tc[3] = { 0.0, 1.0, 2.0 };
    static const double complex vc[3] = { 1.0 + 0.5 * I, -1.0 + 0.5 * I, -1.0 - 0.5 * I };
    gia_beta      s, aff;
    gia_rational  k1 = { 1, 1 }, k2 = { 2, 1 };
    double complex x, y, lo, hi;
    const char   *why = NULL;
    int           all = 1, kk;
    double        t;

    printf("\n[T-MOP-06] sampled boundary conditions\n");
    memset(&s, 0, sizeof s);
    s.kind = GIA_BETA_SAMPLES; s.n = 3; s.t = ta; s.v = va;
    /* Affine through the same points: beta = (1 + 0.5i) + (0.25 - 0.1i) t, p = 1. */
    aff = affine(1.0 + 0.5 * I, 0.25 - 0.1 * I, 1.0);
    for (kk = 1; kk <= 2; kk++)
        for (t = 0.0; t <= 2.0 + 1e-12; t += 0.25) {
            gia_rational k = kk == 1 ? k1 : k2;
            if (gia_mop_couple(&s, k, t, &x, &why) != GIA_OK ||
                gia_mop_couple(&aff, k, t, &y, &why) != GIA_OK ||
                !near_c(x, y, TOL_QUAD)) all = 0;
        }
    ok("samples of an affine beta reproduce N1 to 1e-9 (k = 1, 2)", all);

    s.t = tc; s.v = vc;
    ok("a path across the negative real axis: solves at t = 2",
       gia_mop_couple(&s, k2, 2.0, &x, &why) == GIA_OK);
    (void)gia_mop_couple(&s, k2, 1.5 - 1e-7, &lo, &why);
    (void)gia_mop_couple(&s, k2, 1.5 + 1e-7, &hi, &why);
    ok("  continuously: no jump where it crosses (|da| ~ dt)", cabs(hi - lo) < 1e-5);
    {
        /* Independent oracle: the test unwraps the phase itself, step by step
         * (|d theta| < pi per step), and integrates |beta|^{1/2} e^{i theta/2}
         * by the composite midpoint rule, 400,000 steps. A principal root taken
         * per point flips the integrand's sign past the crossing. */
        const int      steps = 400000;
        double         h = 2.0 / steps, th = carg(vc[0]), prev_arg = carg(vc[0]);
        double complex acc = 0.0, want;
        int            q;
        for (q = 0; q < steps; q++) {
            double         tm = (q + 0.5) * h, w = tm < 1.0 ? tm : tm - 1.0;
            double complex bt = tm < 1.0 ? vc[0] + w * (vc[1] - vc[0]) : vc[1] + w * (vc[2] - vc[1]);
            double         a = carg(bt), d = a - prev_arg;
            if (d > M_PI) d -= 2.0 * M_PI;
            if (d < -M_PI) d += 2.0 * M_PI;
            th += d; prev_arg = a;
            acc += sqrt(cabs(bt)) * cexp(I * th / 2.0) * h;
        }
        want = (acc / 2.0) * (acc / 2.0);         /* alpha = ((1/k) int)^k, k = 2 */
        ok("  and agrees with an independently unwrapped integral (1e-6)",
           near_c(x, want, 1e-6));
    }
    ok("t past the last sample: GIA_E_DOMAIN", gia_mop_couple(&s, k1, 2.5, &x, &why) == GIA_E_DOMAIN);
    {
        static const double         tz[2] = { 0.0, 1.0 };
        static const double complex vz[2] = { 1.0, -1.0 };       /* through 0 at t = 1/2 */
        s.t = tz; s.v = vz; s.n = 2;
        ok("a segment through beta = 0: GIA_E_DOMAIN",
           gia_mop_couple(&s, k1, 1.0, &x, &why) == GIA_E_DOMAIN);
    }
}

/* E(x)'s series, long double: E(x) = sum_{n>=0} binom(q, n+1) x^n / q. */
static long double e_series(long double q, long double x) {
    long double term = q, sum = 0.0L;    /* binom(q, 1) = q */
    int n;
    for (n = 0; n < 60; n++) {
        sum += term * powl(x, (long double)n);
        term *= (q - (long double)(n + 1)) / (long double)(n + 2);
    }
    return sum / q;
}

/* Source: numerics.md N1; probes/first_equation_numerics.py.
 * a = 1, b = 0.25, p = 1, k = 2, so q = 3/2 and alpha = (t/2 E(x))^2, x = bt.
 * The naive closed form loses 5.9e-5 at x = 1e-11; the unified form keeps
 * 1e-12 against the long-double series. */
/* Verifies: NFR-NUM-002, FR-MOP-001 (T-NUM-01) */
static void test_num_cancellation(void) {
    static const double xs[] = { 1e-2, 1e-5, 1e-8, 1e-11 };
    gia_beta     be = affine(1.0, 0.25, 1.0);
    gia_rational k2 = { 2, 1 };
    size_t       i;
    int          all = 1;

    printf("\n[T-NUM-01] no cancellation near t = 0\n");
    for (i = 0; i < 4; i++) {
        double         t = xs[i] / 0.25;
        long double    S = (long double)t / 2.0L * e_series(1.5L, (long double)xs[i]);
        double         want = (double)(S * S);
        double complex al;
        if (gia_mop_couple(&be, k2, t, &al, NULL) != GIA_OK ||
            !(fabs(creal(al) - want) <= TOL_CLOSED * want)) {
            printf("    x = %g: got %.17g want %.17g\n", xs[i], creal(al), want);
            all = 0;
        }
    }
    ok("relative error <= 1e-12 for x = 1e-2, 1e-5, 1e-8, 1e-11", all);
}

/* Verifies: NFR-NUM-003 (T-NUM-02) */
static void test_num_overflow(void) {
    gia_beta       be = affine(1e200, 0.0, 2.0);
    gia_rational   k1 = { 1, 1 }, k3 = { 3, 1 };
    double complex al = SENTINEL_C;
    const char    *why = NULL;
    printf("\n[T-NUM-02] never a non-finite value\n");
    ok("alpha = beta t = 1e400 overflows: GIA_E_RANGE",
       gia_mop_couple(&be, k1, 1.0, &al, &why) == GIA_E_RANGE && al == SENTINEL_C);
    be = affine(1e300, 0.0, 3.0);         /* S = 1e300 t/3 finite, S^3 is not */
    ok("S^k past DBL_MAX (k = 3): GIA_E_RANGE",
       gia_mop_couple(&be, k3, 1.0, &al, &why) == GIA_E_RANGE);
}

static double complex q_sqrt_recip(double t, void *c) { (void)c; return 1.0 / sqrt(t); }
static double complex q_smooth(double t, void *c) { (void)c; return cexp(I * 3.0 * t) * (1.0 + t * t); }

/* Source: numerics N2; NFR-NUM-005. */
/* Verifies: NFR-NUM-005, NFR-NUM-001 (T-NUM-03) */
static void test_num_quadrature(void) {
    double complex I_, exact;
    double         err = -1.0;
    const char    *why = NULL;
    printf("\n[T-NUM-03] the quadrature reports its error, and refuses\n");
    /* int_0^2 e^{3it} (1 + t^2) dt by parts. */
    {
        double complex w = 3.0 * I, e2 = cexp(2.0 * w);
        exact = (e2 - 1.0) / w + (e2 * (4.0 / w - 4.0 / (w * w) + 2.0 / (w * w * w)) - 2.0 / (w * w * w));
    }
    ok("smooth integrand: converges to 1e-10",
       gia_quad_gk15(q_smooth, NULL, 0.0, 2.0, 1e-10, 50, &I_, &err, &why) == GIA_OK &&
       near_c(I_, exact, TOL_QUAD));
    ok("  and its error estimate bounds the true error", err >= cabs(I_ - exact));
    ok("an integrable singularity (1/sqrt t at 0) past 50 levels: GIA_E_CONVERGENCE",
       gia_quad_gk15(q_sqrt_recip, NULL, 0.0, 1.0, 1e-14, 50, &I_, &err, &why) == GIA_E_CONVERGENCE);
}

/* Source: PLAN §1 B6, G6. */
/* Verifies: NFR-NUM-006 (T-NUM-04) */
static void test_num_no_clamp(void) {
    gia_beta       be = affine(-2.0, 0.5, 1.0);
    gia_rational   k1 = { 1, 1 }, k3 = { 3, 1 };
    double complex al;
    printf("\n[T-NUM-04] negative and complex coordinates survive\n");
    /* k = 1: alpha = int_0^t (-2 + 0.5 s) ds = -2t + t^2/4. */
    ok("negative beta: alpha = -2t + t^2/4 = -1.75 at t = 1",
       gia_mop_couple(&be, k1, 1.0, &al, NULL) == GIA_OK && near_c(al, -1.75, TOL_CLOSED));
    /* k = 3, beta = -1: alpha = (beta^{1/3} t/3)^3 = -t^3/27 (principal cube
     * root e^{i pi/3}, cubed: e^{i pi} = -1). */
    be = affine(-1.0, 0.0, 1.0);
    ok("k = 3, beta = -1: alpha = -1/27 at t = 1, sign intact",
       gia_mop_couple(&be, k3, 1.0, &al, NULL) == GIA_OK && near_c(al, -1.0 / 27.0, TOL_CLOSED));
    be = affine(0.5 - 2.0 * I, 0.0, 1.0);
    ok("complex beta, k = 1: alpha = beta t", gia_mop_couple(&be, k1, 2.0, &al, NULL) == GIA_OK &&
                                               near_c(al, 1.0 - 4.0 * I, TOL_CLOSED));
}

/* ------------------------------------------------------------------ *
 * 7.4 Relational algebra
 * ------------------------------------------------------------------ */

static const rel_t RI = { 1, 0, 0 }, RJ = { 0, 1, 0 }, RK = { 0, 0, 1 };

static int rel_eq(rel_t a, rel_t b, double tol) {
    return fabs(a.i - b.i) <= tol && fabs(a.j - b.j) <= tol && fabs(a.k - b.k) <= tol;
}
static rel_t R(double i, double j, double k) { rel_t r; r.i = i; r.j = j; r.k = k; return r; }

/* Source: [23 Eq 5.1.3-5.1.5] (page image p. 3183); PLAN R2. */
/* Verifies: FR-REL-001 (T-REL-01) */
static void test_rel_table(void) {
    rel_t u[3] = { RI, RJ, RK };
    /* The printed table, row x column. */
    rel_t want[3][3] = { { RI, RJ, RK }, { RJ, { -1, 0, 0 }, RK }, { RK, RK, { -1, 0, 0 } } };
    int   a, b, all = 1, bilinear;
    printf("\n[T-REL-01] the relational product, as printed\n");
    for (a = 0; a < 3; a++)
        for (b = 0; b < 3; b++)
            if (!rel_eq(rel_mul(u[a], u[b]), want[a][b], 0.0)) all = 0;
    ok("all nine products of [23 Eq 5.1.3-5.1.5]", all);
    /* Bilinear: (2i - j + 3k) o (0.5i + 4j - k), expanded by hand from the
     * table: i:  2*0.5 - (-1)(4)(-1)... */
    {
        rel_t x = R(2, -1, 3), y = R(0.5, 4, -1), p = rel_mul(x, y);
        /* i: x_i y_i - x_j y_j - x_k y_k = 1 + 4 + 3 = 8
         * j: x_i y_j + x_j y_i = 8 - 0.5 = 7.5
         * k: x_i y_k + x_j y_k + x_k y_i + x_k y_j = -2 + 1 + 1.5 + 12 = 12.5 */
        bilinear = rel_eq(p, R(8, 7.5, 12.5), 1e-15);
    }
    ok("bilinear: (2i - j + 3k) o (0.5i + 4j - k) = 8i + 7.5j + 12.5k by hand", bilinear);
    ok("commutative: x o y = y o x",
       rel_eq(rel_mul(R(2, -1, 3), R(0.5, 4, -1)), rel_mul(R(0.5, 4, -1), R(2, -1, 3)), 0.0));
}

/* Source: [23 Eq 5.1.3-5.1.5]; PLAN R2. */
/* Verifies: FR-REL-004 (T-REL-02) */
static void test_rel_left_to_right(void) {
    printf("\n[T-REL-02] no reassociation\n");
    ok("rel_mul3(j, j, k) = (j o j) o k = -k", rel_eq(rel_mul3(RJ, RJ, RK), R(0, 0, -1), 0.0));
    ok("j o (j o k) = +k: the order is observable", rel_eq(rel_mul(RJ, rel_mul(RJ, RK)), RK, 0.0));
}

/* Source: [23 Eq 5.1.2]; PLAN R2; numerics N8. */
/* Verifies: FR-REL-002 (T-REL-03) */
static void test_rel_exp(void) {
    rel_t       e;
    const char *why = NULL;
    double      a = 0.3, b = 1.2, c = -0.5, rho = sqrt(b * b + c * c), ea = exp(a);
    int         cont = 1, k;
    printf("\n[T-REL-03] the De Moivre exponential\n");
    ok("Exp{a i + b j + c k} = e^a [cos rho, b sin rho/rho, c sin rho/rho]",
       rel_exp(R(a, b, c), &e, &why) == GIA_OK &&
       rel_eq(e, R(ea * cos(rho), ea * b * sin(rho) / rho, ea * c * sin(rho) / rho), 1e-15));
    ok("rho = 0: e^a", rel_exp(R(a, 0, 0), &e, &why) == GIA_OK && rel_eq(e, R(ea, 0, 0), 0.0));
    for (k = 3; k <= 9; k++) {
        double r = pow(10.0, -k);
        if (rel_exp(R(a, r, 0), &e, &why) != GIA_OK ||
            !rel_eq(e, R(ea * cos(r), ea * r * (1.0 - r * r / 6.0), 0), 1e-15)) cont = 0;
    }
    ok("rho = 1e-3 .. 1e-9: continuous into the limit, series sin rho/rho", cont);
    /* Not a power series: Exp(j pi/2) = j, whereas sum (j pi/2)^n/n! under the
     * non-associative table depends on the bracketing. */
    ok("Exp(j pi/2) = j", rel_exp(R(0, M_PI / 2.0, 0), &e, &why) == GIA_OK &&
                          rel_eq(e, RJ, 1e-15));
    ok("e^a past DBL_MAX: GIA_E_RANGE", rel_exp(R(800, 0, 0), &e, &why) == GIA_E_RANGE);
}

/* Source: [23 Eq A2.5-A2.6]; PLAN R3, X10; probes/ordinal_root_power.py. */
/* Verifies: FR-REL-003 (T-REL-04) */
static void test_rel_roots(void) {
    rel_t       r, p;
    const char *why = NULL;
    int         N, l, unity = 1;
    printf("\n[T-REL-04] ordinal roots: angle powers, and the table power apart\n");
    for (N = 3; N <= 7; N++)
        for (l = 1; l < N; l++)
            if (rel_root_pow(N, l, N - 1, &p, &why) != GIA_OK || !rel_eq(p, RI, 1e-12)) unity = 0;
    ok("rel_root_pow(N, l, N - 1) = 1 for N = 3..7, every l", unity);
    ok("rel_root(4, 1) = (cos 2pi/3, sin(2pi/3)/sqrt2, sin(2pi/3)/sqrt2)",
       rel_root(4, 1, &r, &why) == GIA_OK &&
       rel_eq(r, R(cos(2 * M_PI / 3), sin(2 * M_PI / 3) / sqrt(2.0), sin(2 * M_PI / 3) / sqrt(2.0)), 1e-15));
    p = rel_mul_pow(r, 3);
    printf("    table cube of r(4,1): (%.6f, %.6f, %.6f)\n", p.i, p.j, p.k);
    ok("X10: rel_mul_pow(r(4,1), 3) = (0.540721, 0, -0.665721), not 1",
       rel_eq(p, R(0.540721, 0, -0.665721), 1e-6));
    ok("rel_mul_pow(x, 0) = i; rel_root_pow(N, l, 0) = i",
       rel_eq(rel_mul_pow(r, 0), RI, 0.0) && rel_root_pow(4, 1, 0, &p, &why) == GIA_OK &&
       rel_eq(p, RI, 1e-15));
    ok("N < 2, l outside 0..N-1, m < 0: GIA_E_ARG",
       rel_root(1, 0, &r, &why) == GIA_E_ARG && rel_root(4, 4, &r, &why) == GIA_E_ARG &&
       rel_root_pow(4, 1, -1, &p, &why) == GIA_E_ARG);
}

/* Source: PLAN §6 (no division or non-integer power in the relational
 * algebra). k = 1: alpha = int beta componentwise, by hand. */
/* Verifies: FR-MOP-007 (T-MOP-09) */
static void test_mop_relational_couple(void) {
    gia_beta     be[3];
    gia_rational k1 = { 1, 1 }, k2 = { 2, 1 };
    rel_t        al = { 99, 99, 99 };
    const char  *why = NULL;
    printf("\n[T-MOP-09] relational-valued couples\n");
    be[0] = affine(2.0, 0.5, 1.0);      /* int_0^t (2 + s/2) = 2t + t^2/4  */
    be[1] = affine(-1.0, 0.0, 1.0);     /* -t                               */
    be[2] = affine(1.0, 1.0, 2.0);      /* int (1+s)^2 = ((1+t)^3 - 1)/3   */
    ok("k = 1: componentwise, by hand at t = 2",
       gia_mop_couple_rel(be, k1, 2.0, &al, &why) == GIA_OK &&
       rel_eq(al, R(5.0, -2.0, 26.0 / 3.0), 1e-12));
    al = R(99, 99, 99); why = NULL;
    ok("k = 2: GIA_E_UNSUPPORTED, outputs untouched, reason given",
       gia_mop_couple_rel(be, k2, 2.0, &al, &why) == GIA_E_UNSUPPORTED &&
       rel_eq(al, R(99, 99, 99), 0.0) && why != NULL);
}

/* ------------------------------------------------------------------ *
 * The EQS (FR-MOP-006)
 * ------------------------------------------------------------------ */

/* [23 Eq 7.1-7.5] written out as printed, bracket by bracket, for a test to
 * compare against: no relational product involved. */
static void eqs_by_hand(const gia_eqs_params *p, double S0, double F0, double T0, int l,
                        double out[3]) {
    double r2psi = p->psi2 * (p->eps[1] + 2.0 * M_PI * l) / (p->N - 1);
    double B = cos(r2psi), C = sin(r2psi) / sqrt(2.0);
    double E1 = (p->eps[0] + 4.0 * M_PI * l) / (p->N - 1);
    double E2 = (p->eps[1] + 4.0 * M_PI * l) / (p->N - 1);
    double E3 = (p->eps[2] + 4.0 * M_PI * l) / (p->N - 1);
    out[0] = p->A * exp(p->psi1[0] * E1 * (B * S0 - C * (F0 + T0)));
    out[1] = p->psi1[1] * E2 * (B * F0 + C * S0);
    out[2] = p->psi1[2] * E3 * (B * T0 + C * S0 + C * (F0 + T0));
}

/* Source: [23 Eq 7.1-7.5]; PLAN R2, X11.
 * N = 4 with given coordinates and factors, every l, by hand. */
/* Verifies: FR-MOP-006 (T-MOP-08) */
static void test_mop_eqs(void) {
    gia_eqs_params p = { { 0.7, 1.1, -0.4 }, 0.9, { 0.2, 0.3, 0.3 }, 1.5, 4 };
    rel_t          ref = { 0.6, -0.25, 0.8 };
    double         got[3], want[3];
    const char    *why = NULL;
    int            l, all = 1;

    printf("\n[T-MOP-08] the EQS operative form, N = 4\n");
    for (l = 1; l <= 3; l++) {
        eqs_by_hand(&p, ref.i, ref.j, ref.k, l, want);
        if (gia_eqs(&p, ref, l, got, &why) != GIA_OK ||
            !near_c(got[0], want[0], TOL_CLOSED) || !near_c(got[1], want[1], TOL_CLOSED) ||
            !near_c(got[2], want[2], TOL_CLOSED)) {
            printf("    l = %d: got %.15g %.15g %.15g want %.15g %.15g %.15g\n", l,
                   got[0], got[1], got[2], want[0], want[1], want[2]);
            all = 0;
        }
    }
    ok("rho, phi, theta = [23 Eq 7.1-7.5] by hand, l = 1, 2, 3, 1e-12", all);
    p.eps[2] = 0.31; why = NULL;
    ok("X11: eps_2 != eps_3 refused, GIA_E_DOMAIN with a reason",
       gia_eqs(&p, ref, 1, got, &why) == GIA_E_DOMAIN && why && strstr(why, "X11"));
    p.eps[2] = 0.3;
    ok("l outside 1..N-1, N < 3: GIA_E_ARG",
       gia_eqs(&p, ref, 0, got, &why) == GIA_E_ARG && gia_eqs(&p, ref, 4, got, &why) == GIA_E_ARG);
    p.A = 1.0; p.psi1[0] = 800.0; ref.i = -10.0;   /* B < 0 at l = 1, so S > 0 */
    ok("rho overflow: GIA_E_RANGE", gia_eqs(&p, ref, 1, got, &why) == GIA_E_RANGE);
}

/* A fixed 64-bit LCG (vv-plan.md §6), mapped to [lo, hi). No global RNG. */
static double lcg_uniform(unsigned long long *x, double lo, double hi) {
    *x = *x * 6364136223846793005ULL + 1442695040888963407ULL;
    return lo + (hi - lo) * (double)(*x >> 11) / 9007199254740992.0;
}

/* Source: [23 Eq 7.1.1, 7.2, 7.3]; probes/eqs_relational_product.py.
 * Validation: the source's own brackets are the literal table's product
 * of the De Moivre root and the reference coordinates, over 1000 draws. */
/* Verifies: BR-005, BR-008, FR-MOP-006, FR-REL-001 (VAL-02) */
static void test_val_eqs_brackets(void) {
    unsigned long long x = 1;
    double worst = 0.0;
    int    n;
    printf("\n[VAL-02] [23 Eq 7.1-7.3]'s brackets are rel_mul(root, ref)\n");
    for (n = 0; n < 1000; n++) {
        double psi = lcg_uniform(&x, -3, 3), S = lcg_uniform(&x, -2, 2);
        double F = lcg_uniform(&x, -2, 2), T = lcg_uniform(&x, -2, 2);
        double r = sqrt(2.0) * psi, B = cos(r), C = sin(r) / sqrt(2.0);
        rel_t  got = rel_mul(R(B, C, C), R(S, F, T));
        double d0 = fabs(got.i - (B * S - C * (F + T)));
        double d1 = fabs(got.j - (B * F + C * S));
        double d2 = fabs(got.k - (B * T + C * S + C * (F + T)));
        if (d0 > worst) worst = d0;
        if (d1 > worst) worst = d1;
        if (d2 > worst) worst = d2;
    }
    printf("    worst |table product - printed bracket| over 1000 draws: %.3g\n", worst);
    ok("1000 seeded draws: <= 1e-12", worst <= TOL_CLOSED);
}

/* Source: [23 Eq 6.1-6.3]; PLAN R4 (reconstructed oracle), R12.
 * u = A' must solve u' + u^2 = 0 (derivatives by central differences of the
 * returned A); B is specular; {c2, t} is c2 t: A(t) - A(0) = ln(1 + c2 t/c1);
 * c1 + c2 t <= 0 is refused. The mutation ln(c1 + c2 + t) fails the
 * reduction and the Riccati residual. */
/* Verifies: FR-MOP-005 (T-MOP-07) */
static void test_mop_second(void) {
    const double complex a0 = 0.3 + 0.1 * I;
    const double         c1 = 1.0, c2 = 0.5;
    gia_second           s0, sp, sm, s;
    gia_matrioska        r;
    const char          *why = NULL;
    double               t;
    int                  ric = 1, spec = 1, red = 1, j;

    printf("\n[T-MOP-07] the Second Equation's printed solution\n");
    for (t = 0.25; t <= 2.0 + 1e-12; t += 0.25) {
        const double h = 1e-3;
        double complex u, up;
        if (gia_mop_second(a0, c1, c2, 4, t, &s, NULL, &why) != GIA_OK ||
            gia_mop_second(a0, c1, c2, 4, t + h, &sp, NULL, &why) != GIA_OK ||
            gia_mop_second(a0, c1, c2, 4, t - h, &sm, NULL, &why) != GIA_OK ||
            gia_mop_second(a0, c1, c2, 4, 0.0, &s0, NULL, &why) != GIA_OK) {
            ric = spec = red = 0; continue;
        }
        u  = (sp.A - sm.A) / (2.0 * h);
        up = (sp.A - 2.0 * s.A + sm.A) / (h * h);
        if (!(cabs(up + u * u) <= TOL_RESIDUAL * cabs(u * u) * 10.0)) ric = 0;
        if (s.B[0][0] != s.A || s.B[1][1] != s.A || s.B[0][1] != -s.A || s.B[1][0] != -s.A)
            spec = 0;
        if (!near_c(s.A - s0.A, log(1.0 + c2 * t / c1), TOL_CLOSED)) red = 0;
    }
    ok("u = A' solves u' + u^2 = 0 (central differences, 1e-5)", ric);
    ok("B = [[A, -A], [-A, A]], exactly", spec);
    ok("{c2, t} reduces to c2 t: A(t) - A(0) = ln(1 + c2 t / c1)", red);
    ok("A(0) = alpha12(0) w, w = e^{2 pi i/(N - 1)}",
       gia_mop_second(a0, c1, c2, 4, 0.0, &s, NULL, &why) == GIA_OK &&
       near_c(s.A, a0 * cexp(2.0 * M_PI * I / 3.0) + log(c1), TOL_CLOSED));
    {
        const int      solved = gia_mop_second(a0, c1, c2, 5, 1.0, &s, &r, &why) == GIA_OK;
        double complex eB11   = solved ? 1.0 + (cexp(2.0 * s.A) - 1.0) / 2.0 : 0.0;
        int            row = solved, rest = solved;
        ok("solves with a Matrioska row", solved);
        if (solved) for (j = 1; j < 5; j++)
            if (!r.related[j] || !near_c(r.a[j], eB11 * cexp(2.0 * M_PI * I * (j - 1) / 4.0), TOL_CLOSED))
                row = 0;
        if (solved) for (j = 5; j < 25; j++) if (r.related[j]) rest = 0;
        ok("r_1j = (e^B)_11 w^{j-2}, e^B exact: I + (e^{2A} - 1)/2 M", row);
        ok("only row 1 is related", rest && !r.related[0]);
        if (solved) gia_matrioska_free(&r);
    }
    why = NULL;
    ok("c1 + c2 t <= 0 on [0, t]: GIA_E_DOMAIN",
       gia_mop_second(a0, 1.0, -0.5, 4, 3.0, &s, NULL, &why) == GIA_E_DOMAIN && why);
    ok("c1 <= 0: GIA_E_DOMAIN; N < 2: GIA_E_ARG",
       gia_mop_second(a0, 0.0, 1.0, 4, 1.0, &s, NULL, &why) == GIA_E_DOMAIN &&
       gia_mop_second(a0, 1.0, 1.0, 1, 1.0, &s, NULL, &why) == GIA_E_ARG);
    ok("e^{2A} past DBL_MAX (alpha12(0) = 400, N = 2): GIA_E_RANGE",
       gia_mop_second(400.0, 1.0, 1.0, 2, 1.0, &s, NULL, &why) == GIA_E_RANGE);
}

int main(void) {
    printf("=== Giannantoni kernel: verification and validation ===\n");
    test_status_contract();
    test_idc_general_f();
    test_lde2();
    test_lde2_zero_root();
    test_lde2_linear_in_ics();
    test_binary();
    test_binary_linear_in_ics();
    test_riccati();
    test_idc_refusals();
    test_idc_taylor();
    test_val_zoli_2009();
    test_nl1410();
    test_solution_drift();
    test_drift_projection();
    test_em_coproduction();
    test_em_interaction();
    test_em_global_balance();
    test_val_emergy_rules();
    test_em_limits();
    test_em_ordinal_forms();
    test_em_circle_product();
    test_mop_first_residual();
    test_mop_x1();
    test_mop_b_to_zero();
    test_mop_matrioska();
    test_mop_domain();
    test_mop_samples();
    test_num_cancellation();
    test_num_overflow();
    test_num_quadrature();
    test_num_no_clamp();
    test_rel_table();
    test_rel_left_to_right();
    test_rel_exp();
    test_rel_roots();
    test_mop_relational_couple();
    test_mop_eqs();
    test_val_eqs_brackets();
    test_mop_second();

    printf("\n%s\nfailures: %d\n", failures == 0 ? "ALL PASS" : "FAILURES PRESENT",
           failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
