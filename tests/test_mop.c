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

#include "gia_status.h"
#include "idc.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int main(void) {
    printf("=== Giannantoni kernel: verification and validation ===\n");
    test_status_contract();
    test_idc_general_f();
    test_lde2();
    test_lde2_zero_root();
    test_lde2_linear_in_ics();

    printf("\n%s\nfailures: %d\n", failures == 0 ? "ALL PASS" : "FAILURES PRESENT",
           failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
