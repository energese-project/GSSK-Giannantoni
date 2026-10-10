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

int main(void) {
    printf("=== Giannantoni kernel: verification and validation ===\n");
    test_status_contract();
    test_idc_general_f();

    printf("\n%s\nfailures: %d\n", failures == 0 ? "ALL PASS" : "FAILURES PRESENT",
           failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
