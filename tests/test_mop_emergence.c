/* test_mop_emergence.c — the harmony detector (FR-HAR-001, FR-HAR-002, VAL-07).
 *
 * docs/requirements/vv-plan.md §6 (the verdict procedure) and §7 (T-HAR-01..07).
 * This binary is linked WITHOUT harmony.o (FR-HAR-003, T-HAR-01): if the
 * detector or a solver called the constructor gia_harmony_assume_*, it would
 * not link. The oracles are [23 Eq 5.6.5]'s residual evaluated by hand on
 * Matrioskas built here, and the three verdicts PLAN R7 derives from the
 * sources: the First Equation transports harmony, the Second Equation and the
 * EQS impose it.
 */

#define _POSIX_C_SOURCE 200809L

#include "gia_status.h"
#include "mop.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;

static void ok(const char *what, int cond) {
    printf("  %-64s %s\n", what, cond ? "PASS" : "FAIL");
    if (!cond) failures++;
}

/* A Matrioska of size N whose row 1 (index 0) holds row[1..N-1]. */
static gia_matrioska row_matrioska(int N, const double complex *row) {
    gia_matrioska m;
    int           j;
    m.N       = N;
    m.a       = (double complex *)calloc((size_t)N * (size_t)N, sizeof(double complex));
    m.related = (unsigned char *)calloc((size_t)N * (size_t)N, 1);
    if (m.a && m.related)
        for (j = 1; j < N; j++) { m.a[j] = row[j]; m.related[j] = 1; }
    return m;
}

/* vv-plan.md §6's LCG, in [-1, 1]. */
static double lcg(uint64_t *x) {
    *x = 6364136223846793005ULL * *x + 1442695040888963407ULL;
    return 2.0 * ((double)(*x >> 11) * (1.0 / 9007199254740992.0)) - 1.0;
}

/* Verifies: FR-HAR-001 (T-HAR-02, T-HAR-03, T-HAR-04)
 * Source: [23 Eq 5.6.5]; PLAN X9. Oracle: a Matrioska built harmonic by hand
 * (alpha_{1,j+1} = alpha_12 w^{j-1}) has R_H at rounding level; a random
 * complex one does not; a real one cannot get closer than |sin(2 pi/(N-1))|,
 * the imaginary part of the first non-real root it must match. */
static void test_residual(void) {
    static const int Ns[3] = {4, 5, 7};
    const double complex a12 = 0.7 + 0.3 * I;
    double complex row[16];
    const char    *why = NULL;
    double         R;
    int            k, j, all = 1, rnd = 1, real = 1, N;
    uint64_t       x = 7;

    printf("\n[T-HAR-02..04] the harmony residual [23 Eq 5.6.5]\n");
    for (k = 0; k < 3; k++) {
        gia_matrioska m;
        N = Ns[k];
        for (j = 1; j < N; j++) row[j] = a12 * cexp(2.0 * M_PI * I * (double)(j - 1) / (double)(N - 1));
        m = row_matrioska(N, row);
        if (gia_harmony_residual(&m, &R, &why) != GIA_OK || !(R < 1e-9)) all = 0;
        gia_matrioska_free(&m);
    }
    ok("harmonic Matrioska, N = 4, 5, 7: R_H < 1e-9", all);

    for (k = 0; k < 20; k++) {
        gia_matrioska m;
        N = 5;
        for (j = 1; j < N; j++) row[j] = lcg(&x) + I * lcg(&x);
        m = row_matrioska(N, row);
        if (gia_harmony_residual(&m, &R, &why) != GIA_OK || !(R > 1e-2)) rnd = 0;
        gia_matrioska_free(&m);
    }
    ok("20 random complex Matrioskas (LCG), N = 5: R_H > 1e-2", rnd);

    for (N = 4; N <= 8; N++)
        for (k = 0; k < 10; k++) {
            gia_matrioska m;
            for (j = 1; j < N; j++) row[j] = lcg(&x);
            if (row[1] == 0.0) row[1] = 0.5;
            m = row_matrioska(N, row);
            if (gia_harmony_residual(&m, &R, &why) != GIA_OK ||
                !(R >= fabs(sin(2.0 * M_PI / (double)(N - 1))) - 1e-12))
                real = 0;
            gia_matrioska_free(&m);
        }
    ok("real Matrioskas, N = 4..8: R_H >= |sin(2 pi/(N-1))|", real);

    {
        gia_matrioska m;
        int           refused;
        N = 4;
        row[1] = 0.0; row[2] = 1.0; row[3] = 1.0;
        m = row_matrioska(N, row);
        why = NULL;
        refused = gia_harmony_residual(&m, &R, &why) == GIA_E_DOMAIN && why;
        m.a[1] = 1e-320;            /* subnormal: alpha_12 ~ 0 */
        refused = refused && gia_harmony_residual(&m, &R, &why) == GIA_E_DOMAIN;
        ok("alpha_12 = 0 and alpha_12 ~ 1e-320: GIA_E_DOMAIN, with a reason", refused);
        m.a[1] = 1.0; m.related[2] = 0;
        ok("an unrelated row-1 couple: GIA_E_DOMAIN",
           gia_harmony_residual(&m, &R, &why) == GIA_E_DOMAIN);
        m.N = 2;
        ok("N < 3, or NULL: GIA_E_ARG",
           gia_harmony_residual(&m, &R, &why) == GIA_E_ARG &&
           gia_harmony_residual(NULL, &R, &why) == GIA_E_ARG);
        m.N = 4;
        gia_matrioska_free(&m);
    }
}

/* Two constructions with known answers, to show the classifier separates
 * the three verdicts at all: one ignores its input and returns a harmonic
 * row (imposed); one adds a constant to the input's values (absent). */
static gia_status construct_fixed(const gia_beta *in, int N, gia_matrioska *out, void *ctx) {
    double complex row[16];
    int            j;
    (void)in; (void)ctx;
    for (j = 1; j < N; j++) row[j] = cexp(2.0 * M_PI * I * (double)(j - 1) / (double)(N - 1));
    *out = row_matrioska(N, row);
    return GIA_OK;
}

static gia_status construct_shift(const gia_beta *in, int N, gia_matrioska *out, void *ctx) {
    double complex row[16];
    int            j;
    (void)ctx;
    for (j = 1; j < N; j++) row[j] = in[j].a + in[j].b + 0.01;
    *out = row_matrioska(N, row);
    return GIA_OK;
}

static gia_status construct_fails(const gia_beta *in, int N, gia_matrioska *out, void *ctx) {
    (void)in; (void)N; (void)out; (void)ctx;
    return GIA_E_DOMAIN;
}

/* Verifies: FR-HAR-002, BR-007 (T-HAR-05, T-HAR-06, T-HAR-07, VAL-07)
 * Source: PLAN R7; vv-plan.md §6. Oracle: the verdicts R7 derives -- the
 * First Equation carries its input's harmony through (alpha is the integral
 * of beta, so ratios of beta are ratios of alpha: transported); the Second
 * Equation's row and the EQS roots are roots of unity whatever the input
 * (imposed). The catalogue mutation for T-HAR-05, a First Equation
 * construction that returns a constructed (harmonic) Matrioska, reads
 * imposed -- construct_fixed shows the classifier says so. */
static void test_verdicts(void) {
    gia_verdict  v = GIA_H_PRESENT;
    gia_rational k2 = {2, 1};
    const char  *why = NULL;
    int          N, t5 = 1, t6 = 1, t7 = 1;

    printf("\n[T-HAR-05..07, VAL-07] the verdict of each construction\n");
    for (N = 3; N <= 7; N++) {
        if (gia_harmony_verdict(gia_construct_first, NULL, N, &v, &why) != GIA_OK ||
            v != GIA_H_TRANSPORTED) t5 = 0;
        if (gia_harmony_verdict(gia_construct_first, &k2, N, &v, &why) != GIA_OK ||
            v != GIA_H_TRANSPORTED) t5 = 0;
        if (gia_harmony_verdict(gia_construct_second, NULL, N, &v, &why) != GIA_OK ||
            v != GIA_H_IMPOSED) t6 = 0;
        if (gia_harmony_verdict(gia_construct_eqs, NULL, N, &v, &why) != GIA_OK ||
            v != GIA_H_IMPOSED) t7 = 0;
    }
    ok("First Equation (k = 1, 2), N = 3..7: transported", t5);
    ok("Second Equation, N = 3..7: imposed", t6);
    ok("EQS, N = 3..7: imposed", t7);
    ok("a construction that ignores its input: imposed",
       gia_harmony_verdict(construct_fixed, NULL, 5, &v, &why) == GIA_OK && v == GIA_H_IMPOSED);
    ok("a construction that shifts its input: absent",
       gia_harmony_verdict(construct_shift, NULL, 5, &v, &why) == GIA_OK && v == GIA_H_ABSENT);
    ok("a construction's failure is returned",
       gia_harmony_verdict(construct_fails, NULL, 5, &v, &why) == GIA_E_DOMAIN);
    ok("N < 3, NULL construction: GIA_E_ARG",
       gia_harmony_verdict(gia_construct_first, NULL, 2, &v, &why) == GIA_E_ARG &&
       gia_harmony_verdict(NULL, NULL, 5, &v, &why) == GIA_E_ARG);
}

/* Verifies: FR-HAR-002, FR-OUT-003 (T-HAR-02)
 * Oracle: present exactly when R_H <= 1e-9; the verdict strings are
 * FR-OUT-003's four words. */
static void test_observed(void) {
    double complex row[8];
    gia_matrioska  m;
    gia_verdict    v = GIA_H_IMPOSED;
    const char    *why = NULL;
    int            j;

    printf("\n[FR-HAR-002] an observed Matrioska\n");
    for (j = 1; j < 5; j++) row[j] = 2.0 * cexp(2.0 * M_PI * I * (double)(j - 1) / 4.0);
    m = row_matrioska(5, row);
    ok("harmonic: present", gia_harmony_observed(&m, &v, &why) == GIA_OK && v == GIA_H_PRESENT);
    m.a[3] += 1e-6;
    ok("off by 1e-6: absent", gia_harmony_observed(&m, &v, &why) == GIA_OK && v == GIA_H_ABSENT);
    m.a[1] = 0.0;
    ok("alpha_12 = 0: the residual's refusal is returned",
       gia_harmony_observed(&m, &v, &why) == GIA_E_DOMAIN);
    gia_matrioska_free(&m);
    ok("verdict words: imposed, transported, present, absent",
       !strcmp(gia_verdict_str(GIA_H_IMPOSED), "imposed") &&
       !strcmp(gia_verdict_str(GIA_H_TRANSPORTED), "transported") &&
       !strcmp(gia_verdict_str(GIA_H_PRESENT), "present") &&
       !strcmp(gia_verdict_str(GIA_H_ABSENT), "absent") &&
       gia_verdict_str((gia_verdict)99) != NULL);
}

/* Verifies: FR-HAR-002
 * The constructions' own refusals: they need a related beta_12 and N >= 3. */
static void test_constructions(void) {
    gia_beta      in[9];
    gia_matrioska out;
    int           j;

    printf("\n[VAL-07] the constructions refuse what they cannot build\n");
    memset(in, 0, sizeof in);
    for (j = 1; j < 3; j++) { in[j].kind = GIA_BETA_AFFINE; in[j].a = 1.0; in[j].b = 0.25; in[j].p = 1.0; }
    in[1].kind = GIA_BETA_NONE;
    ok("beta_12 unrelated: GIA_E_ARG (second, eqs)",
       gia_construct_second(in, 3, &out, NULL) == GIA_E_ARG &&
       gia_construct_eqs(in, 3, &out, NULL) == GIA_E_ARG);
    in[1].kind = GIA_BETA_AFFINE;
    ok("N < 3: GIA_E_ARG (first, second, eqs)",
       gia_construct_first(in, 2, &out, NULL) == GIA_E_ARG &&
       gia_construct_second(in, 2, &out, NULL) == GIA_E_ARG &&
       gia_construct_eqs(in, 2, &out, NULL) == GIA_E_ARG);
}

int main(void) {
    printf("=== the harmony detector (linked without harmony.o) ===\n");
    test_residual();
    test_verdicts();
    test_observed();
    test_constructions();
    printf("\n%s\nfailures: %d\n", failures == 0 ? "ALL PASS" : "FAILURES PRESENT", failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
