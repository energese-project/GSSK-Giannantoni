/* harmony.c — the harmony constructor: the Harmony Relationships ASSUMED.
 *
 * PLAN E5, W8; srs FR-HAR-004. This builds the N x N matrix of ordinal
 * relationships from one reference couple and the roots of unity. It does not
 * derive harmony from any equation: [23 §8 ii] says the EQS assumes it, and
 * everything this file produces is labelled `assumed` (FR-OUT-001). The
 * detector, which does derive a verdict, is gia_harmony_residual and
 * gia_harmony_verdict in mop.c; no solver and no detector calls anything here
 * (FR-HAR-003, T-HAR-01).
 */

#include "engine.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== *
 * The constructor (FR-HAR-004)
 * ================================================================== */

double complex gia_harmony_assume_root(int roots, int m) {
    double theta;
    if (roots < 1) return 1.0;
    theta = 2.0 * M_PI * (double)m / (double)roots;
    return cos(theta) + I * sin(theta);
}

double complex gia_harmony_assume_reconstruct(int n, double complex alpha_ref,
                                       int i, int j) {
    int m;
    if (n < 2 || i < 0 || j < 0 || i >= n || j >= n) return 0.0;
    if (i == j) return 0.0;  /* no self-relation on the diagonal */
    /* Within row i the partners j != i are indexed in ascending order, so the
     * partner's ordinal position is j, less one if it sits past the diagonal. */
    m = (j < i) ? j : j - 1;
    return alpha_ref * gia_harmony_assume_root(n - 1, m);
}

bool gia_harmony_assume_init(gia_harmony *h, int n, double complex alpha_ref) {
    int i, j;
    if (!h || n < 2) return false;
    h->n         = n;
    h->alpha_ref = alpha_ref;
    h->a         = (double complex *)calloc((size_t)n * (size_t)n,
                                            sizeof(double complex));
    if (!h->a) { h->n = 0; return false; }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            h->a[(size_t)i * (size_t)n + (size_t)j] =
                gia_harmony_assume_reconstruct(n, alpha_ref, i, j);
    return true;
}

void gia_harmony_assume_free(gia_harmony *h) {
    if (!h) return;
    free(h->a);
    h->a = NULL;
    h->n = 0;
}

double complex gia_harmony_assume_at(const gia_harmony *h, int i, int j) {
    if (!h || !h->a || i < 0 || j < 0 || i >= h->n || j >= h->n) return 0.0;
    return h->a[(size_t)i * (size_t)h->n + (size_t)j];
}

double gia_harmony_assume_row_residual(const gia_harmony *h) {
    double worst = 0.0;
    int    i, j;
    if (!h || !h->a) return 0.0;
    for (i = 0; i < h->n; i++) {
        double complex s = 0.0;
        double         mod;
        for (j = 0; j < h->n; j++) s += gia_harmony_assume_at(h, i, j);
        mod = cabs(s);
        if (mod > worst) worst = mod;
    }
    return worst;
}

double gia_harmony_assume_reduction_residual(const gia_harmony *h) {
    double worst = 0.0;
    int    i, j;
    if (!h || !h->a) return 0.0;
    for (i = 0; i < h->n; i++) {
        for (j = 0; j < h->n; j++) {
            double d = cabs(gia_harmony_assume_at(h, i, j) -
                            gia_harmony_assume_reconstruct(h->n, h->alpha_ref, i, j));
            if (d > worst) worst = d;
        }
    }
    return worst;
}

