/* bench_mop.c — NFR-PERF-001 (T-PERF-01): the First Equation for N = 64
 * components (4,032 couples) at 1,000 output times, affine-power beta, in
 * under 2 s.
 *
 * Verifies: NFR-PERF-001 (T-PERF-01)
 *
 * Each couple gets its own (a, b, p), so no evaluation repeats another; k = 2
 * takes N1's general path. CPU time, not wall time, so a loaded runner does not
 * fail the bound for another job's work. Exit 0 within the bound, 1 otherwise.
 */

#include "mop.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N      64
#define TIMES  1000
#define T_END  2.0
#define BOUND  2.0

int main(void) {
    gia_beta      *beta = (gia_beta *)calloc((size_t)N * N, sizeof(gia_beta));
    gia_matrioska  mat;
    gia_rational   k = {2, 1};
    const char    *why = NULL;
    double complex sum = 0.0;
    clock_t        c0;
    double         secs;
    int            i, j, s;

    if (!beta) return 1;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
            gia_beta *b = &beta[i * N + j];
            if (i == j) { b->kind = GIA_BETA_NONE; continue; }
            b->kind = GIA_BETA_AFFINE;
            b->a    = 1.0 + 0.01 * i + 0.001 * j * I;
            b->b    = 0.1 + 0.001 * j;
            b->p    = 0.5 + 0.01 * ((i + j) % 7);
        }
    c0 = clock();
    for (s = 1; s <= TIMES; s++) {
        if (gia_mop_solve(N, beta, k, T_END * s / TIMES, &mat, &why) != GIA_OK) {
            fprintf(stderr, "bench-mop: solve failed: %s\n", why ? why : "?");
            free(beta);
            return 1;
        }
        sum += mat.a[1];
        gia_matrioska_free(&mat);
    }
    secs = (double)(clock() - c0) / CLOCKS_PER_SEC;
    free(beta);
    printf("bench-mop: N = %d (%d couples), %d times: %.3f s CPU (bound %.1f s)  [%g]\n",
           N, N * (N - 1), TIMES, secs, BOUND, creal(sum));
    return secs < BOUND ? 0 : 1;
}
