/* relational.h — the Relational Space algebra of [23 Eq 5.1.1-5.1.5].
 *
 * docs/requirements/icd.md IF-API-003; srs.md §2.4; ADR 0020. A relational
 * element is x = x_i i (+) x_j j (+) x_k k over the spinors of [23 Eq 5.1.1],
 * i being the real unit. The product is the table as printed,
 *
 *     i o i = +1   i o j = j    i o k = k
 *     j o i = j    j o j = -1   j o k = k
 *     k o i = k    k o j = k    k o k = -1,
 *
 * read bilinearly: commutative and NOT associative ((j o j) o k = -k, but
 * j o (j o k) = k). It is adopted literally because it reproduces the EQS
 * brackets of [23 Eq 7.1-7.3] and the associative alternative does not
 * (probes/eqs_relational_product.py). Products of three or more factors are
 * evaluated left to right, never reassociated (FR-REL-004).
 *
 * Pure functions, no allocation, no file-scope state.
 */

#ifndef GIA_RELATIONAL_H
#define GIA_RELATIONAL_H

#include "gia_status.h"

typedef struct { double i, j, k; } rel_t;

/* FR-REL-001: the bilinear product of the table, summed in the fixed order
 * i.i, i.j, i.k, j.i, ..., k.k (numerics N8). */
rel_t rel_mul(rel_t x, rel_t y);

/* FR-REL-004: (x o y) o z, never x o (y o z). */
rel_t rel_mul3(rel_t x, rel_t y, rel_t z);

/* FR-REL-002: the generalised De Moivre exponential [23 Eq 5.1.2],
 * Exp{a i (+) b j (+) c k} = e^a [cos rho + (b j + c k) sin rho / rho],
 * rho = sqrt(b^2 + c^2), with the limit e^a at rho = 0 and sin rho / rho by
 * its series below rho = 1e-4. Never a power series in the algebra. e^a past
 * DBL_MAX is GIA_E_RANGE. */
gia_status rel_exp(rel_t x, rel_t *out, const char **why);

/* FR-REL-003: the canonical ordinal root r_{N,l} (epsilon = 0,
 * sqrt(2) psi = 2 pi l / (N - 1)) [23 Eq A2.5-A2.6]:
 * (cos sqrt2 psi, sin(sqrt2 psi)/sqrt2, sin(sqrt2 psi)/sqrt2). N >= 2,
 * 0 <= l <= N - 1. */
gia_status rel_root(int N, int l, rel_t *out, const char **why);

/* The root raised to the m-th power by multiplying its angle (De Moivre), so
 * rel_root_pow(N, l, N - 1) = 1. The angle is reduced modulo 2 pi in long
 * double. m >= 0. */
gia_status rel_root_pow(int N, int l, int m, rel_t *out, const char **why);

/* The m-fold product x o x o ... o x under the table, left to right; m = 0 is
 * the unit i. Not the root's power: for N = 4, l = 1 the cube is
 * (0.540721, 0, -0.665721), not 1 (PLAN X10). Neither stands in for the other. */
rel_t rel_mul_pow(rel_t x, int m);

#endif /* GIA_RELATIONAL_H */
