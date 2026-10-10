/* relational.c — the Relational Space algebra (include/relational.h).
 *
 * docs/requirements/srs.md §2.4, numerics.md N8, ADR 0020. Pure functions:
 * no allocation, no file-scope mutable state (NFR-REE-001), nothing clamped.
 */

#include "relational.h"

#include <math.h>

#define EXP_MAX 709.78                  /* log(DBL_MAX), numerics N6 */

static gia_status fail(gia_status s, const char *reason, const char **why) {
    if (why) *why = reason;
    return s;
}

/* The nine terms of [23 Eq 5.1.3-5.1.5], in the order i.i, i.j, i.k, j.i,
 * j.j, j.k, k.i, k.j, k.k (N8), so the sum is the same on every toolchain. */
rel_t rel_mul(rel_t x, rel_t y) {
    rel_t r;
    r.i = 0.0; r.j = 0.0; r.k = 0.0;
    r.i += x.i * y.i;          /* i o i = +1 */
    r.j += x.i * y.j;          /* i o j = j  */
    r.k += x.i * y.k;          /* i o k = k  */
    r.j += x.j * y.i;          /* j o i = j  */
    r.i -= x.j * y.j;          /* j o j = -1 */
    r.k += x.j * y.k;          /* j o k = k  */
    r.k += x.k * y.i;          /* k o i = k  */
    r.k += x.k * y.j;          /* k o j = k  */
    r.i -= x.k * y.k;          /* k o k = -1 */
    return r;
}

/* FR-REL-004: left to right, as written. */
rel_t rel_mul3(rel_t x, rel_t y, rel_t z) {
    return rel_mul(rel_mul(x, y), z);
}

/* FR-REL-002 — [23 Eq 5.1.2]. */
gia_status rel_exp(rel_t x, rel_t *out, const char **why) {
    double rho, sinc, ea;
    if (!out) return fail(GIA_E_ARG, "rel_exp: out is NULL", why);
    if (!isfinite(x.i) || !isfinite(x.j) || !isfinite(x.k))
        return fail(GIA_E_DOMAIN, "rel_exp: the element is not finite", why);
    if (x.i > EXP_MAX)
        return fail(GIA_E_RANGE, "rel_exp: e^a overflows (numerics N6)", why);
    rho  = hypot(x.j, x.k);
    sinc = rho < 1e-4 ? 1.0 - rho * rho / 6.0 + rho * rho * rho * rho / 120.0
                      : sin(rho) / rho;
    ea   = exp(x.i);
    out->i = ea * cos(rho);
    out->j = ea * x.j * sinc;
    out->k = ea * x.k * sinc;
    return GIA_OK;
}

/* The root at angle phi = sqrt2 psi: (cos phi, sin phi/sqrt2, sin phi/sqrt2). */
static rel_t root_at(long double phi) {
    rel_t r;
    double s = (double)sinl(phi) / sqrt(2.0);
    r.i = (double)cosl(phi); r.j = s; r.k = s;
    return r;
}

static long double root_angle(int N, int l) {
    return 2.0L * 3.141592653589793238462643383279502884L * (long double)l / (long double)(N - 1);
}

/* FR-REL-003 — [23 Eq A2.5-A2.6], epsilon = 0. */
gia_status rel_root(int N, int l, rel_t *out, const char **why) {
    if (!out) return fail(GIA_E_ARG, "rel_root: out is NULL", why);
    if (N < 2 || l < 0 || l > N - 1)
        return fail(GIA_E_ARG, "rel_root: needs N >= 2 and 0 <= l <= N - 1", why);
    *out = root_at(root_angle(N, l));
    return GIA_OK;
}

/* The m-th power by angle (De Moivre): unity at m = N - 1 (PLAN R3). */
gia_status rel_root_pow(int N, int l, int m, rel_t *out, const char **why) {
    const long double two_pi = 2.0L * 3.141592653589793238462643383279502884L;
    long double phi;
    if (!out) return fail(GIA_E_ARG, "rel_root_pow: out is NULL", why);
    if (N < 2 || l < 0 || l > N - 1 || m < 0)
        return fail(GIA_E_ARG, "rel_root_pow: needs N >= 2, 0 <= l <= N - 1, m >= 0", why);
    /* l m / (N - 1) turns, reduced exactly in integers, then to an angle. */
    phi = two_pi * (long double)(((long long)l * (long long)m) % (long long)(N - 1)) /
          (long double)(N - 1);
    *out = root_at(phi);
    return GIA_OK;
}

/* The table power: x o x o ... o x, left to right (PLAN X10). */
rel_t rel_mul_pow(rel_t x, int m) {
    rel_t r;
    int   n;
    r.i = 1.0; r.j = 0.0; r.k = 0.0;
    for (n = 0; n < m; n++) r = rel_mul(r, x);
    return r;
}
