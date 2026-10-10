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
