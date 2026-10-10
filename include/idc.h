/* idc.h — single-variable Incipient Differential Calculus.
 *
 * docs/requirements/icd.md IF-API-002; srs.md §2.1. Each function implements
 * one requirement, named beside it, from the source equation it cites. The
 * derivative of a single function given pointwise is the incipient one,
 *
 *     (d~/dt)^n f = (f'/f)^n f            [02 Eq 14.9.5], [23 Eq 5.5.2]
 *
 * and of a superposition of exponential terms it is taken termwise (srs §1,
 * PLAN R16).
 *
 * Status of record: this header is the incipient calculus. The kernel's
 * "incipient" method is not (PLAN.md §2 E2, ADR 0022).
 */

#ifndef GIA_IDC_H
#define GIA_IDC_H

#include <complex.h>

#include "gia_status.h"

/* FR-IDC-002. (f'/f)^n * f for f != 0 and integer n >= 0, in C. The power is
 * taken by repeated multiplication, never cpow, so its sign and branch are
 * exact (numerics.md §1). f == 0 is GIA_E_DOMAIN; n < 0 or out == NULL is
 * GIA_E_ARG; a non-finite result is GIA_E_RANGE. */
gia_status gia_idc_of(double complex f, double complex df, int n,
                      double complex *out, const char **why);

#endif /* GIA_IDC_H */
