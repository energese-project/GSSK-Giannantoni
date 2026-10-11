/* gia_status.h — status codes of the Giannantoni library units.
 *
 * docs/requirements/icd.md IF-API-001. Every function of idc.h, relational.h
 * and mop.h returns one of these and takes `const char **why` last. On any
 * status but GIA_OK, a non-NULL `why` is set to a static string naming the
 * cause and its source (an equation, or a PLAN.md §5/§6 row), and the
 * function's outputs are left unmodified. On GIA_OK, `*why` is not written.
 * No function prints, exits or aborts (NFR-ERR-001).
 *
 * This header includes nothing but <stddef.h>, and no Giannantoni header
 * includes gssk.h (NFR-SEP-001).
 */

#ifndef GIA_STATUS_H
#define GIA_STATUS_H

#include <stddef.h>

typedef enum {
    GIA_OK = 0,
    GIA_E_ARG,          /* NULL where not allowed, size or index out of range        */
    GIA_E_DOMAIN,       /* mathematically outside the function's stated domain       */
    GIA_E_RANGE,        /* result would overflow or be non-finite (NFR-NUM-003)      */
    GIA_E_CONVERGENCE,  /* quadrature or root tracking failed its tolerance          */
    GIA_E_UNSUPPORTED,  /* defined by no available source (srs FR-IDC-013, PLAN §6)  */
    GIA_E_LIMIT,        /* a fixed limit exceeded (NFR-LIM-001)                      */
    GIA_E_NOMEM
} gia_status;

/* Never NULL: a status outside the enum reads "unknown status". */
const char *gia_status_str(gia_status s);

#endif /* GIA_STATUS_H */
