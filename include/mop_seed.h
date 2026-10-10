/* mop_seed.h — the seed's `mop` block and the MOP CSV.
 *
 * docs/requirements/icd.md IF-JSON-001, IF-OUT-002; srs.md NFR-ROB-001. The
 * block is parsed strictly: an unknown key anywhere in it, a wrong type, a
 * couple naming a node that is not a component, a couple given twice, or
 * samples whose t does not start at 0 and increase is a load error. A load
 * error is GIA_E_ARG, and `detail` (when non-NULL, `cap` > 0) receives the
 * JSON path of the offending value, so the CLI can name it.
 *
 * This unit sits between the engine's model (engine.h, for component ids) and
 * the First Equation (mop.h). It does not include gssk.h (NFR-SEP-001).
 */

#ifndef GIA_MOP_SEED_H
#define GIA_MOP_SEED_H

#include <complex.h>
#include <stddef.h>

#include "engine.h"
#include "gia_status.h"
#include "mop.h"

/* `"beta": [...]` lists the couples; `"beta": "network"` asks for beta from
 * the network (FR-MOP-008). */
typedef enum { GIA_MOP_BETA_COUPLES, GIA_MOP_BETA_NETWORK } gia_mop_beta_form;

/* One related couple. For SAMPLES, beta.t and beta.v point at t and v, which
 * the seed owns. from and to are node indices into the model. */
typedef struct {
    int             from, to;
    gia_beta        beta;
    double         *t;
    double complex *v;
} gia_mop_couple_spec;

typedef struct {
    int                  present;       /* the seed has a `mop` block */
    gia_rational         k;             /* reduced, den >= 1 */
    int                  ref[2];        /* the reference couple, node indices */
    gia_mop_beta_form    form;
    int                  n_couples;
    gia_mop_couple_spec *couples;       /* sorted by (from id, to id) */
    int                  has_second;    /* FR-MOP-005 */
    double complex       alpha12_0;
    double               c1, c2;
    int                  has_eqs;       /* FR-MOP-006; eqs.N = component count */
    gia_eqs_params       eqs;
} gia_mop_seed;

/* IF-JSON-001. Reads m->root's `mop` object into *out. A seed without one is
 * GIA_OK with out->present = 0. On any failure *out is left zeroed and nothing
 * stays allocated. Free with gia_mop_seed_free. */
gia_status gia_mop_seed_load(const gia_model *m, gia_mop_seed *out, char *detail, size_t cap,
                             const char **why);

void gia_mop_seed_free(gia_mop_seed *s);

/* FR-MOP-008 (PLAN R8) — boundary conditions from the network, at t. For each
 * couple of components (i, j) joined by a direct pathway i -> j that carries
 * quantity, e^{alpha_ij(t)} is the emergy those pathways carry
 * (gia_emergy_carried, summed), so alpha_ij = ln E_ij; and beta_ij follows from
 * [23 Eq 5.5.2] with k = 1, (d~/dt) alpha = alpha', so beta_ij = E_ij'/E_ij,
 * by a central difference (forward at t < h). Every other couple is
 * unrelated. Both Matrioskas are N = n_nodes, row-major by node index, real;
 * either may be NULL. A related couple whose E is 0 near t has no logarithm:
 * GIA_E_RANGE. */
gia_status gia_mop_network(const gia_model *m, double t, gia_matrioska *alpha,
                           gia_matrioska *beta, const char **why);

/* FR-HAR-001 over a model: the reference couple's row, as [23 Eq 5.6.5]
 * reads it. From `full` (N = n_nodes, by node index) this builds a Matrioska
 * over the K components: index 0 is ref[0], index 1 ref[1], then the other
 * components in id order; only row 0 is filled, from ref[0]'s couples, and a
 * couple `full` leaves unrelated stays unrelated. ref NULL means the first two
 * components by id. Allocates *row (free with gia_matrioska_free). GIA_E_ARG
 * for fewer than two components or a ref that is not two distinct ones. */
gia_status gia_mop_reference_row(const gia_model *m, const int *ref, const gia_matrioska *full,
                                 gia_matrioska *row, const char **why);

/* IF-OUT-002. Writes `time`, then `<from>__<to>_re, <from>__<to>_im` per
 * couple in (from, to) id order, on the trajectory CSV's grid
 * t = s t_end / steps, s = 0..steps. Every couple is solved at t_end before
 * the file is opened, so a refusal (the First Equation's GIA_E_DOMAIN, or
 * GIA_E_UNSUPPORTED for the network form with k != 1, PLAN R8) writes nothing;
 * a later failure removes the file. The network form writes alpha = ln E per
 * related couple (gia_mop_network). The last column, R_H, is the harmony residual
 * (FR-HAR-001) of the reference row at that time (gia_mop_reference_row), and
 * is left empty where the residual is undefined: N < 3, a row couple
 * unrelated, or alpha_12 = 0. */
gia_status gia_mop_write_csv(const gia_model *m, const gia_mop_seed *s, const char *path,
                             int steps, const char **why);

#endif /* GIA_MOP_SEED_H */
