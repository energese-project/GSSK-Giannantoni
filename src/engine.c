/* engine.c — Giannantoni Generative Computational Framework.
 *
 * See include/engine.h for the theory this implements and for why this file
 * does not include "gssk.h".
 *
 * There is no main() here. The CLI lives in src/sim_main.c so that the
 * engine can be linked into tests without dragging an entry point along.
 */

#include "engine.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Storage levels are put into exponential form via ln, so they need a floor. */
#define GIA_EPS 1e-12

/* ================================================================== *
 * 1. Exponential form
 * ================================================================== */

double gia_phi_eval(const gia_phi *p, double t) {
    double acc = 0.0, tp = 1.0;
    int k;
    if (!p) return 0.0;
    for (k = 0; k <= p->degree; k++) {
        acc += p->c[k] * tp;
        tp  *= t;
    }
    return acc;
}

double gia_phi_deriv(const gia_phi *p, int order, double t) {
    double acc = 0.0;
    int k;
    if (!p || order < 0) return 0.0;
    if (order == 0) return gia_phi_eval(p, t);

    /* d^n/dt^n (c_k t^k) = c_k * k!/(k-n)! * t^(k-n), zero for k < n. */
    for (k = order; k <= p->degree; k++) {
        double falling = 1.0;
        int    j;
        for (j = 0; j < order; j++) falling *= (double)(k - j);
        acc += p->c[k] * falling * pow(t, (double)(k - order));
    }
    return acc;
}

/* ================================================================== *
 * 2. Incipient vs. traditional, and the drift
 * ================================================================== */

double gia_idc_amplitude(const gia_phi *p, int n, double t) {
    double phi_p;
    if (n < 0) return 0.0;
    if (n == 0) return 1.0;
    /* Persistence of form: the amplitude is a bare power of phi', with no
     * contribution from any higher derivative. That absence is the point. */
    phi_p = gia_phi_deriv(p, 1, t);
    return pow(phi_p, (double)n);
}

double gia_tdc_amplitude(const gia_phi *p, int n, double t) {
    /* Complete Bell polynomial B_n(x_1, ..., x_n) with x_j = phi^(j)(t),
     * via the standard recursion
     *
     *     B_0 = 1,   B_{m+1} = sum_{k=0}^{m} C(m,k) B_{m-k} x_{k+1}
     *
     * which is exactly Faa di Bruno for the outer function exp. */
    double x[GIA_MAX_ORDER + 1];
    double b[GIA_MAX_ORDER + 1];
    double binom[GIA_MAX_ORDER + 1][GIA_MAX_ORDER + 1];
    int    i, j, m, k;

    if (n < 0) return 0.0;
    if (n == 0) return 1.0;
    if (n > GIA_MAX_ORDER) n = GIA_MAX_ORDER;

    for (i = 0; i <= GIA_MAX_ORDER; i++) {
        for (j = 0; j <= GIA_MAX_ORDER; j++) binom[i][j] = 0.0;
        binom[i][0] = 1.0;
        for (j = 1; j <= i; j++)
            binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j];
    }

    for (j = 1; j <= n; j++) x[j] = gia_phi_deriv(p, j, t);

    b[0] = 1.0;
    for (m = 0; m < n; m++) {
        double acc = 0.0;
        for (k = 0; k <= m; k++) acc += binom[m][k] * b[m - k] * x[k + 1];
        b[m + 1] = acc;
    }
    return b[n];
}

double gia_idc_derivative(const gia_phi *p, int n, double t) {
    return gia_idc_amplitude(p, n, t) * exp(gia_phi_eval(p, t));
}

double gia_tdc_derivative(const gia_phi *p, int n, double t) {
    return gia_tdc_amplitude(p, n, t) * exp(gia_phi_eval(p, t));
}

double gia_drift(const gia_phi *p, int n, double t) {
    return gia_tdc_amplitude(p, n, t) - gia_idc_amplitude(p, n, t);
}

bool gia_drift_free(const gia_phi *p) {
    int k;
    if (!p) return true;
    /* An affine phi kills every x_j for j >= 2, collapsing the Bell recursion
     * to B_n = (phi')^n. The two calculi then agree identically, at all
     * orders and all t -- no evaluation required to know it. */
    for (k = 2; k <= p->degree; k++)
        if (p->c[k] != 0.0) return false;
    return true;
}

/* ================================================================== *
 * 3. Binary / duet / n-et
 * ================================================================== */

gia_net gia_incipient_fractional(const gia_phi *p, int num, int den, double t) {
    gia_net out;
    double  phi_p, mag, arg, r, q, e;
    int     k;

    memset(&out, 0, sizeof(out));
    if (den < 1)                den = 1;
    if (den > GIA_MAX_BRANCHES) den = GIA_MAX_BRANCHES;
    out.count = den;

    phi_p = gia_phi_deriv(p, 1, t);
    e     = exp(gia_phi_eval(p, t));
    q     = (double)num / (double)den;

    /* (phi^o)^(num/den) taken over every branch of the den-th root. */
    mag = fabs(phi_p);
    arg = (phi_p < 0.0) ? M_PI : 0.0;
    r   = pow(mag, q);

    for (k = 0; k < den; k++) {
        double theta = (arg + 2.0 * M_PI * (double)k) * q;
        out.branch[k] = r * (cos(theta) + I * sin(theta)) * e;
    }
    return out;
}

bool gia_net_is_binary(const gia_net *n) {
    return n && n->count == 2;
}

double complex gia_net_sum(const gia_net *n) {
    double complex acc = 0.0;
    int k;
    if (!n) return 0.0;
    for (k = 0; k < n->count; k++) acc += n->branch[k];
    return acc;
}

/* ================================================================== *
 * 4. MOP Harmony Relationships
 * ================================================================== */

double complex gia_ordinal_root(int roots, int m) {
    double theta;
    if (roots < 1) return 1.0;
    theta = 2.0 * M_PI * (double)m / (double)roots;
    return cos(theta) + I * sin(theta);
}

double complex gia_harmony_reconstruct(int n, double complex alpha_ref,
                                       int i, int j) {
    int m;
    if (n < 2 || i < 0 || j < 0 || i >= n || j >= n) return 0.0;
    if (i == j) return 0.0;  /* no self-relation on the diagonal */
    /* Within row i the partners j != i are indexed in ascending order, so the
     * partner's ordinal position is j, less one if it sits past the diagonal. */
    m = (j < i) ? j : j - 1;
    return alpha_ref * gia_ordinal_root(n - 1, m);
}

bool gia_harmony_init(gia_harmony *h, int n, double complex alpha_ref) {
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
                gia_harmony_reconstruct(n, alpha_ref, i, j);
    return true;
}

void gia_harmony_free(gia_harmony *h) {
    if (!h) return;
    free(h->a);
    h->a = NULL;
    h->n = 0;
}

double complex gia_harmony_at(const gia_harmony *h, int i, int j) {
    if (!h || !h->a || i < 0 || j < 0 || i >= h->n || j >= h->n) return 0.0;
    return h->a[(size_t)i * (size_t)h->n + (size_t)j];
}

double gia_harmony_row_residual(const gia_harmony *h) {
    double worst = 0.0;
    int    i, j;
    if (!h || !h->a) return 0.0;
    for (i = 0; i < h->n; i++) {
        double complex s = 0.0;
        double         mod;
        for (j = 0; j < h->n; j++) s += gia_harmony_at(h, i, j);
        mod = cabs(s);
        if (mod > worst) worst = mod;
    }
    return worst;
}

double gia_harmony_reduction_residual(const gia_harmony *h) {
    double worst = 0.0;
    int    i, j;
    if (!h || !h->a) return 0.0;
    for (i = 0; i < h->n; i++) {
        for (j = 0; j < h->n; j++) {
            double d = cabs(gia_harmony_at(h, i, j) -
                            gia_harmony_reconstruct(h->n, h->alpha_ref, i, j));
            if (d > worst) worst = d;
        }
    }
    return worst;
}

/* ================================================================== *
 * 4b. The network: flow matrix and matrix exponential
 * ================================================================== */

const char *gia_role_name(gia_role r) {
    switch (r) {
        case GIA_ROLE_ENERGY:      return "energy";
        case GIA_ROLE_CONTROL:     return "control";
        case GIA_ROLE_USED:        return "used";
        case GIA_ROLE_GOODS_IN:    return "goods_in";
        case GIA_ROLE_GOODS_OUT:   return "goods_out";
        case GIA_ROLE_COUNTER_IN:  return "counter_in";
        case GIA_ROLE_COUNTER_OUT: return "counter_out";
        default:                   return "none";
    }
}

/* True for a role that names an OUTGOING leg. Every other role names an input.
 * A work gate's product legs carry no role; its one named output is the used
 * energy, which goes to a sink rather than on to the product. */
static bool role_is_outgoing(gia_role r) {
    return r == GIA_ROLE_GOODS_OUT || r == GIA_ROLE_COUNTER_OUT ||
           r == GIA_ROLE_USED;
}

bool gia_node_is_module(const gia_model *m, int node_idx) {
    if (!m || node_idx < 0 || node_idx >= m->n_nodes) return false;
    return m->nodes[node_idx].is_module;
}

const char *gia_action_name(gia_action a) {
    switch (a) {
        case GIA_ACTION_MULTIPLY: return "multiply";
        case GIA_ACTION_DIVIDE:   return "divide";
        case GIA_ACTION_SUBTRACT: return "subtract";
    }
    return "?";
}

const char *gia_logic_name(gia_logic l) {
    switch (l) {
        case GIA_LOGIC_LINEAR:      return "linear";
        case GIA_LOGIC_REVERSIBLE:  return "reversible";
        case GIA_LOGIC_CONSTANT:    return "constant";
    }
    return "unknown";
}

bool gia_matrix_init(gia_matrix *m, int n) {
    if (!m || n <= 0) return false;
    m->n = n;
    m->a = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    if (!m->a) { m->n = 0; return false; }
    return true;
}

void gia_matrix_free(gia_matrix *m) {
    if (!m) return;
    free(m->a);
    m->a = NULL;
    m->n = 0;
}

double gia_matrix_at(const gia_matrix *m, int i, int j) {
    if (!m || !m->a || i < 0 || j < 0 || i >= m->n || j >= m->n) return 0.0;
    return m->a[(size_t)i * (size_t)m->n + (size_t)j];
}

static void mat_mul(const double *x, const double *y, double *out, int n) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double acc = 0.0;
            for (k = 0; k < n; k++)
                acc += x[(size_t)i * (size_t)n + (size_t)k] *
                       y[(size_t)k * (size_t)n + (size_t)j];
            out[(size_t)i * (size_t)n + (size_t)j] = acc;
        }
    }
}

/* Solve D X = N for X by Gauss-Jordan with partial pivoting. */
static bool mat_solve(double *D, double *N, double *X, int n) {
    int i, j, k, piv;
    double *a = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    double *b = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    if (!a || !b) { free(a); free(b); return false; }
    memcpy(a, D, (size_t)n * (size_t)n * sizeof(double));
    memcpy(b, N, (size_t)n * (size_t)n * sizeof(double));

    for (k = 0; k < n; k++) {
        double best = fabs(a[(size_t)k * (size_t)n + (size_t)k]), p;
        piv = k;
        for (i = k + 1; i < n; i++) {
            p = fabs(a[(size_t)i * (size_t)n + (size_t)k]);
            if (p > best) { best = p; piv = i; }
        }
        if (best < 1e-300) { free(a); free(b); return false; }
        if (piv != k) {
            for (j = 0; j < n; j++) {
                double t;
                t = a[(size_t)k*(size_t)n+(size_t)j];
                a[(size_t)k*(size_t)n+(size_t)j] = a[(size_t)piv*(size_t)n+(size_t)j];
                a[(size_t)piv*(size_t)n+(size_t)j] = t;
                t = b[(size_t)k*(size_t)n+(size_t)j];
                b[(size_t)k*(size_t)n+(size_t)j] = b[(size_t)piv*(size_t)n+(size_t)j];
                b[(size_t)piv*(size_t)n+(size_t)j] = t;
            }
        }
        {
            double d = a[(size_t)k*(size_t)n+(size_t)k];
            for (j = 0; j < n; j++) {
                a[(size_t)k*(size_t)n+(size_t)j] /= d;
                b[(size_t)k*(size_t)n+(size_t)j] /= d;
            }
        }
        for (i = 0; i < n; i++) {
            double f;
            if (i == k) continue;
            f = a[(size_t)i*(size_t)n+(size_t)k];
            if (f == 0.0) continue;
            for (j = 0; j < n; j++) {
                a[(size_t)i*(size_t)n+(size_t)j] -= f * a[(size_t)k*(size_t)n+(size_t)j];
                b[(size_t)i*(size_t)n+(size_t)j] -= f * b[(size_t)k*(size_t)n+(size_t)j];
            }
        }
    }
    memcpy(X, b, (size_t)n * (size_t)n * sizeof(double));
    free(a); free(b);
    return true;
}

/* exp(M) by Pade (3,3) with scaling and squaring:
 *     N = 120I + 60X + 12X^2 + X^3,  D = 120I - 60X + 12X^2 - X^3,  exp = N/D
 * The same order the kernel's IDC path uses, so the two agree on the linear
 * core rather than drifting apart on solver choice. */
bool gia_matrix_exp(const gia_matrix *M, gia_matrix *B) {
    int    n, i, j, sq = 0, s;
    size_t sz;
    double norm = 0.0, scale = 1.0;
    double *X = NULL, *X2 = NULL, *X3 = NULL, *Nm = NULL, *Dm = NULL, *tmp = NULL;
    bool ok = false;

    if (!M || !M->a || !B) return false;
    n  = M->n;
    sz = (size_t)n * (size_t)n * sizeof(double);
    if (!gia_matrix_init(B, n)) return false;

    /* Scale so that ||X||_inf is small, then square back up.
     *
     * The threshold matters more than it looks. Pade (3,3) truncates at
     * O(||X||^7), so scaling only to 1/2 leaves an error near 1e-7 -- which is
     * fine for a stepped solver taking a small dt, and useless here, where the
     * whole claim is that the incipient solution is exact and closed-form. At
     * 2^-6 the truncation term falls below double precision, and the cost is a
     * few extra squarings of a small matrix. */
    for (i = 0; i < n; i++) {
        double row = 0.0;
        for (j = 0; j < n; j++) row += fabs(M->a[(size_t)i*(size_t)n+(size_t)j]);
        if (row > norm) norm = row;
    }
    while (norm * scale > 0.015625) { scale *= 0.5; sq++; }

    /* calloc, not malloc: the Pade numerator and denominator are filled by a
     * loop over a runtime n, which GCC cannot prove covers the whole buffer --
     * it reports 'may be used uninitialized' under -Werror where clang does
     * not. Zeroing costs nothing here and the buffers are small. */
    X  = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    X2 = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    X3 = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    Nm = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    Dm = (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    tmp= (double *)calloc((size_t)n * (size_t)n, sizeof(double));
    if (!X || !X2 || !X3 || !Nm || !Dm || !tmp) goto done;

    for (i = 0; i < n * n; i++) X[i] = M->a[i] * scale;
    mat_mul(X, X, X2, n);
    mat_mul(X2, X, X3, n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            size_t k = (size_t)i * (size_t)n + (size_t)j;
            double id = (i == j) ? 120.0 : 0.0;
            Nm[k] = id + 60.0 * X[k] + 12.0 * X2[k] + X3[k];
            Dm[k] = id - 60.0 * X[k] + 12.0 * X2[k] - X3[k];
        }
    }
    if (!mat_solve(Dm, Nm, B->a, n)) goto done;

    for (s = 0; s < sq; s++) {
        mat_mul(B->a, B->a, tmp, n);
        memcpy(B->a, tmp, sz);
    }
    ok = true;

done:
    free(X); free(X2); free(X3); free(Nm); free(Dm); free(tmp);
    if (!ok) gia_matrix_free(B);
    return ok;
}

/* ---- the flow matrix, assembled from Odum's pathway laws ---- */

/* ---- forcing carried as state ----
 *
 * Each forced component contributes extra rows to the augmented system, laid
 * out as [ Q_0..Q_{n-1} | extras | 1 ]. The extras are the waveform's own
 * generator, so A stays constant and the matrix exponential stays exact.
 *
 * A sine needs two (sin and cos, which generate each other); a ramp and an
 * exponential need one each. */
static int forcing_width(gia_forcing_kind k) {
    switch (k) {
        case GIA_FORCE_SINE:        return 2;
        case GIA_FORCE_RAMP:        return 1;
        case GIA_FORCE_EXPONENTIAL: return 1;
        default:                    return 0;
    }
}

/* Index of a node's first extra state, or -1 if it is not forced. */
static int forcing_slot(const gia_model *m, int node) {
    int i, at = m->n_nodes;
    for (i = 0; i < m->n_nodes; i++) {
        int w = forcing_width(m->nodes[i].forcing.kind);
        if (i == node) return w ? at : -1;
        at += w;
    }
    return -1;
}

static void forcing_decompose(const gia_node *nd, double *konst, double *coeff);

/* Add `g` times the origin's value into `row` of the augmented matrix.
 *
 * Normally the origin's value is a state, so this is one entry in its column.
 * A FORCED held component is different: its value is not the state Q_a, it is
 * the waveform, so the term is split across the phantom column (the constant
 * part) and the waveform's own extra state (the varying part). Reading column
 * a for a forced source would silently use its declared value and ignore the
 * driver entirely. */
static void add_origin_term(const gia_model *m, gia_matrix *out, int dim,
                            int row, int a, double g) {
    const gia_node *nd = &m->nodes[a];
    int slot;

    if (nd->integrates || nd->forcing.kind == GIA_FORCE_NONE) {
        out->a[(size_t)row*(size_t)dim+(size_t)a] += g;
        return;
    }
    slot = forcing_slot(m, a);
    {
        double konst, coeff;
        forcing_decompose(nd, &konst, &coeff);
        out->a[(size_t)row*(size_t)dim+(size_t)(dim-1)] += g * konst;
        if (slot >= 0) out->a[(size_t)row*(size_t)dim+(size_t)slot] += g * coeff;
    }
}

static int forcing_extra_count(const gia_model *m) {
    int i, n = 0;
    for (i = 0; i < m->n_nodes; i++) n += forcing_width(m->nodes[i].forcing.kind);
    return n;
}

double gia_forcing_value(const gia_forcing *f, double base, double t) {
    if (!f) return base;
    switch (f->kind) {
        case GIA_FORCE_SINE:
            return f->offset + f->amplitude * sin(f->rate * t + f->phase);
        case GIA_FORCE_RAMP:
            return f->offset + f->rate * t;
        case GIA_FORCE_EXPONENTIAL:
            return f->offset + f->amplitude * exp(f->rate * t);
        default:
            return base;
    }
}

/* The driver's value written as (constant part, coefficient on its first extra
 * state). A forced node's contribution to a pathway is then
 * k*const on the phantom column plus k*coeff on the extra column. */
static void forcing_decompose(const gia_node *nd, double *konst, double *coeff) {
    const gia_forcing *f = &nd->forcing;
    switch (f->kind) {
        case GIA_FORCE_SINE:
            /* sin(w t + p) = sin(p) cos(w t) + cos(p) sin(w t); the two extras
             * are carried as s = sin(w t + p), c = cos(w t + p), so the value
             * is just the s slot. */
            *konst = f->offset; *coeff = f->amplitude; break;
        case GIA_FORCE_RAMP:
            *konst = f->offset; *coeff = f->rate; break;
        case GIA_FORCE_EXPONENTIAL:
            *konst = f->offset; *coeff = f->amplitude; break;
        default:
            *konst = nd->q0;    *coeff = 0.0; break;
    }
}

/* True when a threshold pathway is open at the operating point. The clamp on
 * `subtract` is the same kind of boundary and is treated the same way. */
bool gia_build_flow_matrix(const gia_model *m, const double *q, double t,
                           gia_matrix *out) {
    int i, n, dim;
    if (!m || !out || m->n_nodes <= 0) return false;
    n   = m->n_nodes;
    dim = n + forcing_extra_count(m) + 1;   /* [ Q | extras | 1 ] */
    if (!gia_matrix_init(out, dim)) return false;

    /* Each forced component's waveform generates itself, so its rows are
     * constant and the whole augmented system stays time-invariant. */
    for (i = 0; i < n; i++) {
        const gia_forcing *f = &m->nodes[i].forcing;
        int slot = forcing_slot(m, i);
        if (slot < 0) continue;
        switch (f->kind) {
            case GIA_FORCE_SINE:
                /* d/dt s = w c,  d/dt c = -w s */
                out->a[(size_t)slot*(size_t)dim+(size_t)(slot+1)]     =  f->rate;
                out->a[(size_t)(slot+1)*(size_t)dim+(size_t)slot]     = -f->rate;
                break;
            case GIA_FORCE_RAMP:
                /* d/dt r = 1, taken from the phantom component */
                out->a[(size_t)slot*(size_t)dim+(size_t)(dim-1)]      =  1.0;
                break;
            case GIA_FORCE_EXPONENTIAL:
                /* d/dt e = lambda e */
                out->a[(size_t)slot*(size_t)dim+(size_t)slot]         =  f->rate;
                break;
            default:
                break;
        }
    }

    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *e = &m->edges[i];
        int    a = e->from, b = e->to;
        double g = 0.0, k_t;
        bool   drain_a, fill_b;

        if (a < 0 || b < 0) continue;

        /* A pathway touching a module carries no law of its own: the module's
         * law reads all of them together, so they are handled in the module
         * pass below. Leaving them here would count the flow twice. */
        if (m->nodes[a].is_module || m->nodes[b].is_module) continue;

        /* ADR 0006 edge attachment: the waveform drives the RATE. Evaluated at
         * t rather than carried as state, because k(t)*Q is bilinear and there
         * is no augmentation that makes it linear again. */
        k_t = gia_forcing_value(&e->forcing, e->weight, t);

        /* Odum SecV: a heat sink absorbs and is never depleted, so a pathway
         * leaving one contributes no drain term. A held component (source or
         * constant, SecII) is likewise never drained. */
        drain_a = m->nodes[a].integrates && m->nodes[a].kind != GIA_NODE_SINK;
        fill_b  = m->nodes[b].integrates;

        switch (e->logic) {
            case GIA_LOGIC_CONSTANT:
                /* F = k, independent of state: an affine term, so it lands in
                 * the augmented column rather than in a conductance. */
                if (drain_a) out->a[(size_t)a*(size_t)dim+(size_t)n] -= k_t;
                if (fill_b)  out->a[(size_t)b*(size_t)dim+(size_t)n] += k_t;
                continue;

            case GIA_LOGIC_LINEAR:
                g = k_t;
                break;

            case GIA_LOGIC_REVERSIBLE:
                g = k_t;
                break;
        }

        if (drain_a) add_origin_term(m, out, dim, a, a, -g);
        if (fill_b)  add_origin_term(m, out, dim, b, a,  g);

        if (e->logic == GIA_LOGIC_REVERSIBLE) {
            if (drain_a) add_origin_term(m, out, dim, a, b,  g);
            if (fill_b)  add_origin_term(m, out, dim, b, b, -g);
        }
    }

    /* ---- modules (ADR 0012) ----
     *
     * A module's flow is computed from every pathway touching it at once, and
     * the pathways say by NAME which is the energy input and which are the
     * controls (ADR 0013). Nothing here consults the order of the edge array,
     * which is the property the whole design exists for.
     *
     * The flow passes THROUGH: the energy origin is drained, the outgoing
     * pathways are filled, and the module's own row stays empty because a work
     * gate is not a stock. */
    for (i = 0; i < n; i++) {
        const gia_node *nd = &m->nodes[i];
        int    j, t, nt, ae = -1, ctrl = -1, col[2];
        double g = 0.0, W = 0.0, konst = 0.0, coef[2];
        bool   affine = false, subtract = false, drain;

        if (!nd->is_module) continue;

        if (nd->kind == GIA_NODE_EXCHANGE) {
            /* Odum 1972 SecXV Eq (103): J_goods = P J_counter, and currency
             * moves opposite to the goods. The four legs are named, so which
             * component pays and which receives is stated rather than inferred
             * — the leg discovery ADR 0013 said module form makes possible.
             *
             * Linear in the goods leaving the seller, so a transactor leaves
             * the flow matrix CONSTANT and psi exactly zero, unlike a work
             * gate. */
            int    gi = -1, go = -1, ci = -1, co = -1;
            double k  = nd->mod_k;
            double P  = (fabs(nd->mod_price) > GIA_EPS) ? nd->mod_price : 1.0;

            for (j = 0; j < m->n_edges; j++) {
                const gia_edge *e = &m->edges[j];
                if (e->to   == i && e->role == GIA_ROLE_GOODS_IN)    gi = e->from;
                if (e->from == i && e->role == GIA_ROLE_GOODS_OUT)   go = e->to;
                if (e->to   == i && e->role == GIA_ROLE_COUNTER_IN)  ci = e->from;
                if (e->from == i && e->role == GIA_ROLE_COUNTER_OUT) co = e->to;
            }
            if (gi < 0 || go < 0 || ci < 0 || co < 0) continue;

            if (m->nodes[gi].integrates && m->nodes[gi].kind != GIA_NODE_SINK)
                out->a[(size_t)gi*(size_t)dim+(size_t)gi] -= k;
            if (m->nodes[go].integrates)
                out->a[(size_t)go*(size_t)dim+(size_t)gi] += k;
            if (m->nodes[ci].integrates && m->nodes[ci].kind != GIA_NODE_SINK)
                out->a[(size_t)ci*(size_t)dim+(size_t)gi] -= k / P;
            if (m->nodes[co].integrates)
                out->a[(size_t)co*(size_t)dim+(size_t)gi] += k / P;
            continue;
        }

        for (j = 0; j < m->n_edges; j++) {
            const gia_edge *e = &m->edges[j];
            if (e->to == i && e->role == GIA_ROLE_ENERGY)  ae   = e->from;
            if (e->to == i && e->role == GIA_ROLE_CONTROL) ctrl = e->from;
            if (e->from == i && e->role != GIA_ROLE_USED) W += fabs(e->weight);
        }
        if (ae < 0 || W <= 0.0) continue;

        switch (nd->kind) {
            case GIA_NODE_INTERACTION:
                switch (nd->mod_action) {
                    case GIA_ACTION_MULTIPLY:
                        /* F = k Q_energy * prod(Q_control). Every control folds
                         * into the conductance, the n-ary case of the same move
                         * the binary work gate already used (ADR 0012). */
                        g = nd->mod_k;
                        for (j = 0; j < m->n_edges; j++) {
                            const gia_edge *e = &m->edges[j];
                            if (e->to == i && e->role == GIA_ROLE_CONTROL &&
                                e->from >= 0 && q)
                                g *= q[e->from];
                        }
                        break;
                    case GIA_ACTION_DIVIDE: {
                        /* Fig. 2.6(d): linear in Q_energy with a conductance
                         * set by the control. ADR 0002's floor keeps it finite
                         * as the control goes to zero. */
                        double qc = (ctrl >= 0 && q) ? q[ctrl] : 1.0;
                        if (qc < GIA_EPS) qc = GIA_EPS;
                        g = nd->mod_k / qc;
                        break;
                    }
                    case GIA_ACTION_SUBTRACT:
                        /* Fig. 2.6(e): linear in two quantities while open, and
                         * nothing while the clamp holds. The crossing is an
                         * event (module_gap), so each side is exact. */
                        if (ctrl < 0 || !q || q[ae] - q[ctrl] <= 0.0) continue;
                        g        = nd->mod_k;
                        subtract = true;
                        break;
                }
                break;
            case GIA_NODE_LOOP_LIMITED: {
                double qe = q ? q[ae] : 0.0;
                double C  = (nd->mod_capacity > GIA_EPS) ? nd->mod_capacity : 1.0;
                g = nd->mod_k * C / (C + qe);
                break;
            }
            case GIA_NODE_GAIN:
                /* Odum SecIX: the control sets the rate, the energy input
                 * supplies the power and is what gets drained. The term
                 * therefore sits in the CONTROL's column. */
                if (ctrl < 0) continue;
                break;
            case GIA_NODE_SWITCH:
                /* Odum SecXI: a fixed rate while the sensor is above the
                 * threshold. Affine, so it lands in the augmented column. */
                if (ctrl < 0 || !q) continue;
                if (q[ctrl] <= nd->mod_threshold) continue;
                affine = true;
                konst  = nd->mod_k;
                break;
            default:
                continue;
        }

        /* F = sum of coef[t] * Q[col[t]]: the energy input's column for a work
         * gate or a cycling receptor, the control's for an amplifier, the
         * augmented column for an open switch's fixed rate, and two columns
         * for a subtracting action, k Q_energy - k Q_control. Every leg below
         * is a multiple of F, so each is the same entries scaled. */
        nt = 1;
        if (affine)                         { col[0] = dim - 1; coef[0] = konst;     }
        else if (nd->kind == GIA_NODE_GAIN) { col[0] = ctrl;    coef[0] = nd->mod_k; }
        else                                { col[0] = ae;      coef[0] = g;         }
        if (subtract) { col[1] = ctrl; coef[1] = -g; nt = 2; }

        drain = m->nodes[ae].integrates && m->nodes[ae].kind != GIA_NODE_SINK;

        for (t = 0; t < nt; t++) {
            /* Drained once, from the energy input. */
            if (drain) out->a[(size_t)ae*(size_t)dim+(size_t)col[t]] -= coef[t];

            /* Filled per product pathway, in proportion to its weight (ADR 0013
             * decision 4). Equal weights give an even split. */
            for (j = 0; j < m->n_edges; j++) {
                const gia_edge *e = &m->edges[j];
                int    b2;
                if (e->from != i || e->to < 0 || e->role == GIA_ROLE_USED) continue;
                b2 = e->to;
                if (!m->nodes[b2].integrates) continue;
                out->a[(size_t)b2*(size_t)dim+(size_t)col[t]] +=
                    coef[t] * fabs(e->weight) / W;
            }

            /* ADR 0017: each control draws s F from its source, and that leaves
             * on the used leg holding the control's carrier -- dissipated,
             * never added to the product. Load guarantees the used leg. */
            for (j = 0; j < m->n_edges; j++) {
                const gia_edge *e = &m->edges[j];
                int c, k, sink = -1;
                if (e->to != i || e->role != GIA_ROLE_CONTROL ||
                    e->use_ratio <= 0.0 || e->from < 0) continue;
                c = e->from;
                for (k = 0; k < m->n_edges; k++)
                    if (m->edges[k].from == i && m->edges[k].role == GIA_ROLE_USED &&
                        gia_node_carrier(m, m->edges[k].to) == gia_node_carrier(m, c))
                        sink = m->edges[k].to;
                if (m->nodes[c].integrates && m->nodes[c].kind != GIA_NODE_SINK)
                    out->a[(size_t)c*(size_t)dim+(size_t)col[t]] -=
                        e->use_ratio * coef[t];
                if (sink >= 0 && m->nodes[sink].integrates)
                    out->a[(size_t)sink*(size_t)dim+(size_t)col[t]] +=
                        e->use_ratio * coef[t];
            }
        }
    }
    return true;
}

bool gia_flow_matrix_is_constant(const gia_model *m) {
    int i;
    if (!m) return true;
    for (i = 0; i < m->n_edges; i++) {
        /* A driven RATE is the one kind of time dependence that cannot be
         * absorbed. A driven node VALUE can be, which is why the node
         * attachment leaves this true -- see the note above. */
        if (m->edges[i].forcing.kind != GIA_FORCE_NONE) return false;
    }
    /* SecIII pathways are constant-coefficient. What makes the matrix depend on
     * the state is a module's law: a work gate folds its controls into the
     * conductance, a cycling receptor saturates, a switch flips regime. */
    for (i = 0; i < m->n_nodes; i++) {
        if (!m->nodes[i].is_module) continue;
        switch (m->nodes[i].kind) {
            case GIA_NODE_INTERACTION:
            case GIA_NODE_LOOP_LIMITED:
            case GIA_NODE_SWITCH:
                return false;
            default:
                break;
        }
    }
    return true;
}

/* Note on forcing and gia_flow_matrix_is_constant: a node-attached waveform
 * does NOT make the matrix non-constant. Its generator rows are constant, so
 * the augmented system stays time-invariant and the incipient solution stays
 * exact -- psi remains exactly zero.
 *
 * That is a result, not an oversight. Drift is about the coefficient depending
 * on the STATE, not about it depending on time: a driver that is its own
 * generator can be absorbed into A, and the two calculi agree. The case that
 * does drift is Odum's other attachment point, forcing on an EDGE RATE, where
 * the flow is k(t) Q -- bilinear in driver and state, so not absorbable. That
 * attachment is not implemented here. */

/* Does this model contain a pathway that switches on or off? Those are solved
 * piecewise: there is no smooth alpha across a crossing, so persistence of form
 * holds on each side and not through it. */
/* A module whose law switches regime at a located crossing: a switch at its
 * threshold, and a subtracting action where the difference reaches zero. */
static bool module_switches(const gia_node *nd) {
    return nd->is_module &&
           (nd->kind == GIA_NODE_SWITCH ||
            (nd->kind == GIA_NODE_INTERACTION &&
             nd->mod_action == GIA_ACTION_SUBTRACT));
}

static bool has_switching(const gia_model *m) {
    int i;
    for (i = 0; i < m->n_nodes; i++)
        if (module_switches(&m->nodes[i])) return true;
    return false;
}

/* Signed distance to a switching module's boundary: > 0 while it conducts.
 * For a switch, the sensor (its control) against the threshold; for a
 * subtracting action, the energy input against the control. Inputs are found
 * by ROLE -- the same rule the law uses, so the crossing located is the one
 * the law acts on. */
static double module_gap(const gia_model *m, int ni, const double *q) {
    int j, ae = -1, ctrl = -1;
    for (j = 0; j < m->n_edges; j++) {
        const gia_edge *e = &m->edges[j];
        if (e->to != ni || e->from < 0) continue;
        if (e->role == GIA_ROLE_CONTROL) ctrl = e->from;
        if (e->role == GIA_ROLE_ENERGY)  ae   = e->from;
    }
    if (ctrl < 0) return 1.0;
    if (m->nodes[ni].kind == GIA_NODE_INTERACTION)
        return (ae >= 0) ? q[ae] - q[ctrl] : 1.0;
    return q[ctrl] - m->nodes[ni].mod_threshold;
}

/* One smooth advance: Q(t+h) = exp(A h) Q(t), with A frozen at `q`.
 * Uses the augmented (n+1) form so constant and open-threshold pathways, which
 * contribute a rate rather than a conductance, are carried exactly. */
static bool advance_from(const gia_model *m, const double *q, double t0,
                         double h, double *out);

/* Full augmented state at time t: components, then each waveform's own
 * state, then the phantom 1. The waveform states are evaluated at t rather
 * than reset, so a run broken at an event resumes the driver where it left
 * off instead of restarting it. */
static void augment_state(const gia_model *m, const double *q, double t,
                          double *x, int dim) {
    int i, n = m->n_nodes;
    for (i = 0; i < dim; i++) x[i] = 0.0;
    for (i = 0; i < n; i++) x[i] = q[i];
    x[dim - 1] = 1.0;
    for (i = 0; i < n; i++) {
        const gia_forcing *f = &m->nodes[i].forcing;
        int slot = forcing_slot(m, i);
        if (slot < 0) continue;
        switch (f->kind) {
            case GIA_FORCE_SINE:
                x[slot]     = sin(f->rate * t + f->phase);
                x[slot + 1] = cos(f->rate * t + f->phase);
                break;
            case GIA_FORCE_RAMP:        x[slot] = t;                     break;
            case GIA_FORCE_EXPONENTIAL: x[slot] = exp(f->rate * t);      break;
            default: break;
        }
    }
}

/* Build A at `at` (state) and `t_at` (time), then apply exp(A h) to `q`. */
static bool advance_from_at(const gia_model *m, const double *q,
                            const double *at, double t_at, double t0, double h,
                            double *out) {
    gia_matrix A, Ah, E;
    double    *x = NULL;
    int        i, j, n = m->n_nodes, dim;
    bool       ok = false;

    memset(&A, 0, sizeof(A)); memset(&Ah, 0, sizeof(Ah)); memset(&E, 0, sizeof(E));
    dim = n + forcing_extra_count(m) + 1;

    if (!gia_build_flow_matrix(m, at, t_at, &A)) goto done;
    if (!gia_matrix_init(&Ah, dim))       goto done;
    for (i = 0; i < dim * dim; i++) Ah.a[i] = A.a[i] * h;
    if (!gia_matrix_exp(&Ah, &E))         goto done;

    x = (double *)calloc((size_t)dim, sizeof(double));
    if (!x) goto done;
    augment_state(m, q, t0, x, dim);

    for (i = 0; i < n; i++) {
        double acc = 0.0;
        for (j = 0; j < dim; j++) acc += gia_matrix_at(&E, i, j) * x[j];
        out[i] = acc;
    }
    ok = true;
done:
    free(x);
    gia_matrix_free(&A); gia_matrix_free(&Ah); gia_matrix_free(&E);
    return ok;
}

/* One midpoint step.
 *
 * advance_from() freezes A at the interval's midpoint IN TIME but at its start
 * IN STATE. That is first order where A depends on the state. Predicting to
 * the midpoint first and rebuilding A there makes the step second order, which
 * matters once the matrix genuinely varies: without it a rate-forced or
 * interaction model was being solved with A frozen across the whole horizon. */
static bool advance_mid(const gia_model *m, const double *q, double t0,
                        double h, double *out) {
    double *mid;
    bool    ok;
    int     n = m->n_nodes;

    if (h <= 0.0) { memcpy(out, q, (size_t)n * sizeof(double)); return true; }

    mid = (double *)malloc((size_t)n * sizeof(double));
    if (!mid) return false;
    if (!advance_from(m, q, t0, 0.5 * h, mid)) { free(mid); return false; }
    /* Rebuild at the predicted midpoint, then take the whole step from q. */
    ok = advance_from_at(m, q, mid, t0 + 0.5 * h, t0, h, out);
    free(mid);
    return ok;
}

/* Compose a span out of `steps` midpoint steps. One step is exact when A is
 * constant, so the composition costs nothing in that case and is skipped. */
static bool advance_span(const gia_model *m, const double *q, double t0,
                         double h, int steps, double *out) {
    double *cur;
    int     i, n = m->n_nodes;
    bool    ok = true;

    if (steps < 1) steps = 1;
    if (steps == 1) return advance_mid(m, q, t0, h, out);

    cur = (double *)malloc((size_t)n * sizeof(double));
    if (!cur) return false;
    memcpy(cur, q, (size_t)n * sizeof(double));
    for (i = 0; i < steps && ok; i++) {
        ok = advance_mid(m, cur, t0 + (double)i * h / steps, h / steps, out);
        if (ok) memcpy(cur, out, (size_t)n * sizeof(double));
    }
    memcpy(out, cur, (size_t)n * sizeof(double));
    free(cur);
    return ok;
}

static bool advance_from(const gia_model *m, const double *q, double t0,
                         double h, double *out) {
    return advance_from_at(m, q, q, t0 + 0.5 * h, t0, h, out);
}

static bool advance(const gia_model *m, const double *q, double h, double *out) {
    return advance_from(m, q, 0.0, h, out);
}

/* Earliest crossing in (0, span] of any switching module, located by the
 * Illinois variant of false position -- the same method the kernel uses. The
 * regime is frozen during the search, which is what makes the bracket valid:
 * inside one regime the trajectory is a single exponential.
 *
 * Returns the crossing time, or `span` if none. */
static double locate_event(const gia_model *m, const double *q0, double span,
                           double *work) {
    int    i;
    double earliest = span;

    /* A switch's boundary is its sensor crossing the threshold, a subtracting
     * action's is the difference reaching zero. Each is located rather than
     * stepped over: there is no smooth alpha across it, so persistence of form
     * holds on each side and not through it. */
    for (i = 0; i < m->n_nodes; i++) {
        double lo, hi, glo, ghi, mid, gmid;
        int    it, side = 0;

        if (!module_switches(&m->nodes[i])) continue;

        glo = module_gap(m, i, q0);
        if (!advance(m, q0, earliest, work)) continue;
        ghi = module_gap(m, i, work);
        if ((glo > 0.0) == (ghi > 0.0)) continue;

        lo = 0.0; hi = earliest;
        for (it = 0; it < 64; it++) {
            double denom = (ghi - glo);
            if (fabs(denom) < 1e-300) break;
            mid = lo - glo * (hi - lo) / denom;
            if (!(mid > lo && mid < hi)) mid = 0.5 * (lo + hi);
            if (!advance(m, q0, mid, work)) break;
            gmid = module_gap(m, i, work);
            if (fabs(gmid) < 1e-14 || (hi - lo) < 1e-12) { hi = mid; break; }
            if ((gmid > 0.0) == (glo > 0.0)) {
                lo = mid; glo = gmid;
                if (side == -1) ghi *= 0.5;
                side = -1;
            } else {
                hi = mid; ghi = gmid;
                if (side == +1) glo *= 0.5;
                side = +1;
            }
        }
        if (hi < earliest) earliest = hi;
    }
    return earliest;
}

#define GIA_MAX_EVENTS 64

/* Subintervals used when the flow matrix is not constant.
 *
 * When A IS constant the single exponential is the exact answer and this is
 * bypassed entirely, so the cost is only paid where there is something to
 * compose. gia_integration_error() reports what the composition cost in
 * accuracy, by comparing against twice as many. */
#define GIA_SUBSTEPS 64

/* Walk [0, t], breaking at each located crossing. Fills `out`; when
 * `out_events` is non-NULL it receives the number of crossings taken. */
static bool run_to_n(const gia_model *m, double t, int substeps,
                     double *out, int *out_events) {
    int     n = m->n_nodes, i, events = 0;
    double *cur = NULL, *nxt = NULL, *work = NULL;
    double  elapsed = 0.0;
    bool    ok = false;

    cur  = (double *)malloc((size_t)n * sizeof(double));
    nxt  = (double *)malloc((size_t)n * sizeof(double));
    work = (double *)malloc((size_t)n * sizeof(double));
    if (!cur || !nxt || !work) goto done;

    for (i = 0; i < n; i++) cur[i] = m->nodes[i].q0;

    if (t <= 0.0) { memcpy(out, cur, (size_t)n * sizeof(double)); ok = true; goto done; }

    while (elapsed < t && events <= GIA_MAX_EVENTS) {
        double span = t - elapsed;
        double h    = has_switching(m) ? locate_event(m, cur, span, work) : span;

        if (h <= 0.0 || h > span) h = span;
        /* One exponential is exact for a constant A; otherwise compose. */
        if (!advance_span(m, cur, elapsed, h,
                          gia_flow_matrix_is_constant(m) ? 1 : substeps,
                          nxt)) goto done;

        memcpy(cur, nxt, (size_t)n * sizeof(double));
        elapsed += h;

        if (h < span) {
            /* A boundary was reached. Step fractionally past it so the regime
             * is re-evaluated on the far side rather than re-locating the same
             * crossing forever. */
            double eps = (t > 0.0) ? t * 1e-9 : 1e-12;
            if (!advance_from(m, cur, elapsed, eps, nxt)) goto done;
            memcpy(cur, nxt, (size_t)n * sizeof(double));
            elapsed += eps;
            events++;
        }
    }
    memcpy(out, cur, (size_t)n * sizeof(double));
    if (out_events) *out_events = events;
    ok = true;
done:
    free(cur); free(nxt); free(work);
    return ok;
}

static bool run_to(const gia_model *m, double t, double *out, int *out_events) {
    return run_to_n(m, t, GIA_SUBSTEPS, out, out_events);
}

double gia_integration_error(const gia_model *m, double t) {
    double *a = NULL, *b = NULL, worst = 0.0;
    int     i, n;

    if (!m || m->n_nodes <= 0) return 0.0;
    /* Nothing is composed when A is constant: the single exponential IS the
     * answer, so there is no integration error to report. */
    if (gia_flow_matrix_is_constant(m)) return 0.0;

    n = m->n_nodes;
    a = (double *)malloc((size_t)n * sizeof(double));
    b = (double *)malloc((size_t)n * sizeof(double));
    if (!a || !b) { free(a); free(b); return 0.0; }

    if (run_to_n(m, t, GIA_SUBSTEPS,     a, NULL) &&
        run_to_n(m, t, GIA_SUBSTEPS * 2, b, NULL)) {
        for (i = 0; i < n; i++) {
            double d = fabs(b[i] - a[i]);
            if (d > worst) worst = d;
        }
        /* The raw difference between N and 2N steps is not the error of the
         * N-step answer -- it is the error times (1 - 2^-p) for a method of
         * order p. The step here is midpoint, so p = 2 and the factor is 3/4;
         * scaling by 4/3 turns the difference into an estimate of the error
         * actually carried. Without this the reported figure understates the
         * error by a quarter, which for a number whose whole job is to be
         * compared against psi is the wrong direction to be wrong in. */
        worst *= 4.0 / 3.0;
    }
    free(a); free(b);
    return worst;
}

int gia_count_events(const gia_model *m, double t) {
    double *q;
    int     events = 0;
    if (!m || m->n_nodes <= 0) return 0;
    q = (double *)malloc((size_t)m->n_nodes * sizeof(double));
    if (!q) return 0;
    (void)run_to(m, t, q, &events);
    free(q);
    return events;
}

bool gia_network_state(const gia_model *m, double t, double *out,
                       double *out_drift) {
    if (!m || !out || m->n_nodes <= 0) return false;
    if (!run_to(m, t, out, NULL)) return false;

    if (out_drift) {
        /* Incipient: (d~/d~t)^n e^alpha = (alpha^o)^n e^alpha, alpha^o = A.
         * For constant A, d/dt exp(At) = A exp(At) -- identical, so psi is
         * exactly zero. Where A depends on Q it is re-evaluated at the evolved
         * point, and the difference is the drift. */
        if (gia_flow_matrix_is_constant(m)) {
            *out_drift = 0.0;
        } else {
            gia_matrix A0, A1;
            double    *q0;
            double     worst = 0.0;
            int        i, j, n = m->n_nodes;

            memset(&A0, 0, sizeof(A0)); memset(&A1, 0, sizeof(A1));
            q0 = (double *)malloc((size_t)n * sizeof(double));
            if (q0) {
                for (i = 0; i < n; i++) q0[i] = m->nodes[i].q0;
                /* A0 at the start, A1 at the end -- so this captures a matrix that
                 * varies with t as well as one that varies with the state. */
                if (gia_build_flow_matrix(m, q0, 0.0, &A0) &&
                    gia_build_flow_matrix(m, out, t, &A1)) {
                    for (i = 0; i < n; i++) {
                        double di = gia_matrix_at(&A0, i, n);
                        double dt_ = gia_matrix_at(&A1, i, n);
                        for (j = 0; j < n; j++) {
                            di  += gia_matrix_at(&A0, i, j) * out[j];
                            dt_ += gia_matrix_at(&A1, i, j) * out[j];
                        }
                        if (fabs(dt_ - di) > worst) worst = fabs(dt_ - di);
                    }
                }
                gia_matrix_free(&A0); gia_matrix_free(&A1); free(q0);
            }
            *out_drift = worst;
        }
    }
    return true;
}

bool gia_system_is_closed(const gia_model *m) {
    int i;
    if (!m) return true;
    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *ed = &m->edges[i];
        int a = ed->from;
        if (a < 0) continue;
        /* ADR 0014: the boundary is crossed by quantity. A control at
         * use_ratio 0 moves nothing; one drawn from a source or a constant
         * brings quantity in (ADR 0017). A module only passes on what its legs
         * brought -- which is where the question is asked. */
        if (ed->role == GIA_ROLE_CONTROL && ed->use_ratio <= 0.0) continue;
        if (m->nodes[a].is_module)        continue;
        if (!m->nodes[a].integrates) return false;  /* forced across a boundary */
    }
    return true;
}

/* Carriers are held by components, and a module is not one -- it is a hyperedge
 * drawn as a symbol, holding nothing. Counting its blank carrier would invent a
 * third carrier class for a two-carrier transaction, and then a residual would
 * be reported for a class no stock belongs to. */
static bool node_holds_carrier(const gia_model *m, int i) {
    return !m->nodes[i].is_module;
}

int gia_carrier_count(const gia_model *m) {
    int i, j, n = 0;
    if (!m) return 1;
    for (i = 0; i < m->n_nodes; i++) {
        const char *c = m->nodes[i].carrier ? m->nodes[i].carrier : "";
        int seen = 0;
        if (!node_holds_carrier(m, i)) continue;
        for (j = 0; j < i; j++) {
            const char *d = m->nodes[j].carrier ? m->nodes[j].carrier : "";
            if (!node_holds_carrier(m, j)) continue;
            if (!strcmp(c, d)) { seen = 1; break; }
        }
        if (!seen) n++;
    }
    return n > 0 ? n : 1;
}

const char *gia_carrier_name(const gia_model *m, int idx) {
    int i, j, n = 0;
    if (!m || idx < 0) return "";
    for (i = 0; i < m->n_nodes; i++) {
        const char *c = m->nodes[i].carrier ? m->nodes[i].carrier : "";
        int seen = 0;
        if (!node_holds_carrier(m, i)) continue;
        for (j = 0; j < i; j++) {
            const char *d = m->nodes[j].carrier ? m->nodes[j].carrier : "";
            if (!node_holds_carrier(m, j)) continue;
            if (!strcmp(c, d)) { seen = 1; break; }
        }
        if (seen) continue;
        if (n == idx) return c;
        n++;
    }
    return "";
}

int gia_node_carrier(const gia_model *m, int node_idx) {
    int k, n;
    const char *c;
    if (!m || node_idx < 0 || node_idx >= m->n_nodes) return 0;
    if (!node_holds_carrier(m, node_idx)) return -1;     /* holds nothing */
    c = m->nodes[node_idx].carrier ? m->nodes[node_idx].carrier : "";
    n = gia_carrier_count(m);
    for (k = 0; k < n; k++)
        if (!strcmp(c, gia_carrier_name(m, k))) return k;
    return 0;
}

double gia_conservation_residual_for(const gia_model *m, double t, int carrier) {
    double *q, s0 = 0.0, st = 0.0;
    int     i;
    if (!m || m->n_nodes <= 0) return 0.0;
    q = (double *)malloc((size_t)m->n_nodes * sizeof(double));
    if (!q) return 0.0;
    if (!gia_network_state(m, t, q, NULL)) { free(q); return 0.0; }
    for (i = 0; i < m->n_nodes; i++) {
        if (!m->nodes[i].integrates) continue;
        if (gia_node_carrier(m, i) != carrier) continue;
        s0 += m->nodes[i].q0;
        st += q[i];
    }
    free(q);
    return fabs(st - s0);
}

double gia_conservation_residual(const gia_model *m, double t) {
    int    k, n;
    double worst = 0.0;
    if (!m) return 0.0;
    /* The maximum, not the sum. A sum would let a surplus in one carrier
     * cancel a deficit in another and report zero -- which is the defect
     * carriers exist to remove, not one to carry forward. */
    n = gia_carrier_count(m);
    for (k = 0; k < n; k++) {
        double r = gia_conservation_residual_for(m, t, k);
        if (r > worst) worst = r;
    }
    return worst;
}

/* ================================================================== *
 * 5. Model loading
 * ================================================================== */

const char *gia_node_kind_name(gia_node_kind k) {
    switch (k) {
        case GIA_NODE_SOURCE:       return "source";
        case GIA_NODE_STORAGE:      return "storage";
        case GIA_NODE_SINK:         return "sink";
        case GIA_NODE_CONSTANT:     return "constant";
        case GIA_NODE_INTERACTION:  return "interaction";
        case GIA_NODE_GAIN:         return "gain";
        case GIA_NODE_LOOP_LIMITED: return "loop_limited";
        case GIA_NODE_SWITCH:       return "switch";
        case GIA_NODE_EXCHANGE:     return "exchange";
        case GIA_NODE_CONSUMER:     return "consumer";
        default:                    return "unknown";
    }
}

static gia_node_kind kind_of(const char *type) {
    if (!type)                              return GIA_NODE_UNKNOWN;
    if (!strcmp(type, "source"))            return GIA_NODE_SOURCE;
    if (!strcmp(type, "storage"))           return GIA_NODE_STORAGE;
    if (!strcmp(type, "store"))             return GIA_NODE_STORAGE;
    if (!strcmp(type, "sink"))              return GIA_NODE_SINK;
    if (!strcmp(type, "constant"))          return GIA_NODE_CONSTANT;
    if (!strcmp(type, "interaction"))       return GIA_NODE_INTERACTION;
    if (!strcmp(type, "gain"))              return GIA_NODE_GAIN;
    if (!strcmp(type, "loop_limited"))      return GIA_NODE_LOOP_LIMITED;
    if (!strcmp(type, "switch"))            return GIA_NODE_SWITCH;
    if (!strcmp(type, "exchange"))          return GIA_NODE_EXCHANGE;
    return GIA_NODE_UNKNOWN;
}

/* Types this engine deliberately refuses rather than accepts as a label.
 *
 * `switch` (Odum 1972 SecXI) needs event location: the flow is discontinuous, so
 * there is no smooth alpha across a crossing and the incipient form holds only
 * piecewise. `exchange` (SecXV) needs two carriers coupled by price, and this
 * engine carries one quantity per component.
 *
 * Accepting either and quietly giving it storage semantics would produce the
 * exact defect this engine was criticised for: a vocabulary that looks
 * meaningful and does nothing. See docs/odum_1972_conformance.md. */
static bool kind_is_refused(const char *type, const char **why) {
    if (!type) return false;
    if (!strcmp(type, "producer")  || !strcmp(type, "consumer") ||
        !strcmp(type, "misc_box")  || !strcmp(type, "system_frame")) {
        *why = "composite types are not expanded here (see ADR 0010)";
        return true;
    }
    return false;
}

static double num_field(const cJSON *obj, const char *key, double fallback) {
    const cJSON *it = cJSON_GetObjectItemCaseSensitive(obj, key);
    return cJSON_IsNumber(it) ? it->valuedouble : fallback;
}

static const char *str_field(const cJSON *obj, const char *key,
                             const char *fallback) {
    const cJSON *it = cJSON_GetObjectItemCaseSensitive(obj, key);
    return (cJSON_IsString(it) && it->valuestring) ? it->valuestring : fallback;
}

/* Put a node into exponential form f(t) = e^phi(t).
 *
 * The assignment is not arbitrary, and the degree of phi is the whole of it:
 *
 *   affine phi (degree <= 1)  -- constant-coefficient process. IDC and TDC
 *       agree identically, at every order. Sources, storages and consumers
 *       land here: a steady inflow, an exponential store and a metabolic
 *       drain are all linear in their own rate.
 *
 *   quadratic phi (degree 2)  -- genuinely variable-coefficient. The Bell
 *       recursion picks up phi'' and the traditional derivative drifts away
 *       from the incipient one. Interactions and emergent regulators land
 *       here, which is the substantive claim: in an Odum graph the drift is
 *       not spread evenly, it is concentrated exactly on the nodes that
 *       transform rather than store, and those are the generative ones.
 */
static void phi_for_node(gia_node *nd, const cJSON *jn) {
    memset(&nd->phi, 0, sizeof(nd->phi));

    switch (nd->kind) {
        case GIA_NODE_SOURCE: {
            double v = num_field(jn, "initial_value", 1.0);
            nd->phi.c[1] = v;
            nd->phi.degree = 1;
            break;
        }
        case GIA_NODE_INTERACTION: {
            double g = num_field(jn, "generativity_factor", 1.0);
            nd->phi.c[2] = 0.5 * g;      /* phi = g t^2 / 2 */
            nd->phi.degree = 2;
            break;
        }
        case GIA_NODE_STORAGE: {
            double level = num_field(jn, "current_level", 1.0);
            double cap   = num_field(jn, "capacity", 0.0);
            double rate  = (cap > GIA_EPS) ? (level / cap) : 0.0;
            if (level < GIA_EPS) level = GIA_EPS;
            nd->phi.c[0] = log(level);   /* so that f(0) = current_level */
            nd->phi.c[1] = rate;         /* fractional fullness as the rate  */
            nd->phi.degree = 1;
            break;
        }
        case GIA_NODE_CONSUMER: {
            double mr = num_field(jn, "metabolic_rate", 0.0);
            nd->phi.c[1] = -mr;
            nd->phi.degree = 1;
            break;
        }
        case GIA_NODE_GAIN: {
            /* Odum SecIX: output is proportional to the control signal, so the
             * amplifier's own form is linear in its gain. */
            double k = num_field(jn, "k", num_field(jn, "gain", 1.0));
            nd->phi.c[1] = k;
            nd->phi.degree = 1;
            break;
        }
        case GIA_NODE_EXCHANGE:
        case GIA_NODE_SWITCH: {
            /* Odum SecXI: on or off. Its own analytic form is a level; the
             * switching itself lives on the threshold pathways it gates. */
            nd->phi.c[0] = log(fmax(num_field(jn, "value", 1.0), GIA_EPS));
            nd->phi.degree = 0;
            break;
        }
        case GIA_NODE_LOOP_LIMITED: {
            /* Odum SecXIII: saturating, so the exponent carries curvature. */
            double k = num_field(jn, "k", 1.0);
            nd->phi.c[2] = 0.5 * k;
            nd->phi.degree = 2;
            break;
        }
        case GIA_NODE_SINK:
        case GIA_NODE_CONSTANT: {
            nd->phi.c[0] = log(fmax(num_field(jn, "value", 1.0), GIA_EPS));
            nd->phi.degree = 0;
            break;
        }
        default: {
            nd->phi.c[1] = num_field(jn, "initial_value", 0.0);
            nd->phi.degree = 1;
            break;
        }
    }
}

/* Map the seed's pathway label onto an Odum law. R2.3: a `flow_type` must mean
 * something or be rejected -- an inert vocabulary that looks meaningful is
 * worse than none. `logic` is preferred when present; `flow_type` is the
 * fallback for the descriptive labels the seed format uses. */
static bool logic_of(const char *s, gia_logic *out) {
    if (!s) return false;
    if (!strcmp(s, "linear"))                 { *out = GIA_LOGIC_LINEAR;      return true; }
    if (!strcmp(s, "reversible"))             { *out = GIA_LOGIC_REVERSIBLE;  return true; }
    if (!strcmp(s, "constant"))               { *out = GIA_LOGIC_CONSTANT;    return true; }
    /* Descriptive labels from the Odum-shaped seed format. */
    if (!strcmp(s, "inflow"))                 { *out = GIA_LOGIC_LINEAR;      return true; }
    if (!strcmp(s, "outflow"))                { *out = GIA_LOGIC_LINEAR;      return true; }
    if (!strcmp(s, "flow"))                   { *out = GIA_LOGIC_LINEAR;      return true; }
    if (!strcmp(s, "ordinal_ascent"))         { *out = GIA_LOGIC_LINEAR;      return true; }
    if (!strcmp(s, "emergent_feedback_loop")) { *out = GIA_LOGIC_LINEAR;      return true; }
    if (!strcmp(s, "diffusion"))              { *out = GIA_LOGIC_REVERSIBLE;  return true; }
    return false;
}

/* The module a law removed from pathways now lives on (ADR 0012 decision 5),
 * or NULL. Named in the load error, so a model written the old way says what
 * to write instead -- the loader never rewrites it, because a gate the author
 * did not draw would vote on ordinality. The two seed labels named the work
 * gate law. */
static const char *module_for_law(const char *s) {
    if (!s) return NULL;
    if (!strcmp(s, "interaction") || !strcmp(s, "generative_production") ||
        !strcmp(s, "ordinal_feedback"))
        return "an `interaction` module";
    if (!strcmp(s, "ratio"))     return "an `interaction` module with action divide";
    if (!strcmp(s, "subtract"))  return "an `interaction` module with action subtract";
    if (!strcmp(s, "limit"))     return "a `loop_limited` module";
    if (!strcmp(s, "gain"))      return "a `gain` module";
    if (!strcmp(s, "threshold")) return "a `switch` module";
    if (!strcmp(s, "exchange"))  return "an `exchange` module with four named legs";
    return NULL;
}

static int find_node(const gia_model *m, const char *id) {
    int i;
    if (!id) return -1;
    for (i = 0; i < m->n_nodes; i++)
        if (m->nodes[i].id && !strcmp(m->nodes[i].id, id)) return i;
    return -1;
}

bool gia_model_load(gia_model *m, cJSON *root) {
    cJSON *jnodes, *jedges, *jparams;
    int    i, n, e;

    if (!m || !root) return false;
    memset(m, 0, sizeof(*m));
    m->root = root;

    jnodes = cJSON_GetObjectItemCaseSensitive(root, "nodes");
    jedges = cJSON_GetObjectItemCaseSensitive(root, "edges");
    if (!cJSON_IsArray(jnodes)) {
        fprintf(stderr, "engine: model has no \"nodes\" array\n");
        return false;
    }

    n = cJSON_GetArraySize(jnodes);
    if (n <= 0) {
        fprintf(stderr, "engine: model has no nodes\n");
        return false;
    }

    m->system_name = str_field(root, "system_name", "(unnamed system)");
    m->nodes = (gia_node *)calloc((size_t)n, sizeof(gia_node));
    if (!m->nodes) return false;
    m->n_nodes = n;

    for (i = 0; i < n; i++) {
        const cJSON *jn = cJSON_GetArrayItem(jnodes, i);
        gia_node    *nd = &m->nodes[i];
        const char *ty  = str_field(jn, "type", NULL);
        const char *why = NULL;

        nd->id    = str_field(jn, "id", NULL);
        nd->label = str_field(jn, "label", nd->id ? nd->id : "(unlabelled)");

        if (kind_is_refused(ty, &why)) {
            fprintf(stderr, "engine: node %d ('%s'): type '%s' is not "
                            "implemented here -- %s\n",
                    i, nd->id ? nd->id : "?", ty, why);
            gia_model_free(m);
            return false;
        }
        nd->kind = kind_of(ty);
        if (nd->kind == GIA_NODE_UNKNOWN) {
            fprintf(stderr, "engine: node %d ('%s'): unknown type '%s'\n",
                    i, nd->id ? nd->id : "?", ty ? ty : "(missing)");
            gia_model_free(m);
            return false;
        }
        phi_for_node(nd, jn);
        nd->quality_input = num_field(jn, "quality_input", 0.0);
        {   /* A module hosts its law over the pathways around it (ADR 0012).
             * Its type and its `module` block go together: a work gate without
             * one would be a storage wearing the name, holding quantity and
             * counting toward ordinality as a component, which is exactly what
             * the example seeds did before they were migrated. */
            const cJSON *mj = cJSON_GetObjectItemCaseSensitive(jn, "module");
            nd->is_module = false;
            if (cJSON_IsObject(mj)) {
                switch (nd->kind) {
                    case GIA_NODE_INTERACTION:
                    case GIA_NODE_GAIN:
                    case GIA_NODE_SWITCH:
                    case GIA_NODE_LOOP_LIMITED:
                    case GIA_NODE_EXCHANGE:
                        break;
                    default:
                        fprintf(stderr,
                                "engine: node %d ('%s'): type '%s' hosts no "
                                "law, so it cannot carry a `module` block "
                                "(ADR 0012)\n", i, nd->id ? nd->id : "?",
                                gia_node_kind_name(nd->kind));
                        gia_model_free(m);
                        return false;
                }
                nd->is_module     = true;
                nd->mod_k         = num_field(mj, "k", 1.0);
                nd->mod_capacity  = num_field(mj, "capacity", 1.0);
                nd->mod_threshold = num_field(mj, "threshold", 0.0);
                {   /* ADR 0016: the interaction glyph's arithmetic. Absent
                     * means multiply, so every work gate written before this
                     * means what it meant. */
                    const cJSON *aj = cJSON_GetObjectItemCaseSensitive(mj, "action");
                    nd->mod_action = GIA_ACTION_MULTIPLY;
                    if (aj && nd->kind != GIA_NODE_INTERACTION) {
                        fprintf(stderr, "engine: module '%s' (%s): `action` "
                                        "belongs to an interaction module "
                                        "(ADR 0016)\n", nd->id ? nd->id : "?",
                                gia_node_kind_name(nd->kind));
                        gia_model_free(m);
                        return false;
                    }
                    if (aj) {
                        const char *a = cJSON_IsString(aj) ? aj->valuestring : NULL;
                        if      (a && !strcmp(a, "multiply")) nd->mod_action = GIA_ACTION_MULTIPLY;
                        else if (a && !strcmp(a, "divide"))   nd->mod_action = GIA_ACTION_DIVIDE;
                        else if (a && !strcmp(a, "subtract")) nd->mod_action = GIA_ACTION_SUBTRACT;
                        else {
                            fprintf(stderr, "engine: module '%s': unknown action "
                                            "(expected multiply, divide or "
                                            "subtract, Odum Fig. 2.6)\n",
                                    nd->id ? nd->id : "?");
                            gia_model_free(m);
                            return false;
                        }
                    }
                }
                /* Neutral name first, the money-specific one as an alias — the
                 * same order PR 3 established, because naming it price is what
                 * made an earlier revision rule out barter. */
                nd->mod_price     = num_field(mj, "exchange_ratio",
                                    num_field(mj, "price", 1.0));
            }
            if (!nd->is_module) {
                switch (nd->kind) {
                    case GIA_NODE_INTERACTION:
                    case GIA_NODE_GAIN:
                    case GIA_NODE_SWITCH:
                    case GIA_NODE_LOOP_LIMITED:
                    case GIA_NODE_EXCHANGE:
                        fprintf(stderr,
                                "engine: node %d ('%s'): type '%s' is a module "
                                "and needs its `module` block (ADR 0012)\n",
                                i, nd->id ? nd->id : "?",
                                gia_node_kind_name(nd->kind));
                        gia_model_free(m);
                        return false;
                    default:
                        break;
                }
            }
        }
        nd->carrier       = str_field(jn, "carrier", "");
        {   /* ADR 0006: the waveform is attached to the element, not held in a
             * model-root block keyed by id -- so a reader of the nodes array
             * can see that a component is driven, and deleting the component
             * cannot orphan its forcing. */
            const cJSON *fj = cJSON_GetObjectItemCaseSensitive(jn, "forcing");
            memset(&nd->forcing, 0, sizeof(nd->forcing));
            nd->forcing.kind = GIA_FORCE_NONE;
            if (cJSON_IsObject(fj)) {
                const char *kind = str_field(fj, "kind", "none");
                if      (!strcmp(kind, "none"))  nd->forcing.kind = GIA_FORCE_NONE;
                else if (!strcmp(kind, "sine"))  nd->forcing.kind = GIA_FORCE_SINE;
                else if (!strcmp(kind, "ramp"))  nd->forcing.kind = GIA_FORCE_RAMP;
                else if (!strcmp(kind, "exponential"))
                                                 nd->forcing.kind = GIA_FORCE_EXPONENTIAL;
                else {
                    /* Refused rather than approximated. A square wave, a
                     * sawtooth or jitter is not its own generator, so it cannot
                     * be carried as state, and the closed form would quietly
                     * become a step-and-hope. */
                    fprintf(stderr,
                            "engine: node %d ('%s'): forcing kind '%s' is not "
                            "implemented here -- only waveforms that generate "
                            "themselves (sine, ramp, exponential) can be "
                            "carried as state and keep the solution exact\n",
                            i, nd->id ? nd->id : "?", kind);
                    gia_model_free(m);
                    return false;
                }
                nd->forcing.amplitude = num_field(fj, "amplitude", 1.0);
                nd->forcing.rate      = num_field(fj, "rate",
                                        num_field(fj, "frequency", 1.0));
                nd->forcing.phase     = num_field(fj, "phase", 0.0);
                nd->forcing.offset    = num_field(fj, "offset", 0.0);
            }
        }

        /* Live quantity for the network solution, and whether it is solved for.
         * Odum holds a source at its value rather than integrating it (1972
         * SecII), and the same is true of a constant. */
        switch (nd->kind) {
            case GIA_NODE_STORAGE:
                nd->q0 = num_field(jn, "current_level",
                                   num_field(jn, "value", 0.0));
                nd->integrates = true;
                break;
            case GIA_NODE_SOURCE:
            case GIA_NODE_CONSTANT:
                /* Odum SecII: a source is a forcing function. Its quantity is
                 * held, never integrated, so a pathway leaving it is not
                 * depleting it. Same for a constant. */
                nd->q0 = num_field(jn, "initial_value",
                                   num_field(jn, "value", 1.0));
                nd->integrates = false;
                break;
            case GIA_NODE_SINK:
                /* Odum SecV: absorbs used energy and is never drained. */
                nd->q0 = num_field(jn, "value", 0.0);
                nd->integrates = true;
                break;
            default:
                nd->q0 = num_field(jn, "value",
                                   num_field(jn, "initial_value", 0.0));
                nd->integrates = true;
                break;
        }
        /* Checked here, after the switch above has set `integrates`. Reading
         * it earlier gave a zeroed field, so the check silently never fired. */
        if (nd->forcing.kind != GIA_FORCE_NONE && nd->integrates) {
            fprintf(stderr,
                    "engine: node %d ('%s'): forcing drives a HELD value "
                    "(Odum's X or N), but this component integrates; attach it "
                    "to a source or constant\n", i, nd->id ? nd->id : "?");
            gia_model_free(m);
            return false;
        }
        if (!nd->id) {
            fprintf(stderr, "engine: node %d has no \"id\"\n", i);
            gia_model_free(m);
            return false;
        }
    }

    e = cJSON_IsArray(jedges) ? cJSON_GetArraySize(jedges) : 0;
    if (e > 0) {
        m->edges = (gia_edge *)calloc((size_t)e, sizeof(gia_edge));
        if (!m->edges) { gia_model_free(m); return false; }
        m->n_edges = e;
        for (i = 0; i < e; i++) {
            const cJSON *je = cJSON_GetArrayItem(jedges, i);
            gia_edge    *ed = &m->edges[i];
            const char *lg;
            ed->from      = find_node(m, str_field(je, "source", NULL));
            ed->to        = find_node(m, str_field(je, "target", NULL));
            ed->flow_type = str_field(je, "flow_type", "flow");
            ed->weight    = num_field(je, "weight", 1.0);
            {   /* ADR 0006's second attachment point: the same waveform
                 * vocabulary, driving this pathway's rate rather than a
                 * component's held value. */
                const cJSON *fj = cJSON_GetObjectItemCaseSensitive(je, "forcing");
                memset(&ed->forcing, 0, sizeof(ed->forcing));
                ed->forcing.kind = GIA_FORCE_NONE;
                if (cJSON_IsObject(fj)) {
                    const char *kind = str_field(fj, "kind", "none");
                    if      (!strcmp(kind, "none")) ed->forcing.kind = GIA_FORCE_NONE;
                    else if (!strcmp(kind, "sine")) ed->forcing.kind = GIA_FORCE_SINE;
                    else if (!strcmp(kind, "ramp")) ed->forcing.kind = GIA_FORCE_RAMP;
                    else if (!strcmp(kind, "exponential"))
                                                    ed->forcing.kind = GIA_FORCE_EXPONENTIAL;
                    else {
                        fprintf(stderr, "engine: edge %d: forcing kind '%s' is "
                                        "not implemented\n", i, kind);
                        gia_model_free(m);
                        return false;
                    }
                    ed->forcing.amplitude = num_field(fj, "amplitude", 1.0);
                    ed->forcing.rate      = num_field(fj, "rate",
                                            num_field(fj, "frequency", 1.0));
                    ed->forcing.phase     = num_field(fj, "phase", 0.0);
                    ed->forcing.offset    = num_field(fj, "offset", 0.0);
                }
            }
            {   /* docs/emergy_synthesis.md 8: partition is the default, because
                 * a split is the ordinary case and co-production is the claim. */
                const char *om = str_field(je, "output_mode", "partition");
                if (!strcmp(om, "replicate"))      ed->out_mode = GIA_OUT_REPLICATE;
                else if (!strcmp(om, "partition")) ed->out_mode = GIA_OUT_PARTITION;
                else {
                    fprintf(stderr, "engine: edge %d: unknown output_mode '%s' "
                                    "(expected 'partition' or 'replicate')\n", i, om);
                    gia_model_free(m);
                    return false;
                }
            }

            {   /* ADR 0013: a pathway entering a module says what it is TO
                 * that module, by name. One spelling, `role`. */
                const char *r = str_field(je, "role", NULL);
                ed->role = GIA_ROLE_NONE;
                if (r) {
                    if      (!strcmp(r, "energy"))      ed->role = GIA_ROLE_ENERGY;
                    else if (!strcmp(r, "control"))     ed->role = GIA_ROLE_CONTROL;
                    else if (!strcmp(r, "used"))        ed->role = GIA_ROLE_USED;
                    else if (!strcmp(r, "goods_in"))    ed->role = GIA_ROLE_GOODS_IN;
                    else if (!strcmp(r, "goods_out"))   ed->role = GIA_ROLE_GOODS_OUT;
                    else if (!strcmp(r, "counter_in"))  ed->role = GIA_ROLE_COUNTER_IN;
                    else if (!strcmp(r, "counter_out")) ed->role = GIA_ROLE_COUNTER_OUT;
                    else {
                        fprintf(stderr, "engine: edge %d: unknown role '%s' "
                                        "(expected energy, control, used, "
                                        "goods_in, goods_out, counter_in or "
                                        "counter_out)\n", i, r);
                        gia_model_free(m);
                        return false;
                    }
                }
            }

            {   /* ADR 0017 decision 4: a control leg says how much it draws.
                 * Required, because the defect it replaces was a silent zero;
                 * 0 is allowed and means the control's energy is neglected.
                 * Anywhere else the field would be read by nothing. */
                const cJSON *ur = cJSON_GetObjectItemCaseSensitive(je, "use_ratio");
                if (ur && ed->role != GIA_ROLE_CONTROL) {
                    fprintf(stderr, "engine: edge %d: use_ratio belongs on a "
                                    "control leg, and this one is %s\n",
                            i, ed->role == GIA_ROLE_NONE
                                   ? "not a module input"
                                   : gia_role_name(ed->role));
                    gia_model_free(m);
                    return false;
                }
                if (ed->role == GIA_ROLE_CONTROL) {
                    if (!ur) {
                        fprintf(stderr, "engine: edge %d: a control leg must "
                                        "state its use_ratio, the share of the "
                                        "module's flow it draws (0 to neglect "
                                        "it; ADR 0017)\n", i);
                        gia_model_free(m);
                        return false;
                    }
                    if (!cJSON_IsNumber(ur) || !(ur->valuedouble >= 0.0) ||
                        !isfinite(ur->valuedouble)) {
                        fprintf(stderr, "engine: edge %d: use_ratio must be a "
                                        "number >= 0\n", i);
                        gia_model_free(m);
                        return false;
                    }
                    ed->use_ratio = ur->valuedouble;
                }
            }

            lg = str_field(je, "logic", ed->flow_type);
            if (!logic_of(lg, &ed->logic)) {
                if (module_for_law(lg))
                    fprintf(stderr, "engine: edge %d: '%s' is a module law and "
                                    "no longer lives on a pathway; write %s with "
                                    "named roles (ADR 0012 decision 5)\n",
                            i, lg, module_for_law(lg));
                else
                    fprintf(stderr, "engine: edge %d: unknown pathway law '%s' "
                                    "(a pathway is linear, reversible or "
                                    "constant)\n", i, lg);
                gia_model_free(m);
                return false;
            }

            if (ed->from < 0 || ed->to < 0)
                fprintf(stderr,
                        "engine: warning: edge %d references an unknown node; "
                        "it will not contribute to ordinality\n", i);
        }
    }

    /* A pathway may not cross carriers. Grain does not become money by flowing
     * along an edge -- it is exchanged for money, and that coupling is what the
     * transaction diamond is for. Anything else moving quantity between
     * incommensurable stocks is a modelling error, so it is named rather than
     * quietly integrated. */
    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *ed = &m->edges[i];
        if (ed->from < 0 || ed->to < 0) continue;
        /* A module has no carrier of its own -- it is a hyperedge, not a stock,
         * so comparing a leg's carrier against the module's says nothing. Which
         * of its legs must agree is the module's own business, checked per
         * module below: a transactor's two pairs each move one carrier, while
         * every other module matches its energy input to its outputs and each
         * drawn control's carrier to a used leg (ADR 0017). */
        if (m->nodes[ed->from].is_module || m->nodes[ed->to].is_module) continue;
        if (gia_node_carrier(m, ed->from) != gia_node_carrier(m, ed->to)) {
            fprintf(stderr,
                    "engine: edge %d ('%s' -> '%s'): carrier '%s' cannot flow "
                    "into carrier '%s'; couple two carriers through an "
                    "`exchange` module (Odum 1972 SecXV)\n",
                    i, m->nodes[ed->from].id, m->nodes[ed->to].id,
                    gia_carrier_name(m, gia_node_carrier(m, ed->from)),
                    gia_carrier_name(m, gia_node_carrier(m, ed->to)));
            gia_model_free(m);
            return false;
        }
    }

    /* A role only means something on a pathway TOUCHING a module, and on the
     * side its name implies: an input role must enter one, an output role must
     * leave one. Anywhere else it is a misplaced field, and ignoring it
     * silently would let a modeller believe a leg was marked when nothing read
     * the mark. */
    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *ed = &m->edges[i];
        int  side;
        bool ok_side;

        if (ed->role == GIA_ROLE_NONE) continue;
        side    = role_is_outgoing(ed->role) ? ed->from : ed->to;
        ok_side = (side >= 0) && m->nodes[side].is_module;

        if (!ok_side) {
            fprintf(stderr,
                    "engine: edge %d: role '%s' must name a pathway %s a "
                    "module, and this one does not (ADR 0013)\n",
                    i, gia_role_name(ed->role),
                    role_is_outgoing(ed->role) ? "leaving" : "entering");
            gia_model_free(m);
            return false;
        }
    }

    /* Each module states the roles it requires (ADR 0013). Anything else is an
     * error naming the module — never a silent default, and never a fallback to
     * position. */
    for (i = 0; i < m->n_nodes; i++) {
        const gia_node *nd = &m->nodes[i];
        int j, n_energy = 0, n_control = 0, n_out = 0, want_control = -1;

        if (!nd->is_module) continue;

        if (nd->kind == GIA_NODE_EXCHANGE) {
            int gi = -1, go = -1, ci = -1, co = -1, dup = 0;

            for (j = 0; j < m->n_edges; j++) {
                const gia_edge *ed = &m->edges[j];
                int touches_in  = (ed->to   == i);
                int touches_out = (ed->from == i);
                if (!touches_in && !touches_out) continue;
                if (ed->role == GIA_ROLE_NONE) {
                    fprintf(stderr, "engine: edge %d touching transactor '%s' "
                                    "declares no role; its four legs are named, "
                                    "never ordered (ADR 0013)\n", j, nd->id);
                    gia_model_free(m);
                    return false;
                }
                switch (ed->role) {
                    case GIA_ROLE_GOODS_IN:    if (gi >= 0) dup = 1; gi = ed->from; break;
                    case GIA_ROLE_GOODS_OUT:   if (go >= 0) dup = 1; go = ed->to;   break;
                    case GIA_ROLE_COUNTER_IN:  if (ci >= 0) dup = 1; ci = ed->from; break;
                    case GIA_ROLE_COUNTER_OUT: if (co >= 0) dup = 1; co = ed->to;   break;
                    default:
                        fprintf(stderr, "engine: edge %d: role '%s' has no "
                                        "meaning for a transactor; its legs are "
                                        "goods_in, goods_out, counter_in and "
                                        "counter_out\n",
                                j, gia_role_name(ed->role));
                        gia_model_free(m);
                        return false;
                }
            }
            if (dup || gi < 0 || go < 0 || ci < 0 || co < 0) {
                fprintf(stderr, "engine: transactor '%s' needs exactly one each "
                                "of goods_in, goods_out, counter_in and "
                                "counter_out\n", nd->id);
                gia_model_free(m);
                return false;
            }
            if (ci == co) {
                fprintf(stderr, "engine: transactor '%s' pays '%s' from itself, "
                                "so nothing moves\n", nd->id, m->nodes[ci].id);
                gia_model_free(m);
                return false;
            }
            /* One leg pair moves one kind of thing. Barter — the goods and
             * counter carriers being the same — stays permitted; what is not
             * is a single pair straddling two carriers. */
            if (gia_node_carrier(m, gi) != gia_node_carrier(m, go) ||
                gia_node_carrier(m, ci) != gia_node_carrier(m, co)) {
                fprintf(stderr, "engine: transactor '%s': each leg pair must "
                                "move one carrier\n", nd->id);
                gia_model_free(m);
                return false;
            }
            continue;                    /* the checks below are for the rest */
        }

        for (j = 0; j < m->n_edges; j++) {
            const gia_edge *ed = &m->edges[j];
            if (ed->to == i) {
                if (ed->role == GIA_ROLE_NONE) {
                    fprintf(stderr,
                            "engine: edge %d entering module '%s' declares no "
                            "role; a module's inputs are named, never ordered "
                            "(ADR 0013)\n", j, nd->id);
                    gia_model_free(m);
                    return false;
                }
                if (ed->role == GIA_ROLE_ENERGY)  n_energy++;
                if (ed->role == GIA_ROLE_CONTROL) n_control++;
            }
            if (ed->from == i && ed->role != GIA_ROLE_USED) n_out++;
        }

        switch (nd->kind) {
            case GIA_NODE_INTERACTION:
                /* A quotient or a difference of three things is not in the
                 * figure, so divide and subtract take exactly one (ADR 0016). */
                want_control = (nd->mod_action == GIA_ACTION_MULTIPLY) ? -1 : 1;
                break;
            case GIA_NODE_GAIN:         want_control =  1; break;
            case GIA_NODE_SWITCH:       want_control =  1; break;
            case GIA_NODE_LOOP_LIMITED: want_control =  0; break;
            default: break;
        }

        if (n_energy != 1) {
            fprintf(stderr, "engine: module '%s' (%s) needs exactly one energy "
                            "input, the one it consumes; it has %d\n",
                    nd->id, gia_node_kind_name(nd->kind), n_energy);
            gia_model_free(m);
            return false;
        }
        if (want_control >= 0 && n_control != want_control) {
            /* A cycling receptor with a surplus input is the case GSSK
             * silently discards. Named here instead. */
            fprintf(stderr, "engine: module '%s' (%s) needs exactly %d control "
                            "input%s; it has %d\n",
                    nd->id, gia_node_kind_name(nd->kind), want_control,
                    want_control == 1 ? "" : "s", n_control);
            gia_model_free(m);
            return false;
        }
        if (n_out < 1) {
            fprintf(stderr, "engine: module '%s' has no outgoing pathway, so "
                            "its output goes nowhere\n", nd->id);
            gia_model_free(m);
            return false;
        }

        /* ADR 0017 decision 5: what the controls are drawn of leaves on a
         * `used` leg, one per carrier, to a sink holding that carrier. Every
         * drawn carrier needs one, and a used leg must have something to
         * carry. A price is the case this catches: money is not dissipated, so
         * a nonzero use_ratio on one demands a money sink. */
        for (j = 0; j < m->n_edges; j++) {
            const gia_edge *ed = &m->edges[j];
            int k, cu, drawn = 0, twin = 0;
            if (ed->from != i || ed->role != GIA_ROLE_USED) continue;
            if (m->nodes[ed->to].kind != GIA_NODE_SINK) {
                fprintf(stderr, "engine: module '%s': its used leg ends at '%s', "
                                "which is a %s; used energy goes to a sink\n",
                        nd->id, m->nodes[ed->to].id,
                        gia_node_kind_name(m->nodes[ed->to].kind));
                gia_model_free(m);
                return false;
            }
            cu = gia_node_carrier(m, ed->to);
            for (k = 0; k < m->n_edges; k++) {
                const gia_edge *o = &m->edges[k];
                if (o->to == i && o->role == GIA_ROLE_CONTROL &&
                    o->use_ratio > 0.0 && gia_node_carrier(m, o->from) == cu)
                    drawn = 1;
                if (k != j && o->from == i && o->role == GIA_ROLE_USED &&
                    gia_node_carrier(m, o->to) == cu)
                    twin = 1;
            }
            if (twin) {
                fprintf(stderr, "engine: module '%s' has two used legs for "
                                "carrier '%s'; it needs one\n",
                        nd->id, gia_carrier_name(m, cu));
                gia_model_free(m);
                return false;
            }
            if (!drawn) {
                fprintf(stderr, "engine: module '%s': its used leg to '%s' "
                                "carries carrier '%s', and no control of that "
                                "carrier has use_ratio > 0\n",
                        nd->id, m->nodes[ed->to].id, gia_carrier_name(m, cu));
                gia_model_free(m);
                return false;
            }
        }
        for (j = 0; j < m->n_edges; j++) {
            const gia_edge *ed = &m->edges[j];
            int k, cc, found = 0;
            if (ed->to != i || ed->role != GIA_ROLE_CONTROL ||
                ed->use_ratio <= 0.0) continue;
            cc = gia_node_carrier(m, ed->from);
            for (k = 0; k < m->n_edges; k++)
                if (m->edges[k].from == i && m->edges[k].role == GIA_ROLE_USED &&
                    gia_node_carrier(m, m->edges[k].to) == cc)
                    found = 1;
            if (!found) {
                fprintf(stderr, "engine: module '%s' draws carrier '%s' from "
                                "control '%s' (use_ratio %g) and has no used "
                                "leg to a sink of that carrier (ADR 0017)\n",
                        nd->id, gia_carrier_name(m, cc),
                        m->nodes[ed->from].id, ed->use_ratio);
                gia_model_free(m);
                return false;
            }
        }

        /* What a module passes through keeps its carrier. The energy leg is
         * what it consumes and the outputs are what it delivers, so those must
         * agree. A control may be any carrier -- a work gate whose rate is set
         * by a price or a population is a perfectly ordinary model -- because
         * what it is drawn of leaves on the used leg, checked above. */
        {
            int cin = -1;
            for (j = 0; j < m->n_edges; j++)
                if (m->edges[j].to == i && m->edges[j].role == GIA_ROLE_ENERGY)
                    cin = gia_node_carrier(m, m->edges[j].from);
            for (j = 0; j < m->n_edges; j++) {
                if (m->edges[j].from != i) continue;
                if (m->edges[j].role == GIA_ROLE_USED) continue;  /* above */
                if (gia_node_carrier(m, m->edges[j].to) == cin) continue;
                fprintf(stderr,
                        "engine: module '%s' takes carrier '%s' and delivers "
                        "carrier '%s' to '%s'; a module passes its input "
                        "through, it does not convert one carrier into "
                        "another (Odum 1972 SecXV is the transactor's job)\n",
                        nd->id, gia_carrier_name(m, cin),
                        gia_carrier_name(m, gia_node_carrier(m, m->edges[j].to)),
                        m->nodes[m->edges[j].to].id);
                gia_model_free(m);
                return false;
            }
        }
    }

    jparams = cJSON_GetObjectItemCaseSensitive(root, "simulation_params");
    m->t_end      = num_field(jparams, "t_val", 1.0);
    m->order      = (int)num_field(jparams, "derivative_order", 1);
    if (m->order < 0)              m->order = 0;
    if (m->order > GIA_MAX_ORDER)  m->order = GIA_MAX_ORDER;
    {
        const cJSON *g = cJSON_GetObjectItemCaseSensitive(jparams,
                                                          "generative_mode");
        m->generative = cJSON_IsBool(g) ? (bool)cJSON_IsTrue(g) : true;
    }
    return true;
}

void gia_model_free(gia_model *m) {
    if (!m) return;
    free(m->nodes);
    free(m->edges);
    m->nodes   = NULL;
    m->edges   = NULL;
    m->n_nodes = 0;
    m->n_edges = 0;
}

/* ================================================================== *
 * 8b. Emergy and transformity — the second accounting
 * ================================================================== */

/* The flow a module drives, at the operating point. A module's pathways carry
 * no law of their own, so asking an edge for its flow would give the pathway's
 * default rather than the module's — which is what the emergy pass would
 * otherwise carry transformity along. */
static double module_flow(const gia_model *m, int ni, const double *q) {
    const gia_node *nd = &m->nodes[ni];
    int    j, ae = -1, ctrl = -1, gi = -1;
    double g;

    for (j = 0; j < m->n_edges; j++) {
        const gia_edge *e = &m->edges[j];
        if (e->to == ni && e->role == GIA_ROLE_ENERGY)   ae   = e->from;
        if (e->to == ni && e->role == GIA_ROLE_CONTROL)  ctrl = e->from;
        if (e->to == ni && e->role == GIA_ROLE_GOODS_IN) gi   = e->from;
    }
    switch (nd->kind) {
        case GIA_NODE_INTERACTION:
            if (ae < 0) return 0.0;
            switch (nd->mod_action) {
                case GIA_ACTION_MULTIPLY:
                    g = nd->mod_k;
                    for (j = 0; j < m->n_edges; j++) {
                        const gia_edge *e = &m->edges[j];
                        if (e->to == ni && e->role == GIA_ROLE_CONTROL &&
                            e->from >= 0)
                            g *= q[e->from];
                    }
                    return g * q[ae];
                case GIA_ACTION_DIVIDE: {
                    double qc = (ctrl >= 0) ? q[ctrl] : 1.0;
                    if (qc < GIA_EPS) qc = GIA_EPS;
                    return nd->mod_k * q[ae] / qc;
                }
                case GIA_ACTION_SUBTRACT:
                    if (ctrl < 0 || q[ae] - q[ctrl] <= 0.0) return 0.0;
                    return nd->mod_k * (q[ae] - q[ctrl]);
            }
            return 0.0;
        case GIA_NODE_LOOP_LIMITED: {
            double C;
            if (ae < 0) return 0.0;
            C = (nd->mod_capacity > GIA_EPS) ? nd->mod_capacity : 1.0;
            return nd->mod_k * q[ae] * C / (C + q[ae]);
        }
        case GIA_NODE_GAIN:
            return (ctrl >= 0) ? nd->mod_k * q[ctrl] : 0.0;
        case GIA_NODE_SWITCH:
            if (ctrl < 0) return 0.0;
            return (q[ctrl] > nd->mod_threshold) ? nd->mod_k : 0.0;
        case GIA_NODE_EXCHANGE:
            return (gi >= 0) ? nd->mod_k * q[gi] : 0.0;
        default:
            return 0.0;
    }
}

double gia_edge_flow(const gia_model *m, const gia_edge *e, const double *q,
                     double t) {
    double qa;
    if (!m || !e || e->from < 0 || e->to < 0 || !q) return 0.0;

    /* A pathway touching a module carries the module's flow, not its own. A
     * control carries the share of it that it is drawn at (ADR 0017), which
     * is nothing at use_ratio 0. */
    if (m->nodes[e->to].is_module) {
        if (e->role == GIA_ROLE_CONTROL)
            return e->use_ratio * module_flow(m, e->to, q);
        if (e->role == GIA_ROLE_COUNTER_IN)
            return module_flow(m, e->to, q) /
                   ((fabs(m->nodes[e->to].mod_price) > GIA_EPS)
                        ? m->nodes[e->to].mod_price : 1.0);
        return module_flow(m, e->to, q);
    }
    if (m->nodes[e->from].is_module) {
        int    ni = e->from, j;
        double W = 0.0, F = module_flow(m, ni, q);
        if (m->nodes[ni].kind == GIA_NODE_EXCHANGE)
            return (e->role == GIA_ROLE_COUNTER_OUT)
                 ? F / ((fabs(m->nodes[ni].mod_price) > GIA_EPS)
                            ? m->nodes[ni].mod_price : 1.0)
                 : F;
        if (e->role == GIA_ROLE_USED) {
            /* Every drawn control of this sink's carrier leaves here. */
            double used = 0.0;
            for (j = 0; j < m->n_edges; j++) {
                const gia_edge *c = &m->edges[j];
                if (c->to == ni && c->role == GIA_ROLE_CONTROL && c->from >= 0 &&
                    gia_node_carrier(m, c->from) == gia_node_carrier(m, e->to))
                    used += c->use_ratio * F;
            }
            return used;
        }
        for (j = 0; j < m->n_edges; j++)
            if (m->edges[j].from == ni && m->edges[j].to >= 0 &&
                m->edges[j].role != GIA_ROLE_USED)
                W += fabs(m->edges[j].weight);
        return (W > 0.0) ? F * fabs(e->weight) / W : 0.0;
    }
    qa = q[e->from];
    /* `q` holds the solved components; a held driven component's q never moves,
     * so its instantaneous value has to come from the waveform. Without this
     * the emergy pass would carry a forced source's declared value forever. */
    if (!m->nodes[e->from].integrates &&
        m->nodes[e->from].forcing.kind != GIA_FORCE_NONE)
        qa = gia_forcing_value(&m->nodes[e->from].forcing, qa, t);

    {   /* The emergy pass carries Tr along F, so F must use the driven rate. */
        double k = gia_forcing_value(&e->forcing, e->weight, t);
        gia_edge driven = *e;
        driven.forcing.kind = GIA_FORCE_NONE;
        driven.weight = k;
        if (e->forcing.kind != GIA_FORCE_NONE)
            return gia_edge_flow(m, &driven, q, t);
    }

    switch (e->logic) {
        case GIA_LOGIC_LINEAR:      return e->weight * qa;
        case GIA_LOGIC_CONSTANT:    return e->weight;
        case GIA_LOGIC_REVERSIBLE:  return e->weight * (qa - q[e->to]);
    }
    return 0.0;
}

/* Mark the edges that close a cycle.
 *
 * Odum's fourth rule: emergy already counted on the way round a loop is not
 * counted again when the loop returns. A depth-first colouring identifies
 * exactly the edges that reach back into the current path, and those carry
 * quantity without re-injecting emergy. Without this, a feedback loop
 * manufactures emergy on every pass and the excess below would measure the
 * loop rather than the co-productions. */
static void mark_back_edges(const gia_model *m, char *colour, bool *is_back,
                            int at) {
    int i;
    colour[at] = 1;                                   /* on the current path */
    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *e = &m->edges[i];
        if (e->from != at || e->to < 0) continue;
        if (colour[e->to] == 1)      is_back[i] = true;   /* reaches back */
        else if (colour[e->to] == 0) mark_back_edges(m, colour, is_back, e->to);
    }
    colour[at] = 2;                                   /* finished */
}

/* Which co-productions lie upstream of each component.
 *
 * `anc` is n x n: anc[i*n + j] is set when component j replicates and can reach
 * component i. Two inflows to the same component whose origins share a set bit
 * are carrying emergy that came from ONE co-production, so adding them would
 * count it twice — Odum's fourth rule. Computed by relaxation, which converges
 * in at most n passes over the acyclic part; back edges are excluded because
 * they carry quantity without re-injecting emergy.
 *
 * Detected rather than declared: whether two inflows share a co-production is a
 * property of the graph, and asking a modeller to mark it would make the
 * accounting depend on their spotting it. */
static void coproduct_ancestry(const gia_model *m, const bool *is_back,
                               char *anc) {
    int i, pass, n = m->n_nodes;

    memset(anc, 0, (size_t)n * (size_t)n);
    for (i = 0; i < m->n_edges; i++) {
        int a = m->edges[i].from;
        if (a < 0 || is_back[i]) continue;
        if (m->edges[i].out_mode == GIA_OUT_REPLICATE)
            anc[(size_t)a * (size_t)n + (size_t)a] = 1;   /* a replicates */
    }
    for (pass = 0; pass < n; pass++) {
        int changed = 0;
        for (i = 0; i < m->n_edges; i++) {
            int a = m->edges[i].from, b = m->edges[i].to, j;
            if (a < 0 || b < 0 || is_back[i]) continue;
            for (j = 0; j < n; j++)
                if (anc[(size_t)a * (size_t)n + (size_t)j] &&
                    !anc[(size_t)b * (size_t)n + (size_t)j]) {
                    anc[(size_t)b * (size_t)n + (size_t)j] = 1;
                    changed = 1;
                }
        }
        if (!changed) break;
    }
}

/* Combine one component's inflows under Odum's fourth rule: maximum within a
 * set sharing a co-production ancestor, sum across sets that do not.
 *
 * `mask` is n entries per contribution. Groups are merged transitively: an
 * inflow bridging two existing groups joins them, because all three then trace
 * to one co-production. */
static double combine_inflows(int n, int count, const char *masks,
                              const double *vals) {
    double gval[GIA_MAX_INFLOWS];
    char   gmask[GIA_MAX_INFLOWS * 64];
    int    ng = 0, i, j, k, wide = (n <= 64) ? n : 64;
    double total = 0.0;

    for (i = 0; i < count && i < GIA_MAX_INFLOWS; i++) {
        const char *m_i = masks + (size_t)i * (size_t)n;
        int         hit = -1;

        for (j = 0; j < ng; j++) {
            for (k = 0; k < wide; k++)
                if (m_i[k] && gmask[(size_t)j * 64 + (size_t)k]) break;
            if (k < wide) { hit = j; break; }
        }
        if (hit < 0) {
            if (ng >= GIA_MAX_INFLOWS) { total += vals[i]; continue; }
            gval[ng] = vals[i];
            for (k = 0; k < wide; k++) gmask[(size_t)ng * 64 + (size_t)k] = m_i[k];
            ng++;
            continue;
        }
        /* Same co-production: take the larger, do not add. */
        if (vals[i] > gval[hit]) gval[hit] = vals[i];
        for (k = 0; k < wide; k++)
            if (m_i[k]) gmask[(size_t)hit * 64 + (size_t)k] = 1;

        /* This inflow may have bridged two groups; merge any that now overlap. */
        for (j = ng - 1; j >= 0; j--) {
            if (j == hit) continue;
            for (k = 0; k < wide; k++)
                if (gmask[(size_t)j * 64 + (size_t)k] &&
                    gmask[(size_t)hit * 64 + (size_t)k]) break;
            if (k == wide) continue;
            if (gval[j] > gval[hit]) gval[hit] = gval[j];
            for (k = 0; k < wide; k++)
                if (gmask[(size_t)j * 64 + (size_t)k])
                    gmask[(size_t)hit * 64 + (size_t)k] = 1;
            for (k = j; k < ng - 1; k++) {
                gval[k] = gval[k + 1];
                memcpy(gmask + (size_t)k * 64, gmask + (size_t)(k + 1) * 64, 64);
            }
            ng--;
            if (hit > j) hit--;
        }
    }
    for (j = 0; j < ng; j++) total += gval[j];
    return total;
}

bool gia_emergy_at(const gia_model *m, double t, double *em, double *tr) {
    double *q = NULL, *flow = NULL, *em_in = NULL, *out_tot = NULL;
    double *em_prev = NULL;
    bool   *is_back = NULL;
    char   *colour = NULL, *anc = NULL, *masks = NULL;
    int     i, pass, n;
    bool    ok = false;

    if (!m || m->n_nodes <= 0) return false;
    /* NFR-LIM-001: refuse, never truncate. */
    if (gia_emergy_check_limits(m, NULL) != GIA_OK) return false;
    n = m->n_nodes;

    q       = (double *)calloc((size_t)n, sizeof(double));
    em_in   = (double *)calloc((size_t)n, sizeof(double));
    em_prev = (double *)calloc((size_t)n, sizeof(double));
    out_tot = (double *)calloc((size_t)n, sizeof(double));
    colour  = (char   *)calloc((size_t)n, sizeof(char));
    anc     = (char   *)calloc((size_t)n * (size_t)n, sizeof(char));
    masks   = (char   *)calloc((size_t)GIA_MAX_INFLOWS * (size_t)n, sizeof(char));
    flow    = (double *)calloc((size_t)(m->n_edges > 0 ? m->n_edges : 1),
                               sizeof(double));
    is_back = (bool   *)calloc((size_t)(m->n_edges > 0 ? m->n_edges : 1),
                               sizeof(bool));
    if (!q || !em_in || !em_prev || !out_tot || !colour || !flow || !is_back ||
        !anc || !masks) goto done;

    if (!gia_network_state(m, t, q, NULL)) goto done;

    for (i = 0; i < m->n_edges; i++) flow[i] = gia_edge_flow(m, &m->edges[i], q, t);

    for (i = 0; i < n; i++)
        if (colour[i] == 0) mark_back_edges(m, colour, is_back, i);

    /* Total energy leaving each component, which is the denominator a partition
     * shares emergy out in proportion to. The used energy is not a product:
     * it carries no emergy onward (ADR 0017), so it takes no share. */
    for (i = 0; i < m->n_edges; i++) {
        int a = m->edges[i].from;
        if (a < 0 || is_back[i] || m->edges[i].role == GIA_ROLE_USED) continue;
        out_tot[a] += fabs(flow[i]);
    }

    coproduct_ancestry(m, is_back, anc);

    /* Propagate from the sources. The acyclic part is at most n levels deep, so
     * n relaxation passes carry emergy all the way through it. */
    for (pass = 0; pass < n + 1; pass++) {
        int b;

        for (i = 0; i < n; i++) em_in[i] = 0.0;

        /* Per component, not per edge: the fourth rule combines a component's
         * inflows against each other, so they have to be gathered first. */
        for (b = 0; b < n; b++) {
            double vals[GIA_MAX_INFLOWS];
            int    count = 0;

            memset(masks, 0, (size_t)GIA_MAX_INFLOWS * (size_t)n);

            for (i = 0; i < m->n_edges && count < GIA_MAX_INFLOWS; i++) {
                const gia_edge *e = &m->edges[i];
                int    a = e->from;
                double f, carried;

                if (a < 0 || e->to != b || is_back[i]) continue;
                if (e->role == GIA_ROLE_USED) continue;
                f = fabs(flow[i]);
                if (f <= 0.0) continue;

                if (!m->nodes[a].integrates && m->nodes[a].quality_input > 0.0) {
                    /* A boundary source injects quality at its declared rate. */
                    carried = f * m->nodes[a].quality_input;
                } else if (e->out_mode == GIA_OUT_REPLICATE) {
                    /* Co-production: each product carries the WHOLE emergy of
                     * the process, because each required all of it. */
                    carried = em_prev[a];
                } else {
                    /* Partition: a share of one kind of flow, so a share of the
                     * emergy. */
                    carried = (out_tot[a] > 0.0)
                            ? em_prev[a] * (f / out_tot[a]) : 0.0;
                }
                vals[count] = carried;
                memcpy(masks + (size_t)count * (size_t)n,
                       anc + (size_t)a * (size_t)n, (size_t)n);
                count++;
            }
            em_in[b] = combine_inflows(n, count, masks, vals);
        }
        memcpy(em_prev, em_in, (size_t)n * sizeof(double));
    }

    /* A boundary source has no inflow, so its em_in is zero -- but its empower
     * is not: it is what the source delivers. Reporting the inflow for a source
     * would show a sun with no power. Report what it emits. */
    for (i = 0; i < m->n_edges; i++) {
        int a = m->edges[i].from;
        if (a < 0 || is_back[i]) continue;
        if (!m->nodes[a].integrates && m->nodes[a].quality_input > 0.0)
            em_in[a] += fabs(flow[i]) * m->nodes[a].quality_input;
    }

    if (em) memcpy(em, em_in, (size_t)n * sizeof(double));
    if (tr) {
        for (i = 0; i < n; i++) {
            /* Transformity is emergy per unit quantity: the component's quality.
             * A source states its own rather than deriving one. */
            if (!m->nodes[i].integrates && m->nodes[i].quality_input > 0.0)
                tr[i] = m->nodes[i].quality_input;
            else
                tr[i] = (fabs(q[i]) > GIA_EPS) ? em_in[i] / q[i] : 0.0;
        }
    }
    ok = true;

done:
    free(q); free(flow); free(em_in); free(em_prev); free(out_tot);
    free(colour); free(is_back); free(anc); free(masks);
    return ok;
}

/* The emergy each node receives (em) and sends onward (out_em) at t, under
 * the same rules as gia_emergy_at. */
static bool emergy_out(const gia_model *m, double t, double *em, double *out_em) {
    double *q = NULL, *flow = NULL, *out_tot = NULL;
    bool   *is_back = NULL;
    char   *colour = NULL;
    int     i, n = m->n_nodes;
    bool    ok = false;

    q       = (double *)calloc((size_t)n, sizeof(double));
    out_tot = (double *)calloc((size_t)n, sizeof(double));
    colour  = (char   *)calloc((size_t)n, sizeof(char));
    flow    = (double *)calloc((size_t)(m->n_edges > 0 ? m->n_edges : 1),
                               sizeof(double));
    is_back = (bool   *)calloc((size_t)(m->n_edges > 0 ? m->n_edges : 1),
                               sizeof(bool));
    if (!q || !out_tot || !colour || !flow || !is_back) goto done;

    if (!gia_emergy_at(m, t, em, NULL))        goto done;
    if (!gia_network_state(m, t, q, NULL))     goto done;

    for (i = 0; i < n; i++) out_em[i] = 0.0;
    for (i = 0; i < m->n_edges; i++) flow[i] = gia_edge_flow(m, &m->edges[i], q, t);
    for (i = 0; i < n; i++)
        if (colour[i] == 0) mark_back_edges(m, colour, is_back, i);
    for (i = 0; i < m->n_edges; i++) {
        int a = m->edges[i].from;
        if (a < 0 || is_back[i] || m->edges[i].role == GIA_ROLE_USED) continue;
        out_tot[a] += fabs(flow[i]);
    }

    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *e = &m->edges[i];
        int    a = e->from;
        double f;
        if (a < 0 || e->to < 0 || is_back[i]) continue;
        if (e->role == GIA_ROLE_USED) continue;     /* carries no emergy */
        f = fabs(flow[i]);
        if (f <= 0.0) continue;
        if (!m->nodes[a].integrates && m->nodes[a].quality_input > 0.0)
            out_em[a] += f * m->nodes[a].quality_input;
        else if (e->out_mode == GIA_OUT_REPLICATE)
            out_em[a] += em[a];
        else
            out_em[a] += (out_tot[a] > 0.0) ? em[a] * (f / out_tot[a]) : 0.0;
    }
    ok = true;
done:
    free(q); free(out_tot); free(colour); free(flow); free(is_back);
    return ok;
}

double gia_emergy_excess(const gia_model *m, double t) {
    double *em = NULL, *out_em = NULL, excess = 0.0;
    int     i, n;

    if (!m || m->n_nodes <= 0) return 0.0;
    n = m->n_nodes;
    em     = (double *)calloc((size_t)n, sizeof(double));
    out_em = (double *)calloc((size_t)n, sizeof(double));
    if (em && out_em && emergy_out(m, t, em, out_em)) {
        /* Emergy created, component by component. A partition contributes
         * nothing: what leaves equals what arrived. A replication contributes
         * the whole inflow again for every product past the first. */
        for (i = 0; i < n; i++) {
            double made = out_em[i] - em[i];
            if (!m->nodes[i].integrates) continue;   /* a boundary source is not creating */
            if (made > 0.0) excess += made;
        }
    }
    free(em); free(out_em);
    return excess;
}

static gia_status em_fail(gia_status st, const char *reason, const char **why) {
    if (why) *why = reason;
    return st;
}

/* NFR-LIM-001. */
gia_status gia_emergy_check_limits(const gia_model *m, const char **why) {
    int b, i, replicates = 0;
    if (!m || m->n_nodes <= 0)
        return em_fail(GIA_E_ARG, "gia_emergy_check_limits: NULL or empty model", why);
    for (i = 0; i < m->n_edges; i++)
        if (m->edges[i].from >= 0 && m->edges[i].out_mode == GIA_OUT_REPLICATE) replicates = 1;
    if (replicates && m->n_nodes > GIA_MAX_INFLOWS)
        return em_fail(GIA_E_LIMIT,
                       "the emergy pass tracks co-production ancestry across at most "
                       "GIA_MAX_INFLOWS (64) nodes; this model has more and a co-production "
                       "(NFR-LIM-001)", why);
    for (b = 0; b < m->n_nodes; b++) {
        int count = 0;
        for (i = 0; i < m->n_edges; i++)
            if (m->edges[i].from >= 0 && m->edges[i].to == b &&
                m->edges[i].role != GIA_ROLE_USED) count++;
        if (count > GIA_MAX_INFLOWS)
            return em_fail(GIA_E_LIMIT,
                           "a node has more than GIA_MAX_INFLOWS (64) inflows, which the "
                           "emergy pass cannot combine under rule 4 (NFR-LIM-001)", why);
    }
    return GIA_OK;
}

/* FR-EM-002, FR-EM-004. */
gia_status gia_emergy_source_term(const gia_model *m, double t, int node,
                                  double *phi, const char **why) {
    double    *em, *out_em;
    gia_status st;
    if (!m || !phi || m->n_nodes <= 0)
        return em_fail(GIA_E_ARG, "gia_emergy_source_term: NULL or empty model", why);
    if (node < 0 || node >= m->n_nodes)
        return em_fail(GIA_E_ARG, "gia_emergy_source_term: node out of range", why);
    if (!m->nodes[node].is_module &&
        (m->nodes[node].kind == GIA_NODE_SOURCE || m->nodes[node].kind == GIA_NODE_SINK ||
         m->nodes[node].kind == GIA_NODE_CONSTANT))
        return em_fail(GIA_E_ARG, "gia_emergy_source_term: a source, sink or constant is "
                       "habitat, not a process ([02 Eq 3.6])", why);
    if ((st = gia_emergy_check_limits(m, why)) != GIA_OK) return st;
    em     = (double *)calloc((size_t)m->n_nodes, sizeof(double));
    out_em = (double *)calloc((size_t)m->n_nodes, sizeof(double));
    if (!em || !out_em || !emergy_out(m, t, em, out_em)) {
        free(em); free(out_em);
        return em_fail(GIA_E_DOMAIN, "gia_emergy_source_term: the network does not solve "
                       "at t", why);
    }
    *phi = out_em[node] - em[node];
    free(em); free(out_em);
    return GIA_OK;
}

/* FR-EM-005 — [22 Eq 6-8]. */
static gia_oform oform(gia_oform_kind kind, int rows, int cols) {
    gia_oform f;
    memset(&f, 0, sizeof(f));
    f.kind = kind; f.rows = rows; f.cols = cols;
    return f;
}

gia_oform gia_oform_scalar(double a) {
    gia_oform f = oform(GIA_OF_SCALAR, 1, 1);
    f.v[0][0] = a;
    return f;
}

gia_oform gia_oform_binary(double em_u) {
    gia_oform f = oform(GIA_OF_BINARY, 2, 1);
    f.v[0][0] = em_u; f.v[1][0] = em_u;
    return f;
}

gia_oform gia_oform_duet(double em_u1, double em_u2) {
    gia_oform f = oform(GIA_OF_DUET, 1, 2);
    f.v[0][0] = em_u1; f.v[0][1] = em_u2;
    return f;
}

gia_oform gia_oform_duet_binary(double a1, double a2) {
    gia_oform f = oform(GIA_OF_DUET_BINARY, 2, 2);
    f.v[0][0] = a1; f.v[0][1] = a2;
    f.v[1][0] = a2; f.v[1][1] = a1;
    return f;
}

/* FR-EM-006 — [06b Eq 2], [02 Eq 14.11.5]; PLAN R12. */
gia_status gia_circle_product(const gia_oform *a, const gia_oform *b, gia_circle *out,
                              const char **why) {
    int i, j;
    if (!a || !b || !out)
        return em_fail(GIA_E_ARG, "gia_circle_product: NULL argument", why);
    if (a->cols != 1 || b->rows != 1 || a->rows < 1 || a->rows > 2 || b->cols < 1 || b->cols > 2)
        return em_fail(GIA_E_ARG, "gia_circle_product: the left factor must be a column "
                       "(a1; a2) and the right a row [b1, b2] ([06b Eq 2])", why);
    memset(out, 0, sizeof(*out));
    out->rows = a->rows; out->cols = b->cols;
    for (i = 0; i < a->rows; i++)
        for (j = 0; j < b->cols; j++) {
            out->pair[i][j][0] = a->v[i][0];
            out->pair[i][j][1] = b->v[0][j];
        }
    return GIA_OK;
}

/* The cardinal reduction: each pair to its product. */
gia_oform gia_circle_reduce(const gia_circle *c) {
    gia_oform f;
    int       i, j;
    memset(&f, 0, sizeof(f));
    if (!c) return f;
    f.kind = (c->rows == 1 && c->cols == 1) ? GIA_OF_SCALAR : GIA_OF_DUET_BINARY;
    f.rows = c->rows; f.cols = c->cols;
    for (i = 0; i < c->rows; i++)
        for (j = 0; j < c->cols; j++)
            f.v[i][j] = c->pair[i][j][0] * c->pair[i][j][1];
    return f;
}

static gia_status balance_sums(const gia_balance_term *in, int n_in,
                               const gia_balance_term *out, int n_out,
                               double *sin, double *sout, const char **why) {
    int i;
    if ((n_in > 0 && !in) || (n_out > 0 && !out) || n_in < 0 || n_out < 0)
        return em_fail(GIA_E_ARG, "gia_emergy balance: NULL terms or negative count", why);
    *sin = *sout = 0.0;
    for (i = 0; i < n_in; i++)  *sin  += in[i].weight * in[i].value;
    for (i = 0; i < n_out; i++) *sout += out[i].weight * out[i].value;
    if (!isfinite(*sin) || !isfinite(*sout))
        return em_fail(GIA_E_RANGE, "gia_emergy balance: a sum is not finite", why);
    return GIA_OK;
}

/* FR-EM-007 — [02 Eq 3.21]. */
gia_status gia_emergy_global_balance(const gia_balance_term *in, int n_in,
                                     const gia_balance_term *out, int n_out,
                                     double *residual, const char **why) {
    double     a, b;
    gia_status st;
    if (!residual) return em_fail(GIA_E_ARG, "gia_emergy_global_balance: residual is NULL", why);
    if ((st = balance_sums(in, n_in, out, n_out, &a, &b, why)) != GIA_OK) return st;
    *residual = a - b;
    return GIA_OK;
}

gia_status gia_emergy_balance_solve(const gia_balance_term *in, int n_in,
                                    const gia_balance_term *out, int n_out,
                                    const double *phi_w, int n_phi, double *phi,
                                    const char **why) {
    double     a, b, w = 0.0, x;
    int        k;
    gia_status st;
    if (!phi_w || !phi || n_phi < 1)
        return em_fail(GIA_E_ARG, "gia_emergy_balance_solve: no source terms to solve for", why);
    if ((st = balance_sums(in, n_in, out, n_out, &a, &b, why)) != GIA_OK) return st;
    for (k = 0; k < n_phi; k++) w += phi_w[k];
    if (w == 0.0 || !isfinite(w))
        return em_fail(GIA_E_DOMAIN, "gia_emergy_balance_solve: the source-term weights sum "
                       "to zero, so the unknown is not determined", why);
    x = (b - a) / w;
    for (k = 0; k < n_phi; k++) phi[k] = phi_w[k] * x;
    return GIA_OK;
}

/* ================================================================== *
 * 6. Ordinality
 * ================================================================== */

/* Depth-first reachability from `at`, looking for `target`. */
/* ADR 0014. A closed pathway is closed by legs that carry quantity, and a
 * module is passed through but is not a component.
 *
 * A control at use_ratio 0 carries nothing, so it is never walked. Since
 * module-hosted laws made controls into edges, walking them let a leg carrying
 * nothing close a loop: a pure accumulator was reported at maximum ordinality
 * and the MOP step declined to generate. A drawn control does carry quantity
 * (ADR 0017), but only through its module to the used leg and so to a sink --
 * never into the product -- so it closes no pathway through the product. */

/* Whether quantity arriving at module `mod` on a leg of role `in` leaves along
 * a leg of role `out`. Quantity passes through a module and does not change
 * kind on the way: a transactor hands goods on as goods and counter-flow on as
 * counter-flow, and never turns one into the other. Every other module passes
 * its energy input to its products, and a drawn control to its used leg. */
static bool module_passes(const gia_node *mod, gia_role in, gia_role out) {
    if (mod->kind == GIA_NODE_EXCHANGE)
        return (in == GIA_ROLE_GOODS_IN   && out == GIA_ROLE_GOODS_OUT) ||
               (in == GIA_ROLE_COUNTER_IN && out == GIA_ROLE_COUNTER_OUT);
    if (out == GIA_ROLE_USED) return in == GIA_ROLE_CONTROL;
    return in == GIA_ROLE_ENERGY;
}

/* A module can be entered on independent streams -- a transactor by goods and
 * by counter-flow, any other module by its energy and by a drawn control -- so
 * it gets one visited slot per stream. A component needs only one. */
#define GIA_STREAMS 3

static int visit_slot(const gia_model *m, int node, gia_role via) {
    int stream = 0;
    if (m->nodes[node].is_module) {
        if (via == GIA_ROLE_COUNTER_IN) stream = 1;
        if (via == GIA_ROLE_CONTROL)    stream = 2;
    }
    return node * GIA_STREAMS + stream;
}

static bool reaches(const gia_model *m, int at, gia_role via, int target,
                    bool *seen) {
    const gia_node *here = &m->nodes[at];
    int i;
    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *ed = &m->edges[i];
        int slot;
        if (ed->from != at || ed->to < 0) continue;
        /* Read, not carried. A drawn control is walked, on its own stream,
         * so reaching a module by it never blocks an energy pathway arriving
         * through the same module later. */
        if (ed->role == GIA_ROLE_CONTROL && ed->use_ratio <= 0.0) continue;
        if (here->is_module && !module_passes(here, via, ed->role)) continue;
        if (ed->to == target) return true;
        slot = visit_slot(m, ed->to, ed->role);
        if (!seen[slot]) {
            seen[slot] = true;
            if (reaches(m, ed->to, ed->role, target, seen)) return true;
        }
    }
    return false;
}

/* Components, not nodes: a module holds nothing, so it is not one. */
static int component_count(const gia_model *m) {
    int i, n = 0;
    for (i = 0; i < m->n_nodes; i++) if (!m->nodes[i].is_module) n++;
    return n;
}

/* The one cycle scan. gia_mark_cycles, gia_ordinality and gia_generate all read
 * it, so the ordinality that is reported and the one that decides emergence
 * cannot disagree. A module's flag is always false. Returns how many components
 * are on a closed pathway. */
static int cycles_into(const gia_model *m, bool *out) {
    int   i, count = 0;
    bool *seen;

    if (!m || m->n_nodes <= 0) return 0;
    seen = (bool *)malloc((size_t)m->n_nodes * GIA_STREAMS * sizeof(bool));
    if (!seen) return 0;

    for (i = 0; i < m->n_nodes; i++) {
        int k;
        out[i] = false;
        if (m->nodes[i].is_module) continue;
        for (k = 0; k < m->n_nodes * GIA_STREAMS; k++) seen[k] = false;
        seen[visit_slot(m, i, GIA_ROLE_NONE)] = true;
        out[i] = reaches(m, i, GIA_ROLE_NONE, i, seen);
        if (out[i]) count++;
    }
    free(seen);
    return count;
}

int gia_mark_cycles(gia_model *m) {
    bool *flags;
    int   i, count;

    if (!m || m->n_nodes <= 0) return 0;
    flags = (bool *)calloc((size_t)m->n_nodes, sizeof(bool));
    if (!flags) return 0;
    count = cycles_into(m, flags);
    for (i = 0; i < m->n_nodes; i++) m->nodes[i].on_cycle = flags[i];
    free(flags);
    return count;
}

int gia_component_count(const gia_model *m) {
    return m ? component_count(m) : 0;
}

static bool is_component(const gia_node *nd);

/* FR-ORD-004: the proxy. Habitat sits outside the system [10], so only
 * components are counted (ADR 0021). */
double gia_closure(gia_model *m) {
    int i, on = 0, n = 0;
    if (!m || m->n_nodes <= 0) return 0.0;
    (void)gia_mark_cycles(m);
    for (i = 0; i < m->n_nodes; i++) {
        if (!is_component(&m->nodes[i])) continue;
        n++;
        if (m->nodes[i].on_cycle) on++;
    }
    return n > 0 ? (double)on / (double)n : 0.0;
}

double gia_ordinality(gia_model *m) {
    return gia_closure(m);
}

/* Does component `from` reach `to` along quantity-carrying legs? The ADR 0014
 * walk of reaches(), from a fresh visited set. */
static bool reaches_from(const gia_model *m, int from, int to, bool *seen) {
    int k;
    for (k = 0; k < m->n_nodes * GIA_STREAMS; k++) seen[k] = false;
    seen[visit_slot(m, from, GIA_ROLE_NONE)] = true;
    return reaches(m, from, GIA_ROLE_NONE, to, seen);
}

/* A leg that carries quantity: everything but a read control. */
static bool carries(const gia_edge *e) {
    return !(e->role == GIA_ROLE_CONTROL && e->use_ratio <= 0.0);
}

/* 2: both feed one interaction module along a quantity-carrying leg. */
static bool feed_one_interaction(const gia_model *m, int a, int b) {
    int i, j;
    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *ea = &m->edges[i];
        if (ea->from != a || ea->to < 0 || !carries(ea)) continue;
        if (!m->nodes[ea->to].is_module || m->nodes[ea->to].kind != GIA_NODE_INTERACTION) continue;
        for (j = 0; j < m->n_edges; j++) {
            const gia_edge *eb = &m->edges[j];
            if (eb->from == b && eb->to == ea->to && carries(eb)) return true;
        }
    }
    return false;
}

/* 1/2: both are products, replicated, of one process. */
static bool coproducts(const gia_model *m, int a, int b) {
    int i, j;
    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *ea = &m->edges[i];
        if (ea->to != a || ea->from < 0 || ea->out_mode != GIA_OUT_REPLICATE) continue;
        for (j = 0; j < m->n_edges; j++) {
            const gia_edge *eb = &m->edges[j];
            if (eb->to == b && eb->from == ea->from && eb->out_mode == GIA_OUT_REPLICATE)
                return true;
        }
    }
    return false;
}

gia_status gia_ordinality_record(const gia_model *m, gia_ordinality_rec *r, const char **why) {
    gia_ordinality_rec out = {0, 0, 0, 0, 0};
    bool              *seen;
    int                i, j;

    if (!m || !r || m->n_nodes < 0) {
        if (why) *why = "gia_ordinality_record: NULL argument";
        return GIA_E_ARG;
    }
    seen = (bool *)malloc((size_t)(m->n_nodes > 0 ? m->n_nodes : 1) * GIA_STREAMS * sizeof(bool));
    if (!seen) {
        if (why) *why = "gia_ordinality_record: out of memory";
        return GIA_E_NOMEM;
    }
    for (i = 0; i < m->n_nodes; i++) {
        if (!is_component(&m->nodes[i])) continue;
        out.k++;
        for (j = i + 1; j < m->n_nodes; j++) {
            if (!is_component(&m->nodes[j])) continue;
            if (reaches_from(m, i, j, seen) && reaches_from(m, j, i, seen)) out.n22++;
            else if (feed_one_interaction(m, i, j))                       out.n2++;
            else if (coproducts(m, i, j))                                  out.nhalf++;
            else                                                           out.nunrelated++;
        }
    }
    free(seen);
    *r = out;
    return GIA_OK;
}

/* FR-ORD-003: every couple 2/2. */
bool gia_at_maximum_ordinality(gia_model *m) {
    gia_ordinality_rec r;
    if (!m || gia_ordinality_record(m, &r, NULL) != GIA_OK || r.k < 2) return false;
    return r.n22 == r.k * (r.k - 1) / 2;
}

/* ================================================================== *
 * 6b. Drift on network trajectories (FR-IDC-011, FR-IDC-014)
 * ================================================================== */

static gia_status drift_fail(gia_status st, const char *reason, const char **why) {
    if (why) *why = reason;
    return st;
}

/* FR-IDC-011. A constant flow matrix makes the solution exp(A t) Q(0), a
 * superposition of exponentials with affine exponents, so the incipient and
 * traditional solutions coincide term by term (PLAN R16) and the drift is
 * zero by [06 §4 (i)] -- stated, not computed. */
gia_status gia_solution_drift(const gia_model *m, int order, double *out,
                              const char **why) {
    int i;
    if (!m || !out || m->n_nodes <= 0)
        return drift_fail(GIA_E_ARG, "gia_solution_drift: NULL or empty model", why);
    if (order < 1 || order > 4)
        return drift_fail(GIA_E_ARG, "gia_solution_drift: order must be 1..4 (FR-IDC-011)", why);
    if (!gia_flow_matrix_is_constant(m))
        return drift_fail(GIA_E_UNSUPPORTED,
                          "solution drift is defined only for a known equation [06 §4]; the "
                          "sources define none for a network whose flow matrix is not "
                          "constant (FR-IDC-011, PLAN §6)", why);
    for (i = 0; i < m->n_nodes; i++) out[i] = 0.0;
    return GIA_OK;
}

/* The flow balance F(Q, t) = A(Q, t) x(Q, t), and, when `second` is not NULL
 * and A is constant, A^2 x as well -- exactly Q'' (numerics.md N7). */
static bool flow_balance(const gia_model *m, const double *q, double t,
                         double *first, double *second) {
    gia_matrix A;
    double    *x, *ax;
    int        i, j, n = m->n_nodes, dim = n + forcing_extra_count(m) + 1;
    bool       ok = false;

    memset(&A, 0, sizeof(A));
    x  = (double *)malloc((size_t)dim * sizeof(double));
    ax = (double *)malloc((size_t)dim * sizeof(double));
    if (!x || !ax || !gia_build_flow_matrix(m, q, t, &A)) goto done;
    augment_state(m, q, t, x, dim);
    for (i = 0; i < dim; i++) {
        double acc = 0.0;
        for (j = 0; j < dim; j++) acc += gia_matrix_at(&A, i, j) * x[j];
        ax[i] = acc;
    }
    for (i = 0; i < n; i++) first[i] = ax[i];
    if (second) {
        for (i = 0; i < n; i++) {
            double acc = 0.0;
            for (j = 0; j < dim; j++) acc += gia_matrix_at(&A, i, j) * ax[j];
            second[i] = acc;
        }
    }
    ok = true;
done:
    gia_matrix_free(&A);
    free(x); free(ax);
    return ok;
}

/* A component in the sense of srs.md §1: not a module, not a boundary node. */
static bool is_component(const gia_node *nd) {
    return !nd->is_module && nd->kind != GIA_NODE_SOURCE &&
           nd->kind != GIA_NODE_SINK && nd->kind != GIA_NODE_CONSTANT;
}

bool gia_node_is_component(const gia_model *m, int i) {
    return m && i >= 0 && i < m->n_nodes && is_component(&m->nodes[i]);
}

/* F at time s along the engine's own solution (numerics.md N7: Q at t +- h
 * comes from the engine, not a re-integration). */
static bool flow_at(const gia_model *m, double s, double *q, double *f) {
    return gia_network_state(m, s, q, NULL) && flow_balance(m, q, s, f, NULL);
}

/* Richardson-extrapolated derivative of F at t with base step h:
 * central where t - 2h >= 0 is not needed, one-sided at t near 0. */
static bool flow_slope(const gia_model *m, double t, double h, double *d,
                       double *q, double *f0, double *f1, double *f2) {
    int    i, n = m->n_nodes, pass;
    double hh;
    bool   central = t - h >= 0.0;
    double *est[2];

    est[0] = d;
    est[1] = (double *)malloc((size_t)n * sizeof(double));
    if (!est[1]) return false;
    for (pass = 0; pass < 2; pass++) {
        hh = pass == 0 ? h : h / 2.0;
        if (central) {
            if (!flow_at(m, t + hh, q, f1) || !flow_at(m, t - hh, q, f2)) { free(est[1]); return false; }
            for (i = 0; i < n; i++) est[pass][i] = (f1[i] - f2[i]) / (2.0 * hh);
        } else {
            if (!flow_at(m, t, q, f0) || !flow_at(m, t + hh, q, f1) ||
                !flow_at(m, t + 2.0 * hh, q, f2)) { free(est[1]); return false; }
            for (i = 0; i < n; i++)
                est[pass][i] = (-3.0 * f0[i] + 4.0 * f1[i] - f2[i]) / (2.0 * hh);
        }
    }
    /* Both forms are second order: (4 D(h/2) - D(h)) / 3. */
    for (i = 0; i < n; i++) d[i] = (4.0 * est[1][i] - est[0][i]) / 3.0;
    free(est[1]);
    return true;
}

gia_status gia_drift_projection(const gia_model *m, double t, double dt,
                                double *out, const char **why) {
    double    *q = NULL, *f = NULL, *qq = NULL, *d1 = NULL, *d2 = NULL, *w0 = NULL,
              *w1 = NULL, *w2 = NULL;
    int        i, n;
    gia_status st = GIA_OK;

    if (!m || !out || m->n_nodes <= 0)
        return drift_fail(GIA_E_ARG, "gia_drift_projection: NULL or empty model", why);
    if (!(t >= 0.0) || !isfinite(t) || !isfinite(dt))
        return drift_fail(GIA_E_ARG, "gia_drift_projection: t must be >= 0 and finite", why);
    n  = m->n_nodes;
    q  = (double *)malloc((size_t)n * sizeof(double));
    f  = (double *)malloc((size_t)n * sizeof(double));
    qq = (double *)malloc((size_t)n * sizeof(double));
    d1 = (double *)malloc((size_t)n * sizeof(double));
    d2 = (double *)malloc((size_t)n * sizeof(double));
    w0 = (double *)malloc((size_t)n * sizeof(double));
    w1 = (double *)malloc((size_t)n * sizeof(double));
    w2 = (double *)malloc((size_t)n * sizeof(double));
    if (!q || !f || !qq || !d1 || !d2 || !w0 || !w1 || !w2) {
        st = drift_fail(GIA_E_NOMEM, "gia_drift_projection: out of memory", why);
        goto done;
    }
    if (!gia_network_state(m, t, q, NULL)) {
        st = drift_fail(GIA_E_DOMAIN, "gia_drift_projection: the network does not solve at t", why);
        goto done;
    }
    if (gia_flow_matrix_is_constant(m)) {
        /* Q'' = A^2 x exactly (N7). */
        if (!flow_balance(m, q, t, f, d1)) {
            st = drift_fail(GIA_E_NOMEM, "gia_drift_projection: flow matrix", why);
            goto done;
        }
    } else {
        double h0 = 1e-3 * (fabs(t) > 1.0 ? fabs(t) : 1.0);
        if (!flow_balance(m, q, t, f, NULL) ||
            !flow_slope(m, t, h0, d1, qq, w0, w1, w2) ||
            !flow_slope(m, t, h0 / 2.0, d2, qq, w0, w1, w2)) {
            st = drift_fail(GIA_E_DOMAIN, "gia_drift_projection: the network does not solve "
                            "near t", why);
            goto done;
        }
        for (i = 0; i < n; i++)
            if (fabs(d1[i] - d2[i]) > 1e-6 * (fabs(d2[i]) + fabs(f[i]) + 1e-300)) {
                st = drift_fail(GIA_E_CONVERGENCE,
                                "gia_drift_projection: Q'' did not converge under "
                                "extrapolation (numerics N7); a crossing or switch near t?",
                                why);
                goto done;
            }
        for (i = 0; i < n; i++) d1[i] = d2[i];
    }
    for (i = 0; i < n; i++) {
        const gia_node *nd = &m->nodes[i];
        double v = 0.0;
        if (is_component(nd) && q[i] != 0.0)
            v = (d1[i] - f[i] * f[i] / q[i]) * dt * dt / 2.0;
        if (!isfinite(v)) {
            st = drift_fail(GIA_E_RANGE, "gia_drift_projection: drift not finite (NFR-NUM-003)", why);
            goto done;
        }
        qq[i] = v;
    }
    for (i = 0; i < n; i++) out[i] = qq[i];
done:
    free(q); free(f); free(qq); free(d1); free(d2); free(w0); free(w1); free(w2);
    return st;
}

/* ================================================================== *
 * 7. Mode 1 — functional trajectories
 * ================================================================== */

bool gia_sample_at(const gia_model *m, double t,
                   double *q, double *idc, double *tdc, double *psi) {
    int i;
    if (!m || m->n_nodes <= 0) return false;

    if (q) {
        if (!gia_network_state(m, t, q, psi)) {
            for (i = 0; i < m->n_nodes; i++) q[i] = m->nodes[i].q0;
            if (psi) *psi = 0.0;
            return false;
        }
    } else if (psi) {
        *psi = 0.0;
    }
    for (i = 0; i < m->n_nodes; i++) {
        const gia_phi *p = &m->nodes[i].phi;
        if (idc) idc[i] = gia_idc_derivative(p, m->order, t);
        if (tdc) tdc[i] = gia_tdc_derivative(p, m->order, t);
    }
    return true;
}

bool gia_write_trajectories(const gia_model *m, const char *path, int steps) {
    FILE   *f;
    int     i, s;
    double  dt;
    double *q, *em, *tr, *dp;

    if (!m || !path) return false;
    if (steps < 1) steps = 1;

    f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "engine: cannot write %s: ", path);
        perror(NULL);
        return false;
    }

    /* IF-OUT-001. Per node: the network quantity Q from the matrix
     * exponential (classical), its empower and transformity, and the
     * output-projection drift of [09 Eq 13] over one output step
     * (FR-IDC-014). The _idc/_tdc/_drift and psi_network columns are gone:
     * they applied the drift identity to a per-node phi the engine invented,
     * decoupled from Q (PLAN E4). */
    fprintf(f, "time");
    for (i = 0; i < m->n_nodes; i++)
        fprintf(f, ",%s_Q,%s_Em,%s_Tr,%s_drift_proj",
                m->nodes[i].id, m->nodes[i].id, m->nodes[i].id, m->nodes[i].id);
    fprintf(f, ",conservation,emergy_excess\n");

    q  = (double *)calloc((size_t)m->n_nodes, sizeof(double));
    em = (double *)calloc((size_t)m->n_nodes, sizeof(double));
    tr = (double *)calloc((size_t)m->n_nodes, sizeof(double));
    dp = (double *)calloc((size_t)m->n_nodes, sizeof(double));
    if (!q || !em || !tr || !dp) {
        free(q); free(em); free(tr); free(dp);
        fclose(f); return false;
    }

    dt = m->t_end / (double)steps;
    for (s = 0; s <= steps; s++) {
        double t = (double)s * dt;
        bool   have_dp, have_em;

        (void)gia_sample_at(m, t, q, NULL, NULL, NULL);
        /* A model past the emergy pass's limits is refused, not truncated
         * (NFR-LIM-001): its emergy cells are left empty. */
        have_em = gia_emergy_at(m, t, em, tr);
        have_dp = gia_drift_projection(m, t, dt, dp, NULL) == GIA_OK;

        fprintf(f, "%.6f", t);
        for (i = 0; i < m->n_nodes; i++) {
            fprintf(f, ",%.10g,", q[i]);
            if (have_em) fprintf(f, "%.10g,%.10g,", em[i], tr[i]); else fprintf(f, ",,");
            if (have_dp) fprintf(f, "%.10g", dp[i]);
        }
        fprintf(f, ",%.10g,", gia_conservation_residual(m, t));
        if (have_em) fprintf(f, "%.10g", gia_emergy_excess(m, t));
        fprintf(f, "\n");
    }

    free(q); free(em); free(tr); free(dp);
    fclose(f);
    return true;
}

/* ================================================================== *
 * 8. Mode 2 — the generative ordinal step
 * ================================================================== */

/* Replace in place where the key already exists, so the field keeps its
 * position in the object. Delete-then-add would move it to the end, which
 * costs nothing semantically but turns a one-line change into two hunks in
 * any text diff of the seed against the output. */
static void set_string(cJSON *obj, const char *key, const char *value) {
    cJSON *item = cJSON_CreateString(value);
    if (!item) return;
    if (cJSON_GetObjectItemCaseSensitive(obj, key)) {
        if (!cJSON_ReplaceItemInObjectCaseSensitive(obj, key, item))
            cJSON_Delete(item);
    } else {
        cJSON_AddItemToObject(obj, key, item);
    }
}

/* Sort node indices by id, so candidates are enumerated in id order and the
 * output is a function of the model rather than of how its file was ordered.
 *
 * qsort's comparator takes no context, so the id travels with the index in
 * each element rather than through a file-scope pointer to the model: that
 * pointer was shared state, and two models generating on two threads could
 * sort by each other's ids (NFR-REE-001; numerics.md §3). */
typedef struct { const char *id; int idx; } id_key;

static int by_id(const void *x, const void *y) {
    return strcmp(((const id_key *)x)->id, ((const id_key *)y)->id);
}

static bool sort_by_id(const gia_model *m, int *ids, int n) {
    id_key *keys;
    int     i;
    if (n < 2) return true;
    keys = (id_key *)malloc((size_t)n * sizeof(id_key));
    if (!keys) return false;
    for (i = 0; i < n; i++) { keys[i].id = m->nodes[ids[i]].id; keys[i].idx = ids[i]; }
    qsort(keys, (size_t)n, sizeof(id_key), by_id);
    for (i = 0; i < n; i++) ids[i] = keys[i].idx;
    free(keys);
    return true;
}

/* The component graph's condensation, enough of it to choose candidates:
 * comps[0..n) are the component indices in id order; scc[a] labels the
 * strongly connected component of comps[a] by its lowest position; src[s] and
 * snk[s] say whether SCC s has no in-edge, no out-edge, from another SCC. */
typedef struct {
    int   n, n_scc, n_src, n_snk;
    int  *comps, *scc;
    bool *src, *snk;
} gia_condensation;

static void condensation_free(gia_condensation *c) {
    free(c->comps); free(c->scc); free(c->src); free(c->snk);
    memset(c, 0, sizeof(*c));
}

static bool condensation(const gia_model *m, gia_condensation *c) {
    bool *reach = NULL, *seen = NULL;
    int   a, b, n = 0;

    memset(c, 0, sizeof(*c));
    c->comps = (int *)malloc((size_t)m->n_nodes * sizeof(int));
    if (!c->comps) return false;
    for (a = 0; a < m->n_nodes; a++) if (is_component(&m->nodes[a])) c->comps[n++] = a;
    c->n   = n;
    c->scc = (int *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    c->src = (bool *)calloc((size_t)(n > 0 ? n : 1), sizeof(bool));
    c->snk = (bool *)calloc((size_t)(n > 0 ? n : 1), sizeof(bool));
    reach  = (bool *)calloc((size_t)(n > 0 ? n * n : 1), sizeof(bool));
    seen   = (bool *)malloc((size_t)m->n_nodes * GIA_STREAMS * sizeof(bool));
    if (!c->scc || !c->src || !c->snk || !reach || !seen || !sort_by_id(m, c->comps, n)) {
        free(reach); free(seen); condensation_free(c);
        return false;
    }
    for (a = 0; a < n; a++)
        for (b = 0; b < n; b++)
            reach[a * n + b] = a == b || reaches_from(m, c->comps[a], c->comps[b], seen);
    for (a = 0; a < n; a++) {
        c->scc[a] = a;
        for (b = 0; b < a; b++)
            if (reach[a * n + b] && reach[b * n + a]) { c->scc[a] = c->scc[b]; break; }
        if (c->scc[a] == a) { c->n_scc++; c->src[a] = c->snk[a] = true; }
    }
    for (a = 0; a < n; a++)
        for (b = 0; b < n; b++)
            if (c->scc[a] != c->scc[b] && reach[a * n + b]) {
                c->snk[c->scc[a]] = false;
                c->src[c->scc[b]] = false;
            }
    for (a = 0; a < n; a++) {
        if (c->scc[a] != a) continue;
        if (c->src[a]) c->n_src++;
        if (c->snk[a]) c->n_snk++;
    }
    free(reach); free(seen);
    return true;
}

/* A `linear` pathway, as the seed's own pathways are written. */
static bool add_pathway(cJSON *root, const char *from, const char *to, double weight) {
    cJSON *edges = cJSON_GetObjectItemCaseSensitive(root, "edges");
    cJSON *ne;
    if (!cJSON_IsArray(edges) && !(edges = cJSON_AddArrayToObject(root, "edges"))) return false;
    ne = cJSON_CreateObject();
    if (!ne) return false;
    cJSON_AddStringToObject(ne, "source", from);
    cJSON_AddStringToObject(ne, "target", to);
    cJSON_AddStringToObject(ne, "logic", "linear");
    cJSON_AddNumberToObject(ne, "weight", weight);
    cJSON_AddStringToObject(ne, "flow_type", "max_empower_pathway");
    cJSON_AddItemToArray(edges, ne);
    return true;
}

/* Total empower at t_end: the components' emergy per unit time, summed.
 * Habitat sits outside the system (ADR 0021), so a sink's accumulation is not
 * counted. False when the emergy pass is unavailable (NFR-LIM-001). */
static bool total_empower(const gia_model *m, double *total) {
    double *em = (double *)calloc((size_t)m->n_nodes, sizeof(double));
    double  s  = 0.0;
    int     i;
    bool    okay;
    if (!em) return false;
    okay = gia_emergy_at(m, m->t_end, em, NULL);
    for (i = 0; okay && i < m->n_nodes; i++) if (is_component(&m->nodes[i])) s += em[i];
    free(em);
    if (!okay || !isfinite(s)) return false;
    *total = s;
    return true;
}

/* The total empower of `root` with the pathway from -> to added. */
static bool trial_empower(const cJSON *root, const char *from, const char *to, double weight,
                          double *total) {
    cJSON    *trial = cJSON_Duplicate(root, 1);
    gia_model tm;
    bool      okay = false;
    if (!trial) return false;
    if (add_pathway(trial, from, to, weight) && gia_model_load(&tm, trial)) {
        okay = total_empower(&tm, total);
        gia_model_free(&tm);
    }
    cJSON_Delete(trial);
    return okay;
}

/* ADR 0021 §2 (FR-ORD-005). Below Maximum Ordinality, repeat: form the
 * condensation of the component graph; the candidates are the pathways from a
 * component in a sink SCC to one in a source SCC (never within one SCC); add
 * the `linear` one, of the seed's mean edge weight, that maximises total
 * empower at t_end -- the discrete form of [02 Eq 5.3] -- ties within 1e-9
 * relative going to the lexicographically first (from, to). Each addition
 * lowers #sources + #sinks of the condensation, so the step ends within that
 * many additions; a model that does not is returned unchanged, as is one the
 * emergy pass cannot evaluate. No component is ever added (ADR 0015's E is
 * retired). */
cJSON *gia_generate(const gia_model *m) {
    cJSON           *out;
    gia_condensation c;
    gia_model        cur;
    double           mean_w = 0.0;
    int              i, added = 0, bound = -1;

    if (!m || !m->root) return NULL;

    /* The output always starts as a faithful copy of the seed. If nothing
     * below fires, it compares equal to the input and the validator will
     * correctly report a functional run. */
    out = cJSON_Duplicate(m->root, 1);
    if (!out) return NULL;
    if (!m->generative) return out;

    for (i = 0; i < m->n_edges; i++) mean_w += m->edges[i].weight;
    mean_w = (m->n_edges > 0) ? mean_w / (double)m->n_edges : 1.0;

    for (;;) {
        double best = 0.0;
        int    a, b, best_a = -1, best_b = -1, n_cand = 0;
        char   from[256], to[256];

        if (!gia_model_load(&cur, out)) break;
        if (!condensation(&cur, &c)) { gia_model_free(&cur); break; }
        if (c.n < 2 || c.n_scc == 1) {
            if (c.n < 2)
                printf("  MOP generative step: fewer than two components, so no couple "
                       "to relate.\n");
            else
                printf("  MOP generative step: at Maximum Ordinality (strongly connected) "
                       "after %d addition%s.\n", added, added == 1 ? "" : "s");
            condensation_free(&c); gia_model_free(&cur);
            break;
        }
        if (bound < 0) bound = c.n_src + c.n_snk;
        if (added >= bound) {
            printf("  MOP generative step: declined. %d additions did not reach Maximum "
                   "Ordinality, so the seed is returned unchanged.\n", added);
            condensation_free(&c); gia_model_free(&cur);
            cJSON_Delete(out);
            return cJSON_Duplicate(m->root, 1);
        }
        for (a = 0; a < c.n; a++) {
            if (!c.snk[c.scc[a]]) continue;
            for (b = 0; b < c.n; b++) {
                double tot;
                if (c.scc[a] == c.scc[b] || !c.src[c.scc[b]]) continue;
                n_cand++;
                if (!trial_empower(out, cur.nodes[c.comps[a]].id, cur.nodes[c.comps[b]].id,
                                   mean_w, &tot))
                    continue;
                if (best_a < 0 || tot > best + 1e-9 * (fabs(best) > 1.0 ? fabs(best) : 1.0)) {
                    best = tot; best_a = a; best_b = b;
                }
            }
        }
        if (best_a >= 0) {
            snprintf(from, sizeof from, "%s", cur.nodes[c.comps[best_a]].id);
            snprintf(to, sizeof to, "%s", cur.nodes[c.comps[best_b]].id);
        }
        condensation_free(&c); gia_model_free(&cur);
        if (best_a < 0 || !add_pathway(out, from, to, mean_w)) {
            printf("  MOP generative step: declined. The emergy pass could not evaluate "
                   "the candidates, so the seed is returned unchanged.\n");
            cJSON_Delete(out);
            return cJSON_Duplicate(m->root, 1);
        }
        added++;
        printf("  MOP generative step %d: %s -> %s, the maximum total empower %.6g of %d "
               "candidate%s ([02 Eq 5.3]).\n", added, from, to, best, n_cand,
               n_cand == 1 ? "" : "s");
    }
    if (added > 0)
        set_string(out, "system_name", "Evolved Self-Organizing Graph (Post-MOP Ordinal Step)");
    return out;
}

/* Same numbers as gia_write_trajectories(), laid out for a terminal.
 *
 * One section per component rather than one wide table: four columns per
 * component times N components does not fit a terminal.
 *
 * Q and the phi columns are different quantities and are labelled so, because
 * confusing them is easy and costly. Q is the simulation: it comes from
 * Q(t) = exp(A t) Q(0) and depends on the whole graph. The phi columns are the
 * single-component analytic form of Sections 1-3, which depends only on that
 * component's own parameters and ignores every edge -- they are what makes the
 * drift theorem checkable against a closed form, and they are not a trajectory
 * of the system.
 *
 * A source shows this most sharply: Odum 1972 SecII holds a source at its value
 * rather than integrating it, so its Q is flat, while its phi form grows.
 */
void gia_print_trajectories(const gia_model *m, int steps) {
    int     i, s;
    double  dt;
    double *q, *dp;

    if (!m || m->n_nodes <= 0) return;
    if (steps < 1) steps = 1;
    dt = m->t_end / (double)steps;

    q  = (double *)malloc((size_t)m->n_nodes * sizeof(double));
    dp = (double *)malloc((size_t)m->n_nodes * sizeof(double));
    if (!q || !dp) { free(q); free(dp); return; }

    for (i = 0; i < m->n_nodes; i++) {
        const gia_node *nd = &m->nodes[i];

        printf("\n  %s  (%s, %s)\n", nd->id, gia_node_kind_name(nd->kind),
               nd->integrates ? "integrated" : "held, not integrated");
        printf("  %12s %16s | %16s\n", "time", "Q (network)", "drift_proj");
        printf("  %12s %16s | %16s\n", "------------", "----------------",
               "----------------");

        for (s = 0; s <= steps; s++) {
            double t = (double)s * dt;
            (void)gia_sample_at(m, t, q, NULL, NULL, NULL);
            if (gia_drift_projection(m, t, dt, dp, NULL) == GIA_OK)
                printf("  %12.4f %16.6g | %16.6g\n", t, q[i], dp[i]);
            else
                printf("  %12.4f %16.6g | %16s\n", t, q[i], "(no estimate)");
        }
    }

    printf("\n  Q is the simulation: Q(t) = exp(A t) Q(0), where A is the flow\n");
    printf("  matrix assembled from the pathway laws -- classical numerics.\n");
    printf("  drift_proj is the output-projection drift of [09 Eq 13] over one\n");
    printf("  output step: (Q'' - Q'^2/Q) dt^2/2 along that same trajectory.\n");

    printf("\n  network summary at horizon t = %g\n", m->t_end);
    printf("  %-20s %10s %16s %16s\n",
           "component", "solved", "Q(0)", "Q(t_end)");
    printf("  %-20s %10s %16s %16s\n",
           "--------------------", "------", "----------------",
           "----------------");

    if (gia_network_state(m, m->t_end, q, NULL)) {
        for (i = 0; i < m->n_nodes; i++)
            printf("  %-20s %10s %16.6g %16.6g\n",
                   m->nodes[i].id, m->nodes[i].integrates ? "yes" : "held",
                   m->nodes[i].q0, q[i]);
    }

    {
        bool   constA = gia_flow_matrix_is_constant(m);
        printf("\n  flow matrix         %s\n",
               constA ? "constant -- solution drift is exactly zero [06 §4 (i)]"
                      : "state-dependent -- no solution drift is defined (FR-IDC-011)");
        double *em = (double *)calloc((size_t)m->n_nodes, sizeof(double));
        double *tr = (double *)calloc((size_t)m->n_nodes, sizeof(double));
        if (em && tr && gia_emergy_at(m, m->t_end, em, tr)) {
            double xs = gia_emergy_excess(m, m->t_end);
            printf("\n  emergy at t = %g   (empower, and transformity Em/Q)\n",
                   m->t_end);
            printf("  %-20s %18s %18s\n", "component", "Em", "Tr");
            printf("  %-20s %18s %18s\n", "--------------------",
                   "------------------", "------------------");
            for (i = 0; i < m->n_nodes; i++)
                printf("  %-20s %18.6g %18.6g\n", m->nodes[i].id, em[i], tr[i]);
            printf("\n  emergy created       %.6g\n", xs);
            if (xs > 0.0) {
                printf("    Emergy is NOT conserved here, and that is correct:\n");
                printf("    a co-production gives each product the whole emergy\n");
                printf("    of the process, because each required all of it.\n");
                printf("    This excess is irreducible to the inputs -- the\n");
                printf("    thing a conservative calculus cannot express, and\n");
                printf("    Giannantoni's reason for needing another one.\n");
            } else {
                printf("    Every bifurcation here is a partition, so emergy is\n");
                printf("    conserved. Mark an edge output_mode \"replicate\" to\n");
                printf("    make it a co-production and this becomes non-zero.\n");
            }
        }
        free(em); free(tr);

        if (gia_system_is_closed(m)) {
            printf("  system              closed\n");
            printf("  conservation        %.3e  (must be ~0)\n",
                   gia_conservation_residual(m, m->t_end));
        } else {
            printf("  system              open -- a pathway leaves a held\n");
            printf("                      component, which delivers quantity\n");
            printf("                      without being depleted (Odum SecII)\n");
            printf("  net boundary inflow %.6g  (expected, not an error)\n",
                   gia_conservation_residual(m, m->t_end));
        }
    }

    free(q); free(dp);
}

