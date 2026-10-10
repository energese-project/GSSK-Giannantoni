/* mop_seed.c — the seed's `mop` block (IF-JSON-001) and the MOP CSV (IF-OUT-002).
 *
 * Strict parsing (AGENTS.md, NFR-ROB-001): every object's keys are checked
 * against the grammar, a key given twice is an error, and every number must be
 * finite. A load error is GIA_E_ARG with a static `why` and the JSON path in
 * the caller's `detail`. No file-scope state, no printing (NFR-REE-001,
 * NFR-ERR-001).
 */

#include "mop_seed.h"

#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The load-error state threaded through the parser. */
typedef struct {
    char        *detail;
    size_t       cap;
    const char **why;
} perr;

static gia_status load_error(const perr *e, const char *reason, const char *fmt, ...) {
    if (e->detail && e->cap > 0) {
        va_list ap;
        size_t  n;
        va_start(ap, fmt);
        vsnprintf(e->detail, e->cap, fmt, ap);
        va_end(ap);
        n = strlen(e->detail);
        if (n + 2 < e->cap) snprintf(e->detail + n, e->cap - n, ": %s", reason);
    }
    if (e->why) *e->why = reason;
    return GIA_E_ARG;
}

/* Every key of obj is one of `allowed` (NULL-terminated) and appears once. */
static gia_status keys_ok(const perr *e, const cJSON *obj, const char *const *allowed,
                          const char *path) {
    const cJSON *c, *d;
    for (c = obj->child; c; c = c->next) {
        const char *const *a;
        int                known = 0;
        if (!c->string) return load_error(e, "a member without a name", "%s", path);
        for (a = allowed; *a; a++) if (strcmp(*a, c->string) == 0) known = 1;
        if (!known)
            return load_error(e, "unknown key (IF-JSON-001)", "%s.%.64s", path, c->string);
        for (d = obj->child; d != c; d = d->next)
            if (strcmp(d->string, c->string) == 0)
                return load_error(e, "key given twice", "%s.%.64s", path, c->string);
    }
    return GIA_OK;
}

static int finite_number(const cJSON *j) {
    return cJSON_IsNumber(j) && isfinite(j->valuedouble);
}

static gia_status get_real(const perr *e, const cJSON *obj, const char *key, const char *path,
                           double *out) {
    const cJSON *j = cJSON_GetObjectItemCaseSensitive(obj, key);
    if (!j) return load_error(e, "missing", "%s.%s", path, key);
    if (!finite_number(j)) return load_error(e, "must be a finite number", "%s.%s", path, key);
    *out = j->valuedouble;
    return GIA_OK;
}

/* A number, or [re, im]. */
static gia_status get_complex(const perr *e, const cJSON *obj, const char *key,
                              const char *path, double complex *out) {
    const cJSON *j = cJSON_GetObjectItemCaseSensitive(obj, key);
    if (!j) return load_error(e, "missing", "%s.%s", path, key);
    if (finite_number(j)) { *out = j->valuedouble; return GIA_OK; }
    if (cJSON_IsArray(j) && cJSON_GetArraySize(j) == 2 && finite_number(j->child) &&
        finite_number(j->child->next)) {
        *out = j->child->valuedouble + j->child->next->valuedouble * I;
        return GIA_OK;
    }
    return load_error(e, "must be a finite number or [re, im]", "%s.%s", path, key);
}

/* An array of exactly n finite numbers. */
static gia_status get_reals(const perr *e, const cJSON *obj, const char *key, const char *path,
                            int n, double *out) {
    const cJSON *j = cJSON_GetObjectItemCaseSensitive(obj, key), *c;
    int          i = 0;
    if (!j) return load_error(e, "missing", "%s.%s", path, key);
    if (!cJSON_IsArray(j) || cJSON_GetArraySize(j) != n)
        return load_error(e, "wrong length or not an array", "%s.%s", path, key);
    for (c = j->child; c; c = c->next, i++) {
        if (!finite_number(c)) return load_error(e, "must be a finite number", "%s.%s[%d]", path, key, i);
        out[i] = c->valuedouble;
    }
    return GIA_OK;
}

static int as_positive_int(const cJSON *j, int *out) {
    if (!finite_number(j) || j->valuedouble != floor(j->valuedouble) ||
        j->valuedouble < 1.0 || j->valuedouble > (double)INT_MAX)
        return 0;
    *out = (int)j->valuedouble;
    return 1;
}

static int gcd(int a, int b) {
    while (b) { int r = a % b; a = b; b = r; }
    return a;
}

static gia_status parse_k(const perr *e, const cJSON *j, gia_rational *k) {
    static const char *const keys[] = {"num", "den", NULL};
    gia_status st;
    int        num, den, g;
    if (!j) return load_error(e, "missing", "mop.k");
    if (cJSON_IsObject(j)) {
        if ((st = keys_ok(e, j, keys, "mop.k")) != GIA_OK) return st;
        if (!as_positive_int(cJSON_GetObjectItemCaseSensitive(j, "num"), &num) ||
            !as_positive_int(cJSON_GetObjectItemCaseSensitive(j, "den"), &den))
            return load_error(e, "num and den must be integers >= 1 (FR-MOP-002)", "mop.k");
    } else if (as_positive_int(j, &num)) {
        den = 1;
    } else {
        return load_error(e, "must be an integer >= 1 or {\"num\", \"den\"} (FR-MOP-002)", "mop.k");
    }
    g = gcd(num, den);
    k->num = num / g; k->den = den / g;
    return GIA_OK;
}

/* The index of the component named by j (a string). */
static gia_status component(const perr *e, const gia_model *m, const cJSON *j, const char *path,
                            int *idx) {
    int i;
    if (!cJSON_IsString(j)) return load_error(e, "must be a component id (a string)", "%s", path);
    for (i = 0; i < m->n_nodes; i++)
        if (m->nodes[i].id && strcmp(m->nodes[i].id, j->valuestring) == 0) {
            if (!gia_node_is_component(m, i))
                return load_error(e, "names a module or habitat, not a component (ADR 0021)",
                                  "%s", path);
            *idx = i;
            return GIA_OK;
        }
    return load_error(e, "no such node", "%s", path);
}

static gia_status parse_samples(const perr *e, const cJSON *j, const char *path,
                                gia_mop_couple_spec *c) {
    const cJSON *row;
    int          n, i = 0;
    if (!cJSON_IsArray(j) || (n = cJSON_GetArraySize(j)) < 2)
        return load_error(e, "must be an array of at least two [t, re, im]", "%s.samples", path);
    c->t = (double *)malloc((size_t)n * sizeof(double));
    c->v = (double complex *)malloc((size_t)n * sizeof(double complex));
    if (!c->t || !c->v) {
        if (e->why) *e->why = "mop: out of memory";
        return GIA_E_NOMEM;
    }
    for (row = j->child; row; row = row->next, i++) {
        const cJSON *t, *re, *im;
        if (!cJSON_IsArray(row) || cJSON_GetArraySize(row) != 3 ||
            !finite_number(t = row->child) || !finite_number(re = t->next) ||
            !finite_number(im = re->next))
            return load_error(e, "must be [t, re, im], finite numbers", "%s.samples[%d]", path, i);
        if (i == 0 ? t->valuedouble != 0.0 : !(t->valuedouble > c->t[i - 1]))
            return load_error(e, "t must start at 0 and strictly increase", "%s.samples[%d]",
                              path, i);
        c->t[i] = t->valuedouble;
        c->v[i] = re->valuedouble + im->valuedouble * I;
    }
    c->beta.kind = GIA_BETA_SAMPLES;
    c->beta.n    = n;
    c->beta.t    = c->t;
    c->beta.v    = c->v;
    return GIA_OK;
}

static gia_status parse_couple(const perr *e, const gia_model *m, const cJSON *j, int i,
                               gia_mop_couple_spec *c) {
    static const char *const keys[] = {"from", "to", "a", "b", "p", "samples", NULL};
    char         path[48];
    const cJSON *samples;
    gia_status   st;
    int          affine;

    snprintf(path, sizeof path, "mop.beta[%d]", i);
    if (!cJSON_IsObject(j)) return load_error(e, "must be an object", "%s", path);
    if ((st = keys_ok(e, j, keys, path)) != GIA_OK) return st;
    if ((st = component(e, m, cJSON_GetObjectItemCaseSensitive(j, "from"), path, &c->from)) != GIA_OK ||
        (st = component(e, m, cJSON_GetObjectItemCaseSensitive(j, "to"), path, &c->to)) != GIA_OK)
        return st;
    if (c->from == c->to) return load_error(e, "from and to are the same component", "%s", path);
    samples = cJSON_GetObjectItemCaseSensitive(j, "samples");
    affine  = cJSON_GetObjectItemCaseSensitive(j, "a") || cJSON_GetObjectItemCaseSensitive(j, "b") ||
              cJSON_GetObjectItemCaseSensitive(j, "p");
    if (samples && affine)
        return load_error(e, "give either a, b, p or samples, not both", "%s", path);
    if (samples) return parse_samples(e, samples, path, c);
    if (!affine) return load_error(e, "needs a, b, p or samples", "%s", path);
    c->beta.kind = GIA_BETA_AFFINE;
    if ((st = get_complex(e, j, "a", path, &c->beta.a)) != GIA_OK ||
        (st = get_complex(e, j, "b", path, &c->beta.b)) != GIA_OK)
        return st;
    return get_real(e, j, "p", path, &c->beta.p);
}

static int couple_cmp(const gia_model *m, const gia_mop_couple_spec *x, const gia_mop_couple_spec *y) {
    int c = strcmp(m->nodes[x->from].id, m->nodes[y->from].id);
    return c ? c : strcmp(m->nodes[x->to].id, m->nodes[y->to].id);
}

static gia_status parse_beta(const perr *e, const gia_model *m, const cJSON *j, gia_mop_seed *s) {
    const cJSON *c;
    int          n, i;
    gia_status   st;

    if (!j) return load_error(e, "missing", "mop.beta");
    if (cJSON_IsString(j)) {
        if (strcmp(j->valuestring, "network") != 0)
            return load_error(e, "the only string form is \"network\" (FR-MOP-008)", "mop.beta");
        s->form = GIA_MOP_BETA_NETWORK;
        return GIA_OK;
    }
    if (!cJSON_IsArray(j)) return load_error(e, "must be an array of couples or \"network\"", "mop.beta");
    s->form = GIA_MOP_BETA_COUPLES;
    n = cJSON_GetArraySize(j);
    if (n == 0) return GIA_OK;
    s->couples = (gia_mop_couple_spec *)calloc((size_t)n, sizeof(gia_mop_couple_spec));
    if (!s->couples) {
        if (e->why) *e->why = "mop: out of memory";
        return GIA_E_NOMEM;
    }
    for (c = j->child, i = 0; c; c = c->next, i++) {
        s->n_couples = i + 1;       /* so that a failure frees what this couple holds */
        if ((st = parse_couple(e, m, c, i, &s->couples[i])) != GIA_OK) return st;
    }
    /* Insertion sort by (from id, to id): n <= N(N-1), and seeds are usually in order. */
    for (i = 1; i < n; i++) {
        gia_mop_couple_spec x = s->couples[i];
        int                 k = i - 1;
        while (k >= 0 && couple_cmp(m, &s->couples[k], &x) > 0) { s->couples[k + 1] = s->couples[k]; k--; }
        s->couples[k + 1] = x;
    }
    /* SAMPLES point at their own arrays, which moved with them. */
    for (i = 0; i < n; i++)
        if (s->couples[i].beta.kind == GIA_BETA_SAMPLES) {
            s->couples[i].beta.t = s->couples[i].t;
            s->couples[i].beta.v = s->couples[i].v;
        }
    for (i = 1; i < n; i++)
        if (couple_cmp(m, &s->couples[i - 1], &s->couples[i]) == 0)
            return load_error(e, "couple given twice", "mop.beta (%.32s -> %.32s)",
                              m->nodes[s->couples[i].from].id, m->nodes[s->couples[i].to].id);
    return GIA_OK;
}

static gia_status parse_reference(const perr *e, const gia_model *m, const cJSON *j, gia_mop_seed *s) {
    gia_status st;
    int        i;
    if (!j) {   /* the first two components by id */
        s->ref[0] = s->ref[1] = -1;
        for (i = 0; i < m->n_nodes; i++) {
            if (!gia_node_is_component(m, i)) continue;
            if (s->ref[0] < 0 || strcmp(m->nodes[i].id, m->nodes[s->ref[0]].id) < 0) {
                s->ref[1] = s->ref[0]; s->ref[0] = i;
            } else if (s->ref[1] < 0 || strcmp(m->nodes[i].id, m->nodes[s->ref[1]].id) < 0) {
                s->ref[1] = i;
            }
        }
        if (s->ref[1] < 0) return load_error(e, "the MOP needs at least two components", "mop");
        return GIA_OK;
    }
    if (!cJSON_IsArray(j) || cJSON_GetArraySize(j) != 2)
        return load_error(e, "must be two component ids", "mop.reference");
    if ((st = component(e, m, j->child, "mop.reference[0]", &s->ref[0])) != GIA_OK ||
        (st = component(e, m, j->child->next, "mop.reference[1]", &s->ref[1])) != GIA_OK)
        return st;
    if (s->ref[0] == s->ref[1]) return load_error(e, "a couple needs two components", "mop.reference");
    return GIA_OK;
}

static gia_status parse_second(const perr *e, const cJSON *j, gia_mop_seed *s) {
    static const char *const keys[] = {"alpha12_0", "c1", "c2", NULL};
    gia_status st;
    if (!cJSON_IsObject(j)) return load_error(e, "must be an object", "mop.second_equation");
    if ((st = keys_ok(e, j, keys, "mop.second_equation")) != GIA_OK ||
        (st = get_complex(e, j, "alpha12_0", "mop.second_equation", &s->alpha12_0)) != GIA_OK ||
        (st = get_real(e, j, "c1", "mop.second_equation", &s->c1)) != GIA_OK ||
        (st = get_real(e, j, "c2", "mop.second_equation", &s->c2)) != GIA_OK)
        return st;
    s->has_second = 1;
    return GIA_OK;
}

static gia_status parse_eqs(const perr *e, const gia_model *m, const cJSON *j, gia_mop_seed *s) {
    static const char *const keys[] = {"psi1", "psi2", "epsilon", "A", NULL};
    gia_status st;
    int        i;
    if (!cJSON_IsObject(j)) return load_error(e, "must be an object", "mop.eqs");
    if ((st = keys_ok(e, j, keys, "mop.eqs")) != GIA_OK ||
        (st = get_reals(e, j, "psi1", "mop.eqs", 3, s->eqs.psi1)) != GIA_OK ||
        (st = get_real(e, j, "psi2", "mop.eqs", &s->eqs.psi2)) != GIA_OK ||
        (st = get_reals(e, j, "epsilon", "mop.eqs", 3, s->eqs.eps)) != GIA_OK ||
        (st = get_real(e, j, "A", "mop.eqs", &s->eqs.A)) != GIA_OK)
        return st;
    s->eqs.N = 0;
    for (i = 0; i < m->n_nodes; i++) if (gia_node_is_component(m, i)) s->eqs.N++;
    s->has_eqs = 1;
    return GIA_OK;
}

gia_status gia_mop_seed_load(const gia_model *m, gia_mop_seed *out, char *detail, size_t cap,
                             const char **why) {
    static const char *const keys[] = {"k", "reference", "beta", "second_equation", "eqs", NULL};
    perr         e;
    gia_mop_seed s;
    const cJSON *j, *x;
    gia_status   st;

    e.detail = detail; e.cap = cap; e.why = why;
    if (detail && cap > 0) detail[0] = '\0';
    if (!m || !out) {
        if (why) *why = "gia_mop_seed_load: NULL argument";
        return GIA_E_ARG;
    }
    memset(out, 0, sizeof(*out));
    memset(&s, 0, sizeof s);
    j = m->root ? cJSON_GetObjectItemCaseSensitive(m->root, "mop") : NULL;
    if (!j) return GIA_OK;
    if (!cJSON_IsObject(j)) return load_error(&e, "must be an object", "mop");
    /* A second "mop" member at the top level would be silently shadowed. */
    for (x = j->next; x; x = x->next)
        if (x->string && strcmp(x->string, "mop") == 0)
            return load_error(&e, "key given twice", "mop");

    if ((st = keys_ok(&e, j, keys, "mop")) == GIA_OK &&
        (st = parse_k(&e, cJSON_GetObjectItemCaseSensitive(j, "k"), &s.k)) == GIA_OK &&
        (st = parse_reference(&e, m, cJSON_GetObjectItemCaseSensitive(j, "reference"), &s)) == GIA_OK &&
        (st = parse_beta(&e, m, cJSON_GetObjectItemCaseSensitive(j, "beta"), &s)) == GIA_OK &&
        (!(x = cJSON_GetObjectItemCaseSensitive(j, "second_equation")) ||
         (st = parse_second(&e, x, &s)) == GIA_OK) &&
        (!(x = cJSON_GetObjectItemCaseSensitive(j, "eqs")) || (st = parse_eqs(&e, m, x, &s)) == GIA_OK)) {
        s.present = 1;
        *out = s;
        return GIA_OK;
    }
    gia_mop_seed_free(&s);
    return st;
}

void gia_mop_seed_free(gia_mop_seed *s) {
    int i;
    if (!s) return;
    for (i = 0; i < s->n_couples && s->couples; i++) {
        free(s->couples[i].t);
        free(s->couples[i].v);
    }
    free(s->couples);
    memset(s, 0, sizeof(*s));
}

/* E_ij(t), the emergy the direct pathways i -> j carry, into E (n*n) and rel. */
static gia_status network_E(const gia_model *m, double t, double *E, unsigned char *rel,
                            const char **why) {
    const int  n = m->n_nodes;
    double    *carried = (double *)malloc((size_t)(m->n_edges > 0 ? m->n_edges : 1) * sizeof(double));
    gia_status st;
    int        e;
    if (!carried) {
        if (why) *why = "gia_mop_network: out of memory";
        return GIA_E_NOMEM;
    }
    memset(E, 0, (size_t)n * (size_t)n * sizeof(double));
    memset(rel, 0, (size_t)n * (size_t)n);
    if ((st = gia_emergy_carried(m, t, carried, why)) == GIA_OK)
        for (e = 0; e < m->n_edges; e++) {
            const gia_edge *ed = &m->edges[e];
            if (ed->from < 0 || ed->to < 0 || ed->from == ed->to) continue;
            if (!gia_node_is_component(m, ed->from) || !gia_node_is_component(m, ed->to)) continue;
            if (ed->role == GIA_ROLE_USED || (ed->role == GIA_ROLE_CONTROL && ed->use_ratio <= 0.0))
                continue;
            rel[ed->from * n + ed->to] = 1;
            E[ed->from * n + ed->to] += carried[e];
        }
    free(carried);
    return st;
}

static gia_status no_log(const char **why) {
    if (why) *why = "gia_mop_network: a pathway carries no emergy near t, so e^alpha = 0 has no "
                    "logarithm (FR-MOP-008, PLAN R8)";
    return GIA_E_RANGE;
}

gia_status gia_mop_network(const gia_model *m, double t, gia_matrioska *alpha,
                           gia_matrioska *beta, const char **why) {
    double        *E[3] = {NULL, NULL, NULL}, h, ts[3];
    unsigned char *rel = NULL;
    double complex *aa = NULL, *bb = NULL;
    unsigned char *ra = NULL, *rb = NULL;
    gia_status     st = GIA_OK;
    int            n, i, s, ns, central;

    if (!m || m->n_nodes <= 0 || !(t >= 0.0) || !isfinite(t)) {
        if (why) *why = "gia_mop_network: no model, or t < 0";
        return GIA_E_ARG;
    }
    n = m->n_nodes;
    /* beta = d ln E / dt: central where t allows, else second-order forward. */
    h       = 1e-5 * (t > 1.0 ? t : 1.0);
    central = t >= h;
    ns      = beta ? 3 : 1;
    ts[0]   = t;
    ts[1]   = central ? t - h : t + h;
    ts[2]   = central ? t + h : t + 2.0 * h;
    rel = (unsigned char *)malloc((size_t)n * (size_t)n);
    for (s = 0; s < ns; s++) E[s] = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    if (alpha) {
        aa = (double complex *)calloc((size_t)n * (size_t)n, sizeof(double complex));
        ra = (unsigned char *)calloc((size_t)n * (size_t)n, 1);
    }
    if (beta) {
        bb = (double complex *)calloc((size_t)n * (size_t)n, sizeof(double complex));
        rb = (unsigned char *)calloc((size_t)n * (size_t)n, 1);
    }
    if (!rel || !E[0] || (beta && (!E[1] || !E[2] || !bb || !rb)) || (alpha && (!aa || !ra))) {
        if (why) *why = "gia_mop_network: out of memory";
        st = GIA_E_NOMEM;
        goto done;
    }
    for (s = 0; s < ns && st == GIA_OK; s++) st = network_E(m, ts[s], E[s], rel, why);
    if (st != GIA_OK) goto done;
    for (i = 0; i < n * n; i++) {
        double l0, l1, l2;
        if (!rel[i]) continue;
        if (!(E[0][i] > 0.0)) { st = no_log(why); goto done; }
        l0 = log(E[0][i]);
        if (alpha) { aa[i] = l0; ra[i] = 1; }
        if (beta) {
            if (!(E[1][i] > 0.0) || !(E[2][i] > 0.0)) { st = no_log(why); goto done; }
            l1 = log(E[1][i]);
            l2 = log(E[2][i]);
            bb[i] = central ? (l2 - l1) / (2.0 * h) : (-3.0 * l0 + 4.0 * l1 - l2) / (2.0 * h);
            rb[i] = 1;
        }
    }
    if (alpha) { alpha->N = n; alpha->a = aa; alpha->related = ra; aa = NULL; ra = NULL; }
    if (beta)  { beta->N = n;  beta->a = bb;  beta->related = rb;  bb = NULL; rb = NULL; }
done:
    for (s = 0; s < 3; s++) free(E[s]);
    free(rel); free(aa); free(ra); free(bb); free(rb);
    return st;
}

gia_status gia_mop_reference_row(const gia_model *m, const int *ref, const gia_matrioska *full,
                                 gia_matrioska *row, const char **why) {
    int  *order, k = 0, i, j, r0, r1;
    const int n = m ? m->n_nodes : 0;

    if (!m || !full || !row || full->N != n || !full->a || !full->related) {
        if (why) *why = "gia_mop_reference_row: NULL argument or a Matrioska of the wrong size";
        return GIA_E_ARG;
    }
    order = (int *)malloc((size_t)n * sizeof(int));
    if (!order) {
        if (why) *why = "gia_mop_reference_row: out of memory";
        return GIA_E_NOMEM;
    }
    for (i = 0; i < n; i++) if (gia_node_is_component(m, i)) order[k++] = i;
    for (i = 1; i < k; i++) {            /* components by id */
        int x = order[i], p = i - 1;
        while (p >= 0 && strcmp(m->nodes[order[p]].id, m->nodes[x].id) > 0) { order[p + 1] = order[p]; p--; }
        order[p + 1] = x;
    }
    r0 = ref ? ref[0] : (k > 0 ? order[0] : -1);
    r1 = ref ? ref[1] : (k > 1 ? order[1] : -1);
    if (k < 2 || r0 == r1 || !gia_node_is_component(m, r0) || !gia_node_is_component(m, r1)) {
        free(order);
        if (why) *why = "gia_mop_reference_row: the reference couple needs two distinct components";
        return GIA_E_ARG;
    }
    /* ref[0], ref[1], then the rest in id order. */
    for (i = 0, j = 2; i < k; i++) {
        const int x = order[i];
        if (x != r0 && x != r1) order[j++] = x;
    }
    order[0] = r0; order[1] = r1;
    row->N       = k;
    row->a       = (double complex *)calloc((size_t)k * (size_t)k, sizeof(double complex));
    row->related = (unsigned char *)calloc((size_t)k * (size_t)k, 1);
    if (!row->a || !row->related) {
        free(order); gia_matrioska_free(row);
        if (why) *why = "gia_mop_reference_row: out of memory";
        return GIA_E_NOMEM;
    }
    for (j = 1; j < k; j++) {
        const int c = r0 * n + order[j];
        if (full->related[c]) { row->a[j] = full->a[c]; row->related[j] = 1; }
    }
    free(order);
    return GIA_OK;
}

/* The R_H cell of one CSV row: the residual of the reference row, or empty. */
static void write_rh(FILE *f, const gia_model *m, const int *ref, const gia_matrioska *full) {
    gia_matrioska row;
    double        R;
    int           good = 0;
    memset(&row, 0, sizeof row);
    if (gia_mop_reference_row(m, ref, full, &row, NULL) == GIA_OK)
        good = gia_harmony_residual(&row, &R, NULL) == GIA_OK;
    gia_matrioska_free(&row);
    if (good) fprintf(f, ",%.17g", R); else fprintf(f, ",");
}

/* IF-OUT-002 for "beta": "network": alpha = ln E for each related couple. */
static gia_status write_network(const gia_model *m, const gia_mop_seed *s, const char *path,
                                int steps, const char **why) {
    gia_matrioska al;
    FILE         *f;
    int          *pairs, np = 0, i, j, r, n = m->n_nodes;
    gia_status    st;
    double        dt;

    if (s->k.num != 1 || s->k.den != 1) {
        if (why) *why = "beta from the network (FR-MOP-008) is defined with k = 1 only "
                        "([23 Eq 5.5.2], PLAN R8)";
        return GIA_E_UNSUPPORTED;
    }
    if ((st = gia_mop_network(m, m->t_end, &al, NULL, why)) != GIA_OK) return st;
    pairs = (int *)malloc((size_t)n * (size_t)n * sizeof(int));
    if (!pairs) {
        gia_matrioska_free(&al);
        if (why) *why = "gia_mop_write_csv: out of memory";
        return GIA_E_NOMEM;
    }
    for (i = 0; i < n * n; i++) {      /* insertion by (from id, to id) */
        int k;
        if (!al.related[i]) continue;
        for (k = np; k > 0; k--) {
            const int p = pairs[k - 1];
            int c = strcmp(m->nodes[p / n].id, m->nodes[i / n].id);
            if (!c) c = strcmp(m->nodes[p % n].id, m->nodes[i % n].id);
            if (c <= 0) break;
            pairs[k] = p;
        }
        pairs[k] = i;
        np++;
    }
    gia_matrioska_free(&al);
    f = fopen(path, "w");
    if (!f) {
        free(pairs);
        if (why) *why = "gia_mop_write_csv: cannot open the output file";
        return GIA_E_ARG;
    }
    fprintf(f, "time");
    for (j = 0; j < np; j++) {
        const char *a = m->nodes[pairs[j] / n].id, *b = m->nodes[pairs[j] % n].id;
        fprintf(f, ",%s__%s_re,%s__%s_im", a, b, a, b);
    }
    fprintf(f, ",R_H\n");
    dt = m->t_end / (double)steps;
    for (r = 0; r <= steps; r++) {
        const double t = (double)r * dt;
        if ((st = gia_mop_network(m, t, &al, NULL, why)) != GIA_OK) {
            fclose(f); remove(path); free(pairs);
            return st;
        }
        fprintf(f, "%.6f", t);
        for (j = 0; j < np; j++)
            fprintf(f, ",%.17g,%.17g", creal(al.a[pairs[j]]), cimag(al.a[pairs[j]]));
        write_rh(f, m, s->ref, &al);
        fprintf(f, "\n");
        gia_matrioska_free(&al);
    }
    free(pairs);
    if (fclose(f) != 0) {
        remove(path);
        if (why) *why = "gia_mop_write_csv: the output file could not be written";
        return GIA_E_ARG;
    }
    return GIA_OK;
}

gia_status gia_mop_write_csv(const gia_model *m, const gia_mop_seed *s, const char *path,
                             int steps, const char **why) {
    FILE          *f;
    double         dt;
    double complex al;
    gia_matrioska  full;
    gia_status     st;
    int            i, r;

    if (!m || !s || !s->present || !path) {
        if (why) *why = "gia_mop_write_csv: no model, no mop block, or no path";
        return GIA_E_ARG;
    }
    if (steps < 1) steps = 1;
    if (s->form == GIA_MOP_BETA_NETWORK) return write_network(m, s, path, steps, why);
    /* Refuse before writing: a couple's domain is decided on [0, t_end]. */
    for (i = 0; i < s->n_couples; i++)
        if ((st = gia_mop_couple(&s->couples[i].beta, s->k, m->t_end, &al, why)) != GIA_OK)
            return st;

    f = fopen(path, "w");
    if (!f) {
        if (why) *why = "gia_mop_write_csv: cannot open the output file";
        return GIA_E_ARG;
    }
    fprintf(f, "time");
    for (i = 0; i < s->n_couples; i++) {
        const char *a = m->nodes[s->couples[i].from].id, *b = m->nodes[s->couples[i].to].id;
        fprintf(f, ",%s__%s_re,%s__%s_im", a, b, a, b);
    }
    fprintf(f, ",R_H\n");
    full.N       = m->n_nodes;
    full.a       = (double complex *)calloc((size_t)m->n_nodes * (size_t)m->n_nodes, sizeof(double complex));
    full.related = (unsigned char *)calloc((size_t)m->n_nodes * (size_t)m->n_nodes, 1);
    if (!full.a || !full.related) {
        gia_matrioska_free(&full); fclose(f); remove(path);
        if (why) *why = "gia_mop_write_csv: out of memory";
        return GIA_E_NOMEM;
    }
    dt = m->t_end / (double)steps;
    for (r = 0; r <= steps; r++) {
        const double t = (double)r * dt;
        fprintf(f, "%.6f", t);
        for (i = 0; i < s->n_couples; i++) {
            const int c = s->couples[i].from * m->n_nodes + s->couples[i].to;
            if ((st = gia_mop_couple(&s->couples[i].beta, s->k, t, &al, why)) != GIA_OK) {
                gia_matrioska_free(&full);
                fclose(f);
                remove(path);
                return st;
            }
            fprintf(f, ",%.17g,%.17g", creal(al), cimag(al));
            full.a[c] = al; full.related[c] = 1;
        }
        write_rh(f, m, s->ref, &full);
        fprintf(f, "\n");
    }
    gia_matrioska_free(&full);
    if (fclose(f) != 0) {
        remove(path);
        if (why) *why = "gia_mop_write_csv: the output file could not be written";
        return GIA_E_ARG;
    }
    return GIA_OK;
}
