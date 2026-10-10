/* sim_main.c — CLI for the Giannantoni generative simulator.
 *
 *   giannantoni_sim [model.json] [--csv PATH] [--out PATH] [--steps N] [--print]
 *                   [--mop-out PATH]
 *
 * Exit status (IF-CLI-001): 0 on success, 1 on a load or validation error, 2
 * when the kernel refuses to compute (FR-OUT-002), with the reason on stderr.
 *
 * The model path defaults to examples/giannantoni/input.json. One run does
 * both modes:
 * the functional trajectories go to CSV, the generative ordinal step goes to
 * JSON, and the validator then decides from the two graphs which of the two
 * actually happened.
 *
 * The entry point is kept out of engine.c so the engine can be linked into
 * tests without an entry point coming with it.
 */

#include "engine.h"
#include "mop_seed.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEFAULT_MODEL "examples/giannantoni/input.json"
#define DEFAULT_CSV   "simulation_output.csv"
#define DEFAULT_OUT   "output.json"
#define DEFAULT_STEPS 10

static const char *RULE =
    "============================================================";

static char *read_file(const char *path) {
    FILE  *f;
    long   len;
    size_t got;
    char  *buf;

    f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "cannot open %s: ", path);
        perror(NULL);
        return NULL;
    }
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    len = ftell(f);
    if (len < 0) { fclose(f); return NULL; }
    rewind(f);

    buf = (char *)malloc((size_t)len + 1);
    if (!buf) { fclose(f); return NULL; }

    got = fread(buf, 1, (size_t)len, f);
    fclose(f);
    if (got != (size_t)len) {
        fprintf(stderr, "short read on %s\n", path);
        free(buf);
        return NULL;
    }
    buf[len] = '\0';
    return buf;
}

/* Derive "<out>.seed.json" from "<out>.json".
 *
 * This exists because a text diff of the ORIGINAL input file against the
 * output is not meaningful: cJSON reformats on the way out (indentation
 * changes, 2.0 prints as 2), so even a run that cJSON_Compare calls identical
 * shows a large textual diff. Re-serialising the seed through the same
 * cJSON_Print path removes that noise, so what a diff tool then shows is
 * exactly the structural change and nothing else. */
static bool derive_seed_path(const char *out, char *buf, size_t cap) {
    const char *slash = strrchr(out, '/');
    const char *dot   = strrchr(out, '.');
    size_t      stem;

    if (dot && (!slash || dot > slash)) {
        stem = (size_t)(dot - out);
        if (stem + strlen(".seed") + strlen(dot) + 1 > cap) return false;
        memcpy(buf, out, stem);
        buf[stem] = '\0';
        strcat(buf, ".seed");
        strcat(buf, dot);
    } else {
        if (strlen(out) + strlen(".seed.json") + 1 > cap) return false;
        strcpy(buf, out);
        strcat(buf, ".seed.json");
    }
    return true;
}

static bool write_file(const char *path, const char *text) {
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "cannot write %s: ", path);
        perror(NULL);
        return false;
    }
    fputs(text, f);
    fputc('\n', f);
    fclose(f);
    return true;
}

/* The reference couple alpha_12 is read off the graph rather than invented:
 * its modulus is the weight of the first relationship in the seed, and its
 * argument is the generativity factor of the component that relationship
 * feeds. A graph with no edges has no couple to reference, and falls back to
 * a real unit. */
static double complex reference_couple(const gia_model *m) {
    double w = 1.0, g = 0.0;
    if (m->n_edges > 0) {
        const gia_edge *e0 = &m->edges[0];
        w = e0->weight;
        if (e0->to >= 0 && e0->to < m->n_nodes) {
            const gia_node *tgt = &m->nodes[e0->to];
            /* phi = g t^2 / 2 for the generative kinds, so g = 2 c[2]. */
            g = 2.0 * tgt->phi.c[2];
        }
    }
    return w * (cos(g) + I * sin(g));
}

/* IF-CLI-001: a refusal -- the kernel will not compute what is asked, and
 * says why -- exits 2; anything else that stops the run exits 1. */
#define EXIT_REFUSED 2

static int is_refusal(gia_status st) {
    return st == GIA_E_DOMAIN || st == GIA_E_UNSUPPORTED || st == GIA_E_RANGE ||
           st == GIA_E_CONVERGENCE || st == GIA_E_LIMIT;
}

/* FR-OUT-002: the feature, then the reason the library gave, which names its
 * source. Returns the exit status. */
static int report_failure(const char *feature, gia_status st, const char *why) {
    fprintf(stderr, "giannantoni_sim: %s %s: %s (%s)\n", feature,
            is_refusal(st) ? "refused" : "failed", why ? why : "no reason given",
            gia_status_str(st));
    return is_refusal(st) ? EXIT_REFUSED : EXIT_FAILURE;
}

static int component_total(const gia_model *m) {
    int i, n = 0;
    for (i = 0; i < m->n_nodes; i++) if (gia_node_is_component(m, i)) n++;
    return n;
}

/* The seed's `mop` block (IF-JSON-001): what it asks for, and the parts that
 * are evaluated at t_end. Returns 0, or the exit status of a refusal. */
static int report_mop(const gia_model *m, const gia_mop_seed *s) {
    gia_status  st;
    const char *why = NULL;

    printf("\n%s\n", RULE);
    printf(" MOP -- THE FUNDAMENTAL EQUATIONS (IF-JSON-001)\n");
    printf("%s\n", RULE);
    printf("  cardinality k     %d/%d\n", s->k.num, s->k.den);
    printf("  reference couple  %s -> %s\n", m->nodes[s->ref[0]].id, m->nodes[s->ref[1]].id);
    if (s->form == GIA_MOP_BETA_NETWORK)
        printf("  beta              from the network (FR-MOP-008)\n");
    else
        printf("  related couples   %d  (First Equation, [23 Eq 5.4.2])\n", s->n_couples);
    if (s->has_second) {
        gia_second sec;
        st = gia_mop_second(s->alpha12_0, s->c1, s->c2, component_total(m), m->t_end, &sec,
                            NULL, &why);
        if (st != GIA_OK) return report_failure("second_equation (FR-MOP-005)", st, why);
        printf("  second equation   A(t_end) = %.10g%+.10gi  ([23 Eq 6.3])\n",
               creal(sec.A), cimag(sec.A));
    }
    if (s->has_eqs) {
        /* The EQS needs the reference couple's relational coordinates, which
         * the seed does not supply; the parameters are checked here, so X11
         * is refused at load rather than silently carried. */
        const rel_t probe = {1.0, 0.0, 0.0};
        double      out[3];
        st = gia_eqs(&s->eqs, probe, 1, out, &why);
        if (st != GIA_OK) return report_failure("eqs (FR-MOP-006)", st, why);
        printf("  eqs               parameters accepted (N = %d); the coordinates need\n"
               "                    the reference couple's {Sigma0, Phi0, Theta0}\n", s->eqs.N);
    }
    return 0;
}

/* FR-OUT-003, FR-HAR-002: the detector's verdict for each construction the
 * sources give, at N = the number of components, and the observed verdict of
 * the network's own Matrioska (FR-MOP-008) on the reference row. With fewer
 * than three components [23 Eq 5.6.5] has no residual, so nothing is
 * evaluated, and the report says so rather than printing a verdict. */
static void report_harmony(const gia_model *m, const gia_mop_seed *mop) {
    static const struct { const char *name; gia_construction c; } cs[3] = {
        {"first_equation", gia_construct_first},
        {"second_equation", gia_construct_second},
        {"eqs", gia_construct_eqs},
    };
    const int     N = component_total(m);
    gia_verdict   v;
    gia_matrioska al, row;
    const char   *why = NULL;
    int           i;

    printf("\n%s\n", RULE);
    printf(" HARMONY -- DETECTED, NOT ASSUMED (FR-HAR-002, [23 Eq 5.6.5])\n");
    printf("%s\n", RULE);
    if (N < 3) {
        printf("  not evaluated: %d component%s; the residual needs N >= 3\n", N,
               N == 1 ? "" : "s");
        return;
    }
    for (i = 0; i < 3; i++) {
        void *ctx = (i == 0 && mop->present) ? (void *)&mop->k : NULL;
        why = NULL;
        if (gia_harmony_verdict(cs[i].c, ctx, N, &v, &why) == GIA_OK)
            printf("harmony.%s: %s\n", cs[i].name, gia_verdict_str(v));
        else
            printf("  %s: not evaluated (%s)\n", cs[i].name, why ? why : "refused");
    }
    /* The network's Matrioska, observed on the reference couple's row. A row
     * couple without a pathway leaves the residual undefined, and a Matrioska
     * that does not relate the couples harmony relates is not harmonic. */
    memset(&al, 0, sizeof al);
    memset(&row, 0, sizeof row);
    if (gia_mop_network(m, m->t_end, &al, NULL, &why) != GIA_OK) {
        printf("  network: not evaluated (%s)\n", why ? why : "refused");
        return;
    }
    if (gia_mop_reference_row(m, mop->present ? mop->ref : NULL, &al, &row, &why) == GIA_OK &&
        gia_harmony_observed(&row, &v, &why) == GIA_OK)
        printf("harmony.network: %s\n", gia_verdict_str(v));
    else
        printf("harmony.network: absent\n  (no residual: %s)\n", why ? why : "?");
    gia_matrioska_free(&row);
    gia_matrioska_free(&al);
}

static void report_model(gia_model *m) {
    gia_ordinality_rec rec;
    const char        *why = NULL;
    int                i;

    printf("%s\n", RULE);
    printf(" GIANNANTONI GENERATIVE FRAMEWORK\n");
    printf("%s\n", RULE);
    printf("  system            %s\n", m->system_name);
    /* Modules are counted apart: a module holds nothing, so it is not a
     * component and does not enter ordinality (ADR 0014). */
    printf("  components        %d\n", gia_component_count(m));
    if (gia_component_count(m) != m->n_nodes)
        printf("  modules           %d\n", m->n_nodes - gia_component_count(m));
    printf("  relationships     %d\n", m->n_edges);
    printf("  horizon t_val     %g\n", m->t_end);
    printf("  derivative order  %d\n", m->order);
    printf("  generative mode   %s\n", m->generative ? "on" : "off");

    /* IF-OUT-003, ADR 0021: the record {k, n22, n2, nhalf, nunrelated}, the
     * verdict (strong connectivity), and the old fraction as a labelled
     * proxy that decides nothing. */
    printf("\n");
    if (gia_ordinality_record(m, &rec, &why) == GIA_OK)
        printf("ordinality: {%d, %d, %d, %d, %d}\n", rec.k, rec.n22, rec.n2, rec.nhalf,
               rec.nunrelated);
    else
        printf("ordinality: not computed (%s)\n", why ? why : "?");
    printf("maximum_ordinality: %s\n", gia_at_maximum_ordinality(m) ? "yes" : "no");
    printf("closure (proxy): %.3f\n", gia_closure(m));

    /* ADR 0017 decision 4: a control at use_ratio 0 is drawn as a pathway and
     * carries nothing, so its energy and its emergy are left out. That is
     * allowed -- a price is not dissipated -- but it is never silent. */
    for (i = 0; i < m->n_edges; i++) {
        const gia_edge *e = &m->edges[i];
        if (e->role != GIA_ROLE_CONTROL || e->use_ratio > 0.0) continue;
        if (e->from < 0 || e->to < 0) continue;
        printf("  control neglected %s -> %s  (use_ratio 0: read only, its "
               "energy and emergy left out)\n",
               m->nodes[e->from].id, m->nodes[e->to].id);
    }

    /* Solution drift (FR-IDC-011): exactly zero for a constant flow matrix,
     * undefined by the sources otherwise. The per-node "calculi agree / TDC
     * drifts" table that stood here judged a phi the engine invented, not the
     * trajectory (PLAN E4), and is gone. */
    {
        double     *zero = (double *)calloc((size_t)m->n_nodes, sizeof(double));
        const char *why  = NULL;
        if (zero && gia_solution_drift(m, m->order >= 1 && m->order <= 4 ? m->order : 2,
                                       zero, &why) == GIA_OK)
            printf("  solution drift    0 for every component  (constant flow "
                   "matrix, [06 §4 (i)])\n");
        else
            printf("  solution drift    not defined: %s\n", why ? why : "out of memory");
        free(zero);
    }
}

/* FR-OUT-001: every output says what produced it. One `label.<name>: <label>`
 * line per CSV column, in header order, then one per reported quantity. The
 * five labels (docs/requirements/srs.md FR-OUT-001):
 *
 *   implemented   the source's construction, verified against it
 *   classical     a correct result of ordinary (TDC) numerics, not of IDC
 *   assumed       built from an assumption the sources state, not derived
 *   proxy         a stand-in for a quantity the sources define differently
 *   illustrative  the source's formula applied to inputs it does not supply
 *
 * PLAN.md §2 is the evidence for each assignment. tests/mop_cli.sh checks the
 * lines against the CSV header the same run writes, so a column cannot be
 * added here or in gia_write_trajectories without the other. */
static void report_labels(const gia_model *m, const gia_mop_seed *mop) {
    int i;

    printf("\n%s\n", RULE);
    printf(" WHAT PRODUCED EACH OUTPUT (FR-OUT-001)\n");
    printf("%s\n", RULE);
    printf("label.time: implemented\n");
    for (i = 0; i < m->n_nodes; i++) {
        const char *id = m->nodes[i].id;
        /* E3: Q(t) = exp(A t) Q(0) is the matrix exponential. */
        printf("label.%s_Q: classical\n", id);
        /* Odum's emergy algebra, [02 p. 23 rules 1-4] (FR-EM-001..003). */
        printf("label.%s_Em: implemented\n", id);
        printf("label.%s_Tr: implemented\n", id);
        /* FR-IDC-014: [09 Eq 13] at k = 2 along the solved trajectory. */
        printf("label.%s_drift_proj: implemented\n", id);
    }
    printf("label.conservation: classical\n");
    printf("label.emergy_excess: implemented\n");

    /* ADR 0021: the record of [23 Eq 3.2] and the verdict of [22 Eq 11.1]. */
    printf("label.ordinality: implemented\n");
    printf("label.maximum_ordinality: implemented\n");
    /* E6: the fraction of components on a cycle, not [22 Eq 11.1]. */
    printf("label.closure: proxy\n");
    /* E5: constructed from roots of unity and checked against the
     * construction; [23 §8 ii] says the EQS assumes it. */
    printf("label.harmony: assumed\n");
    /* E7: ADR 0021 §2, maximum total empower, the discrete [02 Eq 5.3]. */
    printf("label.generative_step: implemented\n");
    /* FR-IDC-011: exactly zero for a constant flow matrix, refused otherwise. */
    printf("label.solution_drift: implemented\n");
    if (mop->present) {
        /* ADR 0019: the First Equation solved from [23 Eq 5.5.6], not the
         * printed Eq 5.5.7-5.5.8 (X1); written by --mop-out. */
        printf("label.mop_alpha: implemented\n");
        /* [23 §8 ii]: the printed solution takes lambda null, an assumption. */
        if (mop->has_second) printf("label.second_equation: assumed\n");
    }
}

int main(int argc, char **argv) {
    const char *model_path = NULL;
    const char *csv_path   = DEFAULT_CSV;
    bool        show_table = false;
    bool        project    = false;
    const char *out_path   = DEFAULT_OUT;
    const char *seed_path  = NULL;
    const char *mop_path   = NULL;
    int         steps      = DEFAULT_STEPS;
    int         i, rc      = EXIT_FAILURE;

    char        *text   = NULL;
    cJSON       *root   = NULL;
    cJSON       *out    = NULL;
    char        *printed = NULL;
    char        *seed_txt = NULL;
    char         seed_buf[1024];
    gia_model    model;
    gia_harmony  harmony;
    gia_mode     mode;
    gia_mop_seed mop;

    memset(&model,   0, sizeof(model));
    memset(&harmony, 0, sizeof(harmony));
    memset(&mop,     0, sizeof(mop));

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--csv") && i + 1 < argc) {
            csv_path = argv[++i];
        } else if (!strcmp(argv[i], "--out") && i + 1 < argc) {
            out_path = argv[++i];
        } else if (!strcmp(argv[i], "--steps") && i + 1 < argc) {
            steps = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed_path = argv[++i];
        } else if (!strcmp(argv[i], "--mop-out") && i + 1 < argc) {
            mop_path = argv[++i];
        } else if (!strcmp(argv[i], "--project")) {
            project = true;
        } else if (!strcmp(argv[i], "--print")) {
            show_table = true;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            printf("usage: %s [model.json] [--csv PATH] [--out PATH] "
                   "[--steps N] [--print] [--seed PATH] [--mop-out PATH]\n", argv[0]);
            printf("  model.json defaults to %s\n", DEFAULT_MODEL);
            printf("  --print  also render the trajectories on stdout\n");
            printf("  --project  read the file as a GSSK-schema model, report\n");
            printf("           what fraction of it this engine can carry, and\n");
            printf("           run the projection\n");
            printf("  --seed   where to write the re-serialised seed for\n");
            printf("           diffing (default: <out> with .seed before .json)\n");
            printf("  --mop-out  write the First Equation's couples (the seed's\n");
            printf("           mop block) as CSV\n");
            printf("  exit status: 0 success, 1 load or validation error,\n");
            printf("           2 refused (the reason is on stderr)\n");
            return EXIT_SUCCESS;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "unknown option: %s\n", argv[i]);
            return EXIT_FAILURE;
        } else if (!model_path) {
            model_path = argv[i];
        } else {
            fprintf(stderr, "unexpected argument: %s\n", argv[i]);
            return EXIT_FAILURE;
        }
    }
    if (!model_path) model_path = DEFAULT_MODEL;
    if (steps < 1)   steps = 1;

    text = read_file(model_path);
    if (!text) goto done;

    root = cJSON_Parse(text);
    if (!root) {
        const char *err = cJSON_GetErrorPtr();
        fprintf(stderr, "%s is not valid JSON", model_path);
        if (err) fprintf(stderr, " (near: %.32s)", err);
        fprintf(stderr, "\n");
        goto done;
    }

    if (project) {
        /* ADR 0011: one-way, and it declares what it drops. The report names
         * the blocking modules rather than publishing a bare percentage,
         * because a score below 100%% is a fact about this model and not a
         * verdict on the framework. */
        cJSON        *projected = NULL;
        gia_coverage  cov;
        if (!gia_project(root, &projected, &cov)) {
            fprintf(stderr, "%s does not look like a GSSK-schema model\n",
                    model_path);
            goto done;
        }
        gia_report_coverage(&cov);
        cJSON_Delete(root);
        root = projected;          /* run what could be carried */
    }

    if (!gia_model_load(&model, root)) goto done;

    {
        char        detail[256];
        const char *why = NULL;
        gia_status  st  = gia_mop_seed_load(&model, &mop, detail, sizeof detail, &why);
        if (st != GIA_OK) {
            fprintf(stderr, "%s: %s\n", model_path, detail[0] ? detail : (why ? why : "mop"));
            goto done;
        }
    }
    if (mop_path && !mop.present) {
        fprintf(stderr, "--mop-out: %s has no mop block (IF-JSON-001)\n", model_path);
        goto done;
    }

    report_model(&model);

    if (mop.present) {
        int refused = report_mop(&model, &mop);
        if (refused) { rc = refused; goto done; }
    }
    if (mop_path) {
        const char *why = NULL;
        gia_status  st  = gia_mop_write_csv(&model, &mop, mop_path, steps, &why);
        if (st != GIA_OK) { rc = report_failure("--mop-out", st, why); goto done; }
        printf("  wrote %s  (%d couples, %d steps to t = %g)\n", mop_path, mop.n_couples,
               steps, model.t_end);
    }

    /* ---- Mode 1: functional. Numbers move, the graph does not. ---- */
    printf("\n%s\n", RULE);
    printf(" MODE 1 -- FUNCTIONAL: TRAJECTORIES OVER FIXED TOPOLOGY\n");
    printf("%s\n", RULE);
    if (!gia_write_trajectories(&model, csv_path, steps)) goto done;
    printf("  wrote %s  (%d steps to t = %g)\n", csv_path, steps, model.t_end);
    printf("  columns per node: _Q (network), _Em (empower), _Tr\n");
    printf("  (transformity), _drift_proj (output-projection drift, [09 Eq 13]);\n");
    printf("  plus conservation and emergy_excess for the run as a whole\n");
    if (show_table) gia_print_trajectories(&model, steps);

    /* ---- The harmony CONSTRUCTOR: assumed, and labelled so (FR-HAR-004) ---- */
    if (model.n_nodes >= 2) {
        if (!gia_harmony_assume_init(&harmony, model.n_nodes,
                              reference_couple(&model))) {
            fprintf(stderr, "cannot build harmony matrix\n");
            goto done;
        }
        gia_validate_harmony(&harmony, 1e-9);
    }

    /* ---- The harmony DETECTOR: one verdict per construction (FR-OUT-003) ---- */
    report_harmony(&model, &mop);

    /* ---- Mode 2: generative. The graph itself may change. ---- */
    printf("\n%s\n", RULE);
    printf(" MODE 2 -- GENERATIVE: MOP ORDINAL STEP\n");
    printf("%s\n", RULE);
    out = gia_generate(&model);
    if (!out) {
        fprintf(stderr, "generative step failed\n");
        goto done;
    }
    printed = cJSON_Print(out);
    if (!printed || !write_file(out_path, printed)) goto done;

    if (!seed_path) {
        if (!derive_seed_path(out_path, seed_buf, sizeof(seed_buf))) {
            fprintf(stderr, "output path too long to derive a seed path\n");
            goto done;
        }
        seed_path = seed_buf;
    }
    seed_txt = cJSON_Print(root);
    if (!seed_txt || !write_file(seed_path, seed_txt)) goto done;

    printf("\n  input   %s\n", model_path);
    printf("  seed    %s\n", seed_path);
    printf("  output  %s\n", out_path);
    printf("\n  compare  code --diff %s %s\n", seed_path, out_path);
    printf("           (diff the seed, not the input file: cJSON reformats on\n");
    printf("            output, so an input-vs-output text diff shows\n");
    printf("            serialisation noise on top of the real change)\n");

    /* ---- Which mode was it, really? Decided from the two graphs. ---- */
    mode = gia_validate_mode(root, out);
    gia_report_mode(mode, root, out);

    report_labels(&model, &mop);

    rc = EXIT_SUCCESS;

done:
    gia_harmony_assume_free(&harmony);
    gia_mop_seed_free(&mop);
    gia_model_free(&model);
    if (printed)  free(printed);
    if (seed_txt) free(seed_txt);
    if (out)  cJSON_Delete(out);
    if (root) cJSON_Delete(root);
    if (text) free(text);
    return rc;
}
