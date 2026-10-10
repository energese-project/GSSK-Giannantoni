/* sim_main.c — CLI for the Giannantoni generative simulator.
 *
 *   giannantoni_sim [model.json] [--csv PATH] [--out PATH] [--steps N] [--print]
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

static void report_model(gia_model *m) {
    double ordinality;
    int    i;

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

    ordinality = gia_ordinality(m);
    printf("\n  ordinality        %.3f  (%s)\n", ordinality,
           ordinality >= 1.0 ? "at maximum" : "below maximum");

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
static void report_labels(const gia_model *m) {
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

    /* E6: the fraction of components on a cycle, not [22 Eq 11.1]. */
    printf("label.ordinality: proxy\n");
    /* E5: constructed from roots of unity and checked against the
     * construction; [23 §8 ii] says the EQS assumes it. */
    printf("label.harmony: assumed\n");
    /* E7: the ADR 0015 graph heuristic, not [02 Eq 5.3] / [22 Eq 2]. */
    printf("label.generative_step: illustrative\n");
    /* FR-IDC-011: exactly zero for a constant flow matrix, refused otherwise. */
    printf("label.solution_drift: implemented\n");
}

int main(int argc, char **argv) {
    const char *model_path = NULL;
    const char *csv_path   = DEFAULT_CSV;
    bool        show_table = false;
    bool        project    = false;
    const char *out_path   = DEFAULT_OUT;
    const char *seed_path  = NULL;
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

    memset(&model,   0, sizeof(model));
    memset(&harmony, 0, sizeof(harmony));

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--csv") && i + 1 < argc) {
            csv_path = argv[++i];
        } else if (!strcmp(argv[i], "--out") && i + 1 < argc) {
            out_path = argv[++i];
        } else if (!strcmp(argv[i], "--steps") && i + 1 < argc) {
            steps = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--seed") && i + 1 < argc) {
            seed_path = argv[++i];
        } else if (!strcmp(argv[i], "--project")) {
            project = true;
        } else if (!strcmp(argv[i], "--print")) {
            show_table = true;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            printf("usage: %s [model.json] [--csv PATH] [--out PATH] "
                   "[--steps N] [--print] [--seed PATH]\n", argv[0]);
            printf("  model.json defaults to %s\n", DEFAULT_MODEL);
            printf("  --print  also render the trajectories on stdout\n");
            printf("  --project  read the file as a GSSK-schema model, report\n");
            printf("           what fraction of it this engine can carry, and\n");
            printf("           run the projection\n");
            printf("  --seed   where to write the re-serialised seed for\n");
            printf("           diffing (default: <out> with .seed before .json)\n");
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

    report_model(&model);

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

    /* ---- MOP harmony relationships over the N components ---- */
    if (model.n_nodes >= 2) {
        if (!gia_harmony_init(&harmony, model.n_nodes,
                              reference_couple(&model))) {
            fprintf(stderr, "cannot build harmony matrix\n");
            goto done;
        }
        gia_validate_harmony(&harmony, 1e-9);
    }

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

    report_labels(&model);

    rc = EXIT_SUCCESS;

done:
    gia_harmony_free(&harmony);
    gia_model_free(&model);
    if (printed)  free(printed);
    if (seed_txt) free(seed_txt);
    if (out)  cJSON_Delete(out);
    if (root) cJSON_Delete(root);
    if (text) free(text);
    return rc;
}
