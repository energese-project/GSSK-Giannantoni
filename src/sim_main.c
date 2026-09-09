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
    printf("  components        %d\n", m->n_nodes);
    printf("  relationships     %d\n", m->n_edges);
    printf("  horizon t_val     %g\n", m->t_end);
    printf("  derivative order  %d\n", m->order);
    printf("  generative mode   %s\n", m->generative ? "on" : "off");

    ordinality = gia_ordinality(m);
    printf("\n  ordinality        %.3f  (%s)\n", ordinality,
           ordinality >= 1.0 ? "at maximum" : "below maximum");

    printf("\n  component            kind          phi degree   calculi\n");
    for (i = 0; i < m->n_nodes; i++) {
        const gia_node *nd = &m->nodes[i];
        bool            df = gia_drift_free(&nd->phi);
        printf("  %-20s %-13s %6d       %s\n",
               nd->id, gia_node_kind_name(nd->kind), nd->phi.degree,
               df ? "agree" : "TDC drifts");
    }
    printf("    'agree' means phi is affine, so the incipient and\n");
    printf("    traditional derivatives coincide identically. 'TDC drifts'\n");
    printf("    means phi carries curvature and the Bell expansion picks up\n");
    printf("    terms the incipient derivative does not.\n");
}

static void report_duet(const gia_model *m) {
    int i;
    for (i = 0; i < m->n_nodes; i++) {
        const gia_node *nd = &m->nodes[i];
        gia_net         d;
        if (nd->kind != GIA_NODE_INTERACTION) continue;

        /* The half-order incipient derivative of an interaction does not
         * return a function but a couple -- the binary/duet of the 2006
         * paper, carried as one state rather than two candidates. */
        d = gia_incipient_fractional(&nd->phi, 1, 2, m->t_end);
        printf("\n  duet at '%s', order 1/2, t = %g\n", nd->id, m->t_end);
        printf("    branch +   %+.6f %+.6fi\n",
               creal(d.branch[0]), cimag(d.branch[0]));
        printf("    branch -   %+.6f %+.6fi\n",
               creal(d.branch[1]), cimag(d.branch[1]));
        printf("    sum        %+.3e %+.3ei  (branches cancel)\n",
               creal(gia_net_sum(&d)), cimag(gia_net_sum(&d)));
    }
}

int main(int argc, char **argv) {
    const char *model_path = NULL;
    const char *csv_path   = DEFAULT_CSV;
    bool        show_table = false;
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
        } else if (!strcmp(argv[i], "--print")) {
            show_table = true;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            printf("usage: %s [model.json] [--csv PATH] [--out PATH] "
                   "[--steps N] [--print] [--seed PATH]\n", argv[0]);
            printf("  model.json defaults to %s\n", DEFAULT_MODEL);
            printf("  --print  also render the trajectories on stdout\n");
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

    if (!gia_model_load(&model, root)) goto done;

    report_model(&model);
    report_duet(&model);

    /* ---- Mode 1: functional. Numbers move, the graph does not. ---- */
    printf("\n%s\n", RULE);
    printf(" MODE 1 -- FUNCTIONAL: TRAJECTORIES OVER FIXED TOPOLOGY\n");
    printf("%s\n", RULE);
    if (!gia_write_trajectories(&model, csv_path, steps)) goto done;
    printf("  wrote %s  (%d steps to t = %g)\n", csv_path, steps, model.t_end);
    printf("  columns per component: _Q (network), _Em (empower), _Tr\n");
    printf("  (transformity), then _idc, _tdc, _drift; plus psi_network,\n");
    printf("  conservation and emergy_excess for the run as a whole\n");
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
