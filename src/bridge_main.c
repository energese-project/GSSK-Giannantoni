/* bridge_main.c — one GSSK-schema model through both engines (TODO.md, the
 * bridge item; srs FR-BRG-001).
 *
 *   gia_bridge model.json [--method euler|rk4|expm|auto|adaptive] [--csv PATH]
 *
 * The kernel (gssk.h) steps the model with the chosen method; the Giannantoni
 * engine (engine.h) solves the projection of the same model, gia_project, at
 * the kernel's step times, by Q(t) = exp(A t) Q(0). Every node both engines
 * hold is compared, by id.
 *
 * Both trajectories are classical: the Giannantoni engine's network solution
 * is the matrix exponential (E3), not an incipient one. So the difference this
 * reports is the kernel integrator's error against the exact exponential --
 * not the incipient drift -- and the report labels it so.
 *
 * Only a model the projection carries whole is compared. A partial projection
 * is a different model (the dropped pathways change the equations of the nodes
 * they touched), so its difference would be the missing pathways, not
 * integration error: it is refused, exit 2, with the coverage report.
 *
 * This is the one binary that links both engines; neither engine's library
 * units include the other's header (ADR 0011, NFR-SEP-001).
 *
 * Exit: 0 compared, 1 load or argument error, 2 refused.
 */

#include "engine.h"
#include "gssk.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EXIT_REFUSED 2

static char *read_file(const char *path) {
    FILE *f = fopen(path, "rb");
    long  n;
    char *buf;
    if (!f) return NULL;
    if (fseek(f, 0, SEEK_END) != 0 || (n = ftell(f)) < 0) { fclose(f); return NULL; }
    rewind(f);
    buf = (char *)malloc((size_t)n + 1);
    if (buf && fread(buf, 1, (size_t)n, f) != (size_t)n) { free(buf); buf = NULL; }
    if (buf) buf[n] = '\0';
    fclose(f);
    return buf;
}

static int method_known(const char *m) {
    static const char *const ok[] = {"euler", "rk4", "expm", "incipient", "auto", "adaptive", NULL};
    int i;
    for (i = 0; ok[i]; i++) if (!strcmp(ok[i], m)) return 1;
    return 0;
}

/* |a - b| relative to the larger magnitude; 0 when both are 0. */
static double rel_diff(double a, double b) {
    const double s = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
    return s > 0.0 ? fabs(a - b) / s : 0.0;
}

int main(int argc, char **argv) {
    const char    *model_path = NULL, *csv_path = NULL, *method = NULL;
    char          *text = NULL, *kjson = NULL;
    cJSON         *root = NULL, *mop = NULL, *cfg;
    GSSK_Instance *inst = NULL;
    gia_model      gm;
    gia_coverage   cov;
    FILE          *csv = NULL;
    int           *gidx = NULL, ncmp = 0, i, s, steps, rc = EXIT_FAILURE, loaded = 0;
    double        *q = NULL, *max_abs = NULL, *max_rel = NULL, *max_mag = NULL;
    double         t0, t_end, dt, worst = 0.0, worst_scaled = 0.0;
    size_t         nk;

    memset(&gm, 0, sizeof gm);
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--method") && i + 1 < argc) method = argv[++i];
        else if (!strcmp(argv[i], "--csv") && i + 1 < argc) csv_path = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            printf("usage: %s model.json [--method euler|rk4|expm|auto|adaptive] [--csv PATH]\n"
                   "  runs a GSSK-schema model through the kernel and, projected, through the\n"
                   "  Giannantoni engine, and compares every node both hold. exit 0 compared,\n"
                   "  1 load or argument error, 2 refused (the projection is partial)\n", argv[0]);
            return EXIT_SUCCESS;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "gia_bridge: unknown option %s\n", argv[i]);
            return EXIT_FAILURE;
        } else if (!model_path) model_path = argv[i];
        else { fprintf(stderr, "gia_bridge: unexpected argument %s\n", argv[i]); return EXIT_FAILURE; }
    }
    if (!model_path) { fprintf(stderr, "gia_bridge: no model given (--help)\n"); return EXIT_FAILURE; }
    if (method && !method_known(method)) {
        fprintf(stderr, "gia_bridge: unknown --method '%s' (euler, rk4, expm, auto, adaptive)\n", method);
        return EXIT_FAILURE;
    }
    if (!(text = read_file(model_path))) {
        fprintf(stderr, "gia_bridge: cannot read %s\n", model_path);
        goto done;
    }
    if (!(root = cJSON_Parse(text))) {
        fprintf(stderr, "gia_bridge: %s is not valid JSON\n", model_path);
        goto done;
    }

    /* ---- The Giannantoni side: project, and refuse a partial projection. ---- */
    if (!gia_project(root, &mop, &cov)) {
        fprintf(stderr, "gia_bridge: %s does not look like a GSSK-schema model\n", model_path);
        goto done;
    }
    gia_report_coverage(&cov);
    printf("bridge.coverage: %.1f%%\n", 100.0 * gia_coverage_fraction(&cov));
    if (gia_coverage_fraction(&cov) < 1.0) {
        fprintf(stderr, "gia_bridge: refused: the projection carries %.1f%% of %s (coverage report "
                        "above). A partial projection is a different model, so a difference would "
                        "measure the missing pathways, not integration error.\n",
                100.0 * gia_coverage_fraction(&cov), model_path);
        rc = EXIT_REFUSED;
        goto done;
    }
    if (!gia_model_load(&gm, mop)) {
        fprintf(stderr, "gia_bridge: the projection does not load in the Giannantoni engine\n");
        goto done;
    }
    loaded = 1;

    /* ---- The kernel side, with the method asked for. ---- */
    cfg = cJSON_GetObjectItemCaseSensitive(root, "config");
    if (method) {
        if (!cJSON_IsObject(cfg)) cfg = cJSON_AddObjectToObject(root, "config");
        if (!cfg) goto done;
        cJSON_DeleteItemFromObjectCaseSensitive(cfg, "method");
        cJSON_AddStringToObject(cfg, "method", method);
    }
    if (!(kjson = cJSON_PrintUnformatted(root)) || GSSK_Init(kjson, &inst) != GSSK_SUCCESS) {
        fprintf(stderr, "gia_bridge: the kernel does not load %s\n", model_path);
        goto done;
    }
    t0    = GSSK_GetTStart(inst);
    t_end = GSSK_GetTEnd(inst);
    dt    = GSSK_GetDt(inst);
    steps = (int)floor((t_end - t0) / dt + 0.5);
    nk    = GSSK_GetStateSize(inst);

    /* Nodes both engines hold, in the kernel's order. */
    gidx    = (int *)malloc(nk * sizeof(int));
    max_abs = (double *)calloc(nk, sizeof(double));
    max_rel = (double *)calloc(nk, sizeof(double));
    max_mag = (double *)calloc(nk, sizeof(double));
    q       = (double *)calloc((size_t)gm.n_nodes, sizeof(double));
    if (!gidx || !max_abs || !max_rel || !max_mag || !q) goto done;
    for (i = 0; i < (int)nk; i++) {
        const char *id = GSSK_GetNodeID(inst, (size_t)i);
        int         j;
        gidx[i] = -1;
        for (j = 0; j < gm.n_nodes && id; j++)
            if (gm.nodes[j].id && !strcmp(gm.nodes[j].id, id) && !gm.nodes[j].is_module) gidx[i] = j;
        if (gidx[i] >= 0) ncmp++;
    }

    if (csv_path) {
        if (!(csv = fopen(csv_path, "w"))) {
            fprintf(stderr, "gia_bridge: cannot write %s\n", csv_path);
            goto done;
        }
        fprintf(csv, "time");
        for (i = 0; i < (int)nk; i++) {
            const char *id = GSSK_GetNodeID(inst, (size_t)i);
            if (gidx[i] >= 0) fprintf(csv, ",%s_kernel,%s_gia,%s_diff", id, id, id);
        }
        fprintf(csv, "\n");
    }
    for (s = 0; s <= steps; s++) {
        const double  t = s * dt;
        const double *qk;
        if (s > 0 && GSSK_Step(inst, dt) != GSSK_SUCCESS) {
            fprintf(stderr, "gia_bridge: the kernel failed at step %d\n", s);
            goto done;
        }
        qk = GSSK_GetState(inst);
        if (!gia_network_state(&gm, t, q, NULL)) {
            fprintf(stderr, "gia_bridge: the Giannantoni engine did not solve at t = %g\n", t);
            goto done;
        }
        if (csv) fprintf(csv, "%.6f", t0 + t);
        /* A forced held node (a source driven by its waveform, ADR 0006) is
         * reported by the kernel at its forced value; this engine holds the
         * declared value in its state and applies the waveform where the
         * node is read. Compare like with like. */
        for (i = 0; i < (int)nk; i++) {
            int j = gidx[i];
            if (j >= 0 && !gm.nodes[j].integrates)
                q[j] = gia_forcing_value(&gm.nodes[j].forcing, q[j], t);
        }
        for (i = 0; i < (int)nk; i++) {
            double a, r;
            if (gidx[i] < 0) continue;
            a = fabs(qk[i] - q[gidx[i]]);
            r = rel_diff(qk[i], q[gidx[i]]);
            if (a > max_abs[i]) max_abs[i] = a;
            if (fabs(q[gidx[i]]) > max_mag[i]) max_mag[i] = fabs(q[gidx[i]]);
            if (r > max_rel[i]) max_rel[i] = r;
            if (r > worst) worst = r;
            if (csv) fprintf(csv, ",%.17g,%.17g,%.17g", qk[i], q[gidx[i]], qk[i] - q[gidx[i]]);
        }
        if (csv) fprintf(csv, "\n");
    }

    printf("\nbridge.model: %s\n", model_path);
    printf("bridge.kernel_method: %s\n", method ? method : "(the model's own)");
    printf("bridge.compared: %d nodes, %d steps of dt = %g\n", ncmp, steps, dt);
    for (i = 0; i < (int)nk; i++) {
        if (gidx[i] < 0) continue;
        printf("bridge.max_abs_difference.%s: %.6e\n", GSSK_GetNodeID(inst, (size_t)i), max_abs[i]);
        printf("bridge.max_rel_difference.%s: %.17g\n", GSSK_GetNodeID(inst, (size_t)i), max_rel[i]);
        /* Scaled by the node's own size over the run: a store that is
         * exactly 0 on one side and 1e-9 on the other is a relative
         * difference of 1, and a scaled one of nearly nothing. */
        {
            const double sc = max_mag[i] > 0.0 ? max_abs[i] / max_mag[i] : max_abs[i];
            printf("bridge.max_scaled_difference.%s: %.17g\n", GSSK_GetNodeID(inst, (size_t)i), sc);
            if (sc > worst_scaled) worst_scaled = sc;
        }
    }
    printf("bridge.max_rel_difference: %.17g\n", worst);
    printf("bridge.max_scaled_difference: %.17g\n", worst_scaled);
    /* E3: both trajectories are classical. */
    printf("label.bridge_difference: classical\n");
    printf("  (the kernel's integration error against the exact exponential, not incipient drift:\n"
           "   the Giannantoni engine's network trajectories are exp(A t), PLAN.md §2 E3)\n");
    rc = EXIT_SUCCESS;

done:
    if (csv && fclose(csv) != 0 && rc == EXIT_SUCCESS) rc = EXIT_FAILURE;
    if (rc != EXIT_SUCCESS && csv_path && csv) remove(csv_path);
    if (inst) GSSK_Free(inst);
    if (loaded) gia_model_free(&gm);
    free(gidx); free(max_abs); free(max_rel); free(max_mag); free(q);
    if (kjson) free(kjson);
    if (mop) cJSON_Delete(mop);
    if (root) cJSON_Delete(root);
    free(text);
    return rc;
}
