/* projection.c — a GSSK model read into MOP terms, declaring its losses.
 *
 * ADR 0011 decisions 2 and 3. See include/engine.h section 8c for why the
 * direction is one-way and what the coverage number is for.
 *
 * The rule this file exists to enforce: nothing is silently substituted. Where
 * a component or pathway cannot be carried, it is dropped and a finding names
 * the module that stopped it. A projection that quietly turned a switch into a
 * linear pathway would produce a model that runs, gives numbers, and is not
 * the model anyone wrote.
 */

#include "engine.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

static void finding(gia_coverage *cov, const char *fmt, ...) {
    va_list ap;
    if (!cov) return;
    if (cov->n_findings >= GIA_MAX_FINDINGS) { cov->n_elided++; return; }
    va_start(ap, fmt);
    vsnprintf(cov->finding[cov->n_findings], sizeof(cov->finding[0]), fmt, ap);
    va_end(ap);
    cov->n_findings++;
}

static const char *sfield(const cJSON *o, const char *k, const char *dflt) {
    const cJSON *it = cJSON_GetObjectItemCaseSensitive(o, k);
    return (cJSON_IsString(it) && it->valuestring) ? it->valuestring : dflt;
}

static double nfield(const cJSON *o, const char *k, double dflt) {
    const cJSON *it = cJSON_GetObjectItemCaseSensitive(o, k);
    return cJSON_IsNumber(it) ? it->valuedouble : dflt;
}

/* The MOP engine's node vocabulary is the nine GSSK primitives. Composites are
 * not expanded here (ADR 0010), and a user archetype is a composite by another
 * name. */
static bool node_type_carried(const char *t, const char **why) {
    static const char *prim[] = { "storage", "source", "sink", "constant",
                                  "interaction", "gain", "loop_limited",
                                  "exchange", "switch" };
    size_t i;
    if (!t) { *why = "no type"; return false; }
    for (i = 0; i < sizeof(prim) / sizeof(prim[0]); i++)
        if (!strcmp(t, prim[i])) return true;

    if (!strcmp(t, "producer") || !strcmp(t, "consumer") ||
        !strcmp(t, "misc_box") || !strcmp(t, "system_frame"))
        *why = "a built-in composite; composites are not expanded (ADR 0010)";
    else
        *why = "a user archetype, which is a composite by another name "
               "(ADR 0010)";
    return false;
}

/* Pathway laws. The names match the kernel's deliberately, so the two
 * vocabularies cannot drift; what differs is which the MOP engine implements. */
static bool logic_carried(const char *l, const char **why) {
    static const char *ok[] = { "linear", "interaction", "limit", "threshold",
                                "ratio", "reversible", "subtract", "constant" };
    size_t i;
    if (!l) l = "linear";                     /* the schema's default */
    for (i = 0; i < sizeof(ok) / sizeof(ok[0]); i++)
        if (!strcmp(l, ok[i])) return true;
    *why = "not a pathway law this engine implements";
    return false;
}

bool gia_project(const cJSON *gssk, cJSON **out_mop, gia_coverage *cov) {
    const cJSON *nodes, *edges, *it;
    cJSON       *out, *onodes, *oedges, *params;

    if (!gssk || !out_mop || !cov) return false;
    memset(cov, 0, sizeof(*cov));
    *out_mop = NULL;

    nodes = cJSON_GetObjectItemCaseSensitive(gssk, "nodes");
    edges = cJSON_GetObjectItemCaseSensitive(gssk, "edges");
    if (!cJSON_IsArray(nodes)) return false;

    out = cJSON_CreateObject();
    if (!out) return false;
    cJSON_AddStringToObject(out, "system_name",
        sfield(cJSON_GetObjectItemCaseSensitive(gssk, "metadata"), "name",
               "(projected from a GSSK model)"));
    onodes = cJSON_AddArrayToObject(out, "nodes");
    oedges = cJSON_AddArrayToObject(out, "edges");
    if (!onodes || !oedges) { cJSON_Delete(out); return false; }

    cJSON_ArrayForEach(it, nodes) {
        const char *id   = sfield(it, "id", NULL);
        const char *type = sfield(it, "type", NULL);
        const char *why  = NULL;
        cJSON      *n;

        cov->nodes_total++;
        if (!node_type_carried(type, &why)) {
            finding(cov, "node '%s': type '%s' not carried — %s",
                    id ? id : "?", type ? type : "(none)", why);
            continue;
        }
        /* GSSK configures a PROCESSING node's law in its own `params` block
         * (Phase 7: k, C, threshold, price). This engine puts laws on
         * pathways, not on components, so those parameters have nowhere to
         * land — and a node emitted without them is not a work gate or a
         * transactor, it is a storage wearing the name. That is the silent
         * substitution this projection exists to refuse, so it is reported
         * and the node is dropped. */
        if (!strcmp(type, "interaction") || !strcmp(type, "gain") ||
            !strcmp(type, "loop_limited") || !strcmp(type, "exchange") ||
            !strcmp(type, "switch")) {
            const cJSON *np = cJSON_GetObjectItemCaseSensitive(it, "params");
            if (cJSON_IsObject(np) && cJSON_GetArraySize(np) > 0) {
                if (cJSON_IsString(cJSON_GetObjectItemCaseSensitive(np,
                                                            "price_node")))
                    finding(cov, "node '%s': '%s' with an endogenous price "
                                 "resolved from a node — this engine takes a "
                                 "constant exchange ratio (ADR 0001)",
                            id ? id : "?", type);
                else
                    finding(cov, "node '%s': '%s' carries its law in node "
                                 "params, and this engine puts laws on "
                                 "pathways — emitting it without them would "
                                 "make it a storage wearing the name",
                            id ? id : "?", type);
                continue;
            }
        }

        n = cJSON_CreateObject();
        if (!n) continue;
        cJSON_AddStringToObject(n, "id", id ? id : "?");
        cJSON_AddStringToObject(n, "type", type);
        /* GSSK's `value` is the MOP engine's initial quantity. Storage reads
         * current_level, everything else reads value; write both so neither
         * has to guess. */
        cJSON_AddNumberToObject(n, "value", nfield(it, "value", 0.0));
        if (!strcmp(type, "storage"))
            cJSON_AddNumberToObject(n, "current_level", nfield(it, "value", 0.0));
        {   /* Carriers exist on both sides and mean the same thing. */
            const char *c = sfield(it, "carrier", NULL);
            if (c) cJSON_AddStringToObject(n, "carrier", c);
        }
        {   /* Only self-generating waveforms can be carried as state. */
            const cJSON *f = cJSON_GetObjectItemCaseSensitive(it, "forcing");
            if (cJSON_IsObject(f)) {
                const char *k = sfield(f, "kind", "none");
                if (!strcmp(k, "sine") || !strcmp(k, "ramp") ||
                    !strcmp(k, "exponential") || !strcmp(k, "none")) {
                    cJSON_AddItemToObject(n, "forcing", cJSON_Duplicate(f, 1));
                } else {
                    finding(cov, "node '%s': forcing '%s' dropped — not its own "
                                 "generator, so it cannot be carried as state",
                            id ? id : "?", k);
                }
            }
        }
        cJSON_AddItemToArray(onodes, n);
        cov->nodes_carried++;
    }

    cJSON_ArrayForEach(it, edges) {
        const char *eid   = sfield(it, "id", NULL);
        const char *logic = sfield(it, "logic", "linear");
        const char *why   = NULL;
        cJSON      *e;

        cov->edges_total++;
        if (!logic_carried(logic, &why)) {
            finding(cov, "edge '%s': logic '%s' not carried — %s",
                    eid ? eid : "?", logic, why);
            continue;
        }
        params = cJSON_GetObjectItemCaseSensitive((cJSON *)it, "params");

        /* An n-ary work gate is a real arity, not a spelling (ADR 0008), and
         * this engine reads one control. Carrying it with the first control
         * silently would change the model's law. */
        {
            const cJSON *cn = cJSON_GetObjectItemCaseSensitive(params,
                                                               "control_nodes");
            if (cJSON_IsArray(cn) && cJSON_GetArraySize(cn) > 1) {
                finding(cov, "edge '%s': %d control nodes — this engine reads "
                             "one, and an n-ary product is a different law "
                             "(ADR 0008)", eid ? eid : "?",
                        cJSON_GetArraySize(cn));
                continue;
            }
        }
        /* An endogenous price is resolved from a node each step; this engine
         * takes a constant (ADR 0001). */
        if (cJSON_IsString(cJSON_GetObjectItemCaseSensitive(params,
                                                            "price_node"))) {
            finding(cov, "edge '%s': price resolved from a node — this engine "
                         "takes a constant exchange ratio (ADR 0001)",
                    eid ? eid : "?");
            continue;
        }

        e = cJSON_CreateObject();
        if (!e) continue;
        cJSON_AddStringToObject(e, "source", sfield(it, "origin", "?"));
        cJSON_AddStringToObject(e, "target", sfield(it, "target", "?"));
        cJSON_AddStringToObject(e, "logic",  logic);
        cJSON_AddNumberToObject(e, "weight", nfield(params, "k", 1.0));
        if (cJSON_IsNumber(cJSON_GetObjectItemCaseSensitive(params, "C")))
            cJSON_AddNumberToObject(e, "capacity", nfield(params, "C", 0.0));
        if (cJSON_IsNumber(cJSON_GetObjectItemCaseSensitive(params, "threshold")))
            cJSON_AddNumberToObject(e, "threshold",
                                    nfield(params, "threshold", 0.0));
        {
            const char *cn = sfield(params, "control_node", NULL);
            if (cn) cJSON_AddStringToObject(e, "control_node", cn);
        }
        cJSON_AddItemToArray(oedges, e);
        cov->edges_carried++;
    }

    {   /* Horizon, so the projected model is runnable rather than a fragment. */
        const cJSON *cfg = cJSON_GetObjectItemCaseSensitive(gssk, "config");
        cJSON *sp = cJSON_AddObjectToObject(out, "simulation_params");
        if (sp) {
            cJSON_AddNumberToObject(sp, "t_val", nfield(cfg, "t_end", 1.0));
            cJSON_AddNumberToObject(sp, "derivative_order", 1);
            cJSON_AddBoolToObject(sp, "generative_mode", 0);
        }
    }

    *out_mop = out;
    return true;
}

double gia_coverage_fraction(const gia_coverage *cov) {
    int total, carried;
    if (!cov) return 0.0;
    total   = cov->nodes_total  + cov->edges_total;
    carried = cov->nodes_carried + cov->edges_carried;
    if (total <= 0) return 1.0;
    return (double)carried / (double)total;
}

void gia_report_coverage(const gia_coverage *cov) {
    int i;
    if (!cov) return;

    printf("\n============================================================\n");
    printf(" MOP COVERAGE OF THIS MODEL\n");
    printf("============================================================\n");
    printf("  components     %d / %d\n", cov->nodes_carried, cov->nodes_total);
    printf("  pathways       %d / %d\n", cov->edges_carried, cov->edges_total);
    printf("  coverage       %.1f%%\n", 100.0 * gia_coverage_fraction(cov));

    if (cov->n_findings == 0) {
        printf("\n  Everything in this model has a MOP representation.\n");
    } else {
        printf("\n  What could not be carried:\n");
        for (i = 0; i < cov->n_findings; i++)
            printf("    - %s\n", cov->finding[i]);
        if (cov->n_elided)
            printf("    (%d further findings not listed)\n", cov->n_elided);
        printf("\n  A score below 100%% says THIS MODEL uses modules the\n");
        printf("  engine cannot carry. It does not say the framework is that\n");
        printf("  fraction correct, which is why the modules are named.\n");
    }
    printf("============================================================\n");
}
