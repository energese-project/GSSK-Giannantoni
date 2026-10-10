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

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void added(gia_coverage *cov, const char *fmt, ...) {
    va_list ap;
    if (!cov || cov->n_added >= GIA_MAX_FINDINGS) return;
    va_start(ap, fmt);
    vsnprintf(cov->added[cov->n_added], sizeof(cov->added[0]), fmt, ap);
    va_end(ap);
    cov->n_added++;
}

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
    static const char *const prim[] = { "storage", "source", "sink", "constant",
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

/* GSSK's edge laws. The names match the kernel's deliberately, so the two
 * vocabularies cannot drift. Odum SecIII stays on the pathway; the rest are
 * module laws here (ADR 0012 decision 5) and are carried by writing the gate
 * GSSK draws implicitly -- in the new document, and declared as added. */
static bool logic_is_pathway(const char *l) {
    return !strcmp(l, "linear") || !strcmp(l, "reversible") ||
           !strcmp(l, "constant");
}

static bool logic_is_module(const char *l) {
    return !strcmp(l, "interaction") || !strcmp(l, "ratio") ||
           !strcmp(l, "subtract")    || !strcmp(l, "limit") ||
           !strcmp(l, "threshold");
}

static bool is_carried(const cJSON *onodes, const char *id) {
    const cJSON *it;
    if (!id) return false;
    cJSON_ArrayForEach(it, onodes) {
        const cJSON *v = cJSON_GetObjectItemCaseSensitive(it, "id");
        if (cJSON_IsString(v) && !strcmp(v->valuestring, id)) return true;
    }
    return false;
}

static void add_leg(cJSON *oedges, const char *from, const char *to,
                    const char *role) {
    cJSON *e = cJSON_CreateObject();
    if (!e) return;
    cJSON_AddStringToObject(e, "source", from);
    cJSON_AddStringToObject(e, "target", to);
    if (role) cJSON_AddStringToObject(e, "role", role);
    /* GSSK reads a control and never draws it, so the faithful translation
     * is use_ratio 0 -- which the MOP run report then names (ADR 0017). */
    if (role && !strcmp(role, "control"))
        cJSON_AddNumberToObject(e, "use_ratio", 0.0);
    cJSON_AddItemToArray(oedges, e);
}

/* The kernel's forcing vocabulary (`waveform`, `mean`, `period`, `phase` in
 * time, `t_on`, `v0`, `slope`, `min`/`max`) translated into this engine's
 * (`kind`, `offset`, `amplitude`, `rate`, `phase` in radians), on this
 * engine's clock s = t - t0, where [t0, t1] is the kernel's run. Only
 * waveforms that generate themselves can be carried as state (sine, ramp,
 * exponential), and only where the kernel's waveform is that one formula over
 * the whole run: a ramp or exponential that switches on after t0 is a step,
 * and a min/max clamp that can bind is a kink. Anything else returns false
 * with the reason in `why`, and the caller reports it and does not count the
 * element as carried -- a forcing silently read as a constant was the failure
 * this replaces. *out is NULL when there is no forcing. */
static bool project_forcing(const cJSON *f, double t0, double t1, cJSON **out,
                            char *why, size_t cap) {
    const char *w;
    double      lo, hi, t_on;
    bool        has_min, has_max;
    double      vmin, vmax;
    cJSON      *o;

    *out = NULL;
    if (!cJSON_IsObject(f)) return true;
    w       = sfield(f, "waveform", "none");
    t_on    = nfield(f, "t_on", t0);
    has_min = cJSON_IsNumber(cJSON_GetObjectItemCaseSensitive(f, "min"));
    has_max = cJSON_IsNumber(cJSON_GetObjectItemCaseSensitive(f, "max"));
    vmin    = nfield(f, "min", 0.0);
    vmax    = nfield(f, "max", 0.0);
    if (!(o = cJSON_CreateObject())) { snprintf(why, cap, "out of memory"); return false; }

    if (!strcmp(w, "sine")) {
        const double mean = nfield(f, "mean", 0.0), amp = nfield(f, "amplitude", 0.0);
        const double per  = nfield(f, "period", 0.0), ph = nfield(f, "phase", 0.0);
        if (!(per > 0.0)) { snprintf(why, cap, "a sine with no period"); cJSON_Delete(o); return false; }
        lo = mean - fabs(amp); hi = mean + fabs(amp);
        /* mean + amp sin(2 pi (t - t_on - phase)/period), t = s + t0. */
        cJSON_AddStringToObject(o, "kind", "sine");
        cJSON_AddNumberToObject(o, "offset", mean);
        cJSON_AddNumberToObject(o, "amplitude", amp);
        cJSON_AddNumberToObject(o, "rate", 2.0 * M_PI / per);
        cJSON_AddNumberToObject(o, "phase", 2.0 * M_PI * (t0 - t_on - ph) / per);
    } else if (!strcmp(w, "ramp") || !strcmp(w, "exponential")) {
        const double v0 = nfield(f, "v0", 0.0);
        if (t_on > t0) {
            snprintf(why, cap, "a %s that switches on at t = %g, after the run starts: a step, "
                               "not its own generator", w, t_on);
            cJSON_Delete(o);
            return false;
        }
        if (!strcmp(w, "ramp")) {
            const double slope = nfield(f, "slope", 0.0);
            const double a = v0 + slope * (t0 - t_on), b = v0 + slope * (t1 - t_on);
            lo = a < b ? a : b; hi = a < b ? b : a;
            cJSON_AddStringToObject(o, "kind", "ramp");
            cJSON_AddNumberToObject(o, "offset", a);
            cJSON_AddNumberToObject(o, "rate", slope);
        } else {
            const double r = nfield(f, "rate", 0.0), amp = v0 * exp(r * (t0 - t_on));
            const double b = amp * exp(r * (t1 - t0));
            lo = amp < b ? amp : b; hi = amp < b ? b : amp;
            cJSON_AddStringToObject(o, "kind", "exponential");
            cJSON_AddNumberToObject(o, "offset", 0.0);
            cJSON_AddNumberToObject(o, "amplitude", amp);
            cJSON_AddNumberToObject(o, "rate", r);
        }
    } else {
        snprintf(why, cap, "waveform '%s' is not its own generator, so it cannot be carried as "
                           "state", w);
        cJSON_Delete(o);
        return false;
    }
    if ((has_min && lo < vmin) || (has_max && hi > vmax)) {
        snprintf(why, cap, "its %s clamp binds during the run, and a clamped waveform is not its "
                           "own generator", (has_min && lo < vmin) ? "min" : "max");
        cJSON_Delete(o);
        return false;
    }
    *out = o;
    return true;
}

bool gia_project(const cJSON *gssk, cJSON **out_mop, gia_coverage *cov) {
    const cJSON *nodes, *edges, *it;
    cJSON       *out, *onodes, *oedges, *params;
    double       t0, t1;

    if (!gssk || !out_mop || !cov) return false;
    memset(cov, 0, sizeof(*cov));
    *out_mop = NULL;

    nodes = cJSON_GetObjectItemCaseSensitive(gssk, "nodes");
    edges = cJSON_GetObjectItemCaseSensitive(gssk, "edges");
    if (!cJSON_IsArray(nodes)) return false;
    {   /* The kernel's run, [t0, t1]; this engine's clock is s = t - t0. */
        const cJSON *cfg = cJSON_GetObjectItemCaseSensitive(gssk, "config");
        t0 = nfield(cfg, "t_start", 0.0);
        t1 = nfield(cfg, "t_end", t0 + 1.0);
    }

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
         * (Phase 7: k, C, threshold, price) and finds its legs by POSITION.
         * Here a module names its legs by role (ADR 0013), and guessing which
         * pathway is the energy and which the control is the silent
         * substitution this projection exists to refuse. With no params it
         * fares no better: a module type without its law is a storage wearing
         * the name (ADR 0012). Either way it is reported and dropped. */
        if (!strcmp(type, "interaction") || !strcmp(type, "gain") ||
            !strcmp(type, "loop_limited") || !strcmp(type, "exchange") ||
            !strcmp(type, "switch")) {
            const cJSON *np = cJSON_GetObjectItemCaseSensitive(it, "params");
            if (cJSON_IsObject(np) &&
                cJSON_IsString(cJSON_GetObjectItemCaseSensitive(np, "price_node")))
                finding(cov, "node '%s': '%s' with an endogenous price "
                             "resolved from a node — this engine takes a "
                             "constant exchange ratio (ADR 0001)",
                        id ? id : "?", type);
            else if (cJSON_IsObject(np) && cJSON_GetArraySize(np) > 0)
                finding(cov, "node '%s': '%s' carries its law in node params "
                             "and finds its legs by position; this engine names "
                             "a module's legs by role (ADR 0013)",
                        id ? id : "?", type);
            else
                finding(cov, "node '%s': '%s' with no law of its own would be "
                             "a storage wearing the name (ADR 0012)",
                        id ? id : "?", type);
            continue;
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
        {   /* Only self-generating waveforms can be carried as state. A node
             * whose forcing cannot be is still written (its edges need it),
             * held at its value, and NOT counted as carried: its dynamics
             * are not the GSSK model's. */
            cJSON *pf = NULL;
            char   why_f[160];
            if (!project_forcing(cJSON_GetObjectItemCaseSensitive(it, "forcing"), t0, t1, &pf,
                                 why_f, sizeof why_f)) {
                finding(cov, "node '%s': forcing not carried — %s", id ? id : "?", why_f);
                cJSON_AddItemToArray(onodes, n);
                continue;
            }
            if (pf) cJSON_AddItemToObject(n, "forcing", pf);
        }
        cJSON_AddItemToArray(onodes, n);
        cov->nodes_carried++;
    }

    cJSON_ArrayForEach(it, edges) {
        const char   *eid    = sfield(it, "id", NULL);
        const char   *logic  = sfield(it, "logic", "linear");
        const char   *origin = sfield(it, "origin", NULL);
        const char   *target = sfield(it, "target", NULL);
        const cJSON  *cns;
        const char   *ctrl;
        int           n_ctrl;
        double        k;
        cJSON        *e;

        cov->edges_total++;
        if (!logic_is_pathway(logic) && !logic_is_module(logic)) {
            finding(cov, "edge '%s': logic '%s' not carried — not a law this "
                         "engine implements", eid ? eid : "?", logic);
            continue;
        }
        params = cJSON_GetObjectItemCaseSensitive((cJSON *)it, "params");
        k      = nfield(params, "k", 1.0);

        /* An endogenous price is resolved from a node each step; this engine
         * takes a constant (ADR 0001). */
        if (cJSON_IsString(cJSON_GetObjectItemCaseSensitive(params,
                                                            "price_node"))) {
            finding(cov, "edge '%s': price resolved from a node — this engine "
                         "takes a constant exchange ratio (ADR 0001)",
                    eid ? eid : "?");
            continue;
        }

        /* Controls: one named, several (ADR 0008), or by GSSK's default the
         * target -- the autocatalytic reading. */
        cns    = cJSON_GetObjectItemCaseSensitive(params, "control_nodes");
        ctrl   = sfield(params, "control_node", target);
        n_ctrl = cJSON_IsArray(cns) ? cJSON_GetArraySize(cns) : 1;
        if (n_ctrl > 1 && (!strcmp(logic, "ratio") || !strcmp(logic, "subtract"))) {
            finding(cov, "edge '%s': %d control nodes on '%s' — a quotient or "
                         "difference of more than two is not in Odum Fig. 2.6 "
                         "(ADR 0016)", eid ? eid : "?", n_ctrl, logic);
            continue;
        }

        /* Every endpoint must have been carried. An edge to a dropped node
         * would be ignored on load, and a gate missing a leg would not load at
         * all, so it is dropped here and named. */
        {
            bool ends = is_carried(onodes, origin) && is_carried(onodes, target);
            if (logic_is_module(logic) && strcmp(logic, "limit") &&
                strcmp(logic, "threshold")) {
                if (cJSON_IsArray(cns)) {
                    const cJSON *c;
                    cJSON_ArrayForEach(c, cns)
                        if (!cJSON_IsString(c) || !is_carried(onodes, c->valuestring))
                            ends = false;
                } else if (!is_carried(onodes, ctrl)) {
                    ends = false;
                }
            }
            if (!ends) {
                finding(cov, "edge '%s': it touches a node that was not "
                             "carried", eid ? eid : "?");
                continue;
            }
        }

        {   /* An edge's forcing drives its rate k (ADR 0006). */
            const cJSON *ef = cJSON_GetObjectItemCaseSensitive(it, "forcing");
            cJSON       *pf = NULL;
            char         why_f[160];
            if (cJSON_IsObject(ef) && !logic_is_pathway(logic)) {
                finding(cov, "edge '%s': forcing on a '%s' law not carried — this engine "
                             "forces a pathway's rate, not a module's", eid ? eid : "?", logic);
                continue;
            }
            if (!project_forcing(ef, t0, t1, &pf, why_f, sizeof why_f)) {
                finding(cov, "edge '%s': forcing not carried — %s", eid ? eid : "?", why_f);
                continue;
            }
            if (logic_is_pathway(logic)) {
                e = cJSON_CreateObject();
                if (!e) { cJSON_Delete(pf); continue; }
                cJSON_AddStringToObject(e, "source", origin);
                cJSON_AddStringToObject(e, "target", target);
                cJSON_AddStringToObject(e, "logic",  logic);
                cJSON_AddNumberToObject(e, "weight", k);
                if (pf) cJSON_AddItemToObject(e, "forcing", pf);
                cJSON_AddItemToArray(oedges, e);
                cov->edges_carried++;
                continue;
            }
        }

        {   /* A module law: write the gate GSSK draws implicitly on the edge,
             * with its legs named. */
            char        gid[160];
            const char *type, *action = NULL;
            cJSON      *gate, *mod;

            if (eid) snprintf(gid, sizeof(gid), "%s__gate", eid);
            else     snprintf(gid, sizeof(gid), "%s_%s__gate", origin, target);

            if      (!strcmp(logic, "limit"))     type = "loop_limited";
            else if (!strcmp(logic, "threshold")) type = "switch";
            else {
                type = "interaction";
                if (!strcmp(logic, "ratio"))    action = "divide";
                if (!strcmp(logic, "subtract")) action = "subtract";
            }

            gate = cJSON_CreateObject();
            mod  = cJSON_CreateObject();
            if (!gate || !mod) { cJSON_Delete(gate); cJSON_Delete(mod); continue; }
            cJSON_AddStringToObject(gate, "id", gid);
            cJSON_AddStringToObject(gate, "type", type);
            cJSON_AddNumberToObject(mod, "k", k);
            if (action) cJSON_AddStringToObject(mod, "action", action);
            if (!strcmp(type, "loop_limited"))
                cJSON_AddNumberToObject(mod, "capacity", nfield(params, "C", 1.0));
            if (!strcmp(type, "switch"))
                cJSON_AddNumberToObject(mod, "threshold",
                                        nfield(params, "threshold", 0.0));
            cJSON_AddItemToObject(gate, "module", mod);
            cJSON_AddItemToArray(onodes, gate);

            add_leg(oedges, origin, gid, "energy");
            if (!strcmp(type, "interaction")) {
                if (cJSON_IsArray(cns)) {
                    const cJSON *c;
                    cJSON_ArrayForEach(c, cns) add_leg(oedges, c->valuestring, gid, "control");
                } else {
                    add_leg(oedges, ctrl, gid, "control");
                }
            }
            /* GSSK's threshold reads the level of the stock it drains, so the
             * switch's sensor is its own energy input. */
            if (!strcmp(type, "switch")) add_leg(oedges, origin, gid, "control");
            add_leg(oedges, gid, target, NULL);

            added(cov, "edge '%s': GSSK '%s' on the edge becomes %s module "
                       "'%s', a node the GSSK author did not write (ADR 0012)",
                  eid ? eid : "?", logic, type, gid);
            cov->edges_carried++;
        }
    }

    {   /* Horizon, so the projected model is runnable rather than a fragment. */
        cJSON *sp = cJSON_AddObjectToObject(out, "simulation_params");
        if (sp) {
            /* This engine's clock starts at 0: the kernel's t_start is s = 0. */
            cJSON_AddNumberToObject(sp, "t_val", t1 - t0);
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

    if (cov->n_added > 0) {
        /* Carried, not lost -- but these nodes move the component and module
         * counts, so a reader comparing ordinality needs to see them. */
        printf("\n  Added in translation:\n");
        for (i = 0; i < cov->n_added; i++)
            printf("    + %s\n", cov->added[i]);
    }

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
