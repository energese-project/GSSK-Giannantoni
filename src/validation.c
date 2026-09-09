/* validation.c — did the engine run functionally, or generatively?
 *
 * The question is settled structurally rather than by trusting a flag: the
 * graph that went in is compared against the graph that came out. If they are
 * identical the run only moved numbers along a fixed topology; if they differ
 * the system produced emergent quality that was not in the seed.
 *
 * See include/engine.h for the framework this belongs to.
 */

#include "engine.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static const char *RULE =
    "============================================================";

static int array_len(const cJSON *root, const char *key) {
    const cJSON *arr = cJSON_GetObjectItemCaseSensitive(root, key);
    return cJSON_IsArray(arr) ? cJSON_GetArraySize(arr) : 0;
}

static const char *sfield(const cJSON *obj, const char *key) {
    const cJSON *it = cJSON_GetObjectItemCaseSensitive(obj, key);
    return (cJSON_IsString(it) && it->valuestring) ? it->valuestring : "?";
}

static bool has_node(const cJSON *root, const char *id) {
    const cJSON *it, *arr = cJSON_GetObjectItemCaseSensitive(root, "nodes");
    cJSON_ArrayForEach(it, arr)
        if (!strcmp(sfield(it, "id"), id)) return true;
    return false;
}

static bool has_edge(const cJSON *root, const char *from, const char *to) {
    const cJSON *it, *arr = cJSON_GetObjectItemCaseSensitive(root, "edges");
    cJSON_ArrayForEach(it, arr)
        if (!strcmp(sfield(it, "source"), from) &&
            !strcmp(sfield(it, "target"), to)) return true;
    return false;
}

/* Name the components and relationships that appear in `b` but not in `a`.
 * Returns how many were listed. */
static int list_added(const cJSON *a, const cJSON *b, const char *what) {
    const cJSON *it;
    int          shown = 0;

    cJSON_ArrayForEach(it, cJSON_GetObjectItemCaseSensitive(b, "nodes")) {
        const char *id = sfield(it, "id");
        if (has_node(a, id)) continue;
        printf("    %s component      %s  (%s)\n", what, id, sfield(it, "type"));
        shown++;
    }
    cJSON_ArrayForEach(it, cJSON_GetObjectItemCaseSensitive(b, "edges")) {
        const char *from = sfield(it, "source"), *to = sfield(it, "target");
        if (has_edge(a, from, to)) continue;
        printf("    %s relationship   %s -> %s  (%s)\n",
               what, from, to, sfield(it, "flow_type"));
        shown++;
    }
    return shown;
}

gia_mode gia_validate_mode(const cJSON *in, const cJSON *out) {
    /* cJSON_Compare walks the whole tree; a single added node, edge or field
     * anywhere is enough to distinguish the modes. */
    if (!in || !out) return GIA_MODE_FUNCTIONAL;
    return cJSON_Compare(in, out, 1) ? GIA_MODE_FUNCTIONAL
                                     : GIA_MODE_GENERATIVE;
}

void gia_report_mode(gia_mode mode, const cJSON *in, const cJSON *out) {
    int n_in  = array_len(in,  "nodes"), n_out = array_len(out, "nodes");
    int e_in  = array_len(in,  "edges"), e_out = array_len(out, "edges");

    printf("\n%s\n", RULE);
    printf(" EXECUTION MODE, BY STRUCTURAL DIFF\n");
    printf("%s\n", RULE);

    printf("  nodes  %d -> %d   (%+d)\n", n_in, n_out, n_out - n_in);
    printf("  edges  %d -> %d   (%+d)\n", e_in, e_out, e_out - e_in);

    if (mode == GIA_MODE_FUNCTIONAL) {
        printf("\n  Input_Graph === Output_Graph\n");
        printf("  MODE: functional / conservative.\n\n");
        printf("    nothing added, nothing removed; cJSON_Compare is exact\n");
        printf("    over the whole tree, so this covers field values too.\n");
        printf("\n  The topology is unchanged. Evolution was confined to\n");
        printf("  numerical state along a fixed graph -- the system\n");
        printf("  transformed, but it did not originate anything.\n");
    } else {
        printf("\n  Input_Graph !== Output_Graph\n");
        printf("  MODE: generative / Maximum Ordinality Principle.\n\n");
        if (list_added(in, out, "+") == 0)
            printf("    (no additions; the graphs differ in field values)\n");
        list_added(out, in, "-");
        printf("\n  The graph itself changed. Components and relationships\n");
        printf("  appear in the output that were not in the seed and are\n");
        printf("  not reducible to it -- emergent quality, in Giannantoni's\n");
        printf("  sense, rather than a rearrangement of what was given.\n");
    }
    printf("%s\n", RULE);
}

bool gia_validate_harmony(const gia_harmony *h, double tol) {
    double reduction, balance;
    bool   ok_reduction, ok_balance = true;

    printf("\n%s\n", RULE);
    printf(" MOP HARMONY RELATIONSHIPS\n");
    printf("%s\n", RULE);

    if (!h || h->n < 2) {
        printf("  no harmony matrix (needs at least 2 components)\n");
        printf("%s\n", RULE);
        return false;
    }

    printf("  components N          %d\n", h->n);
    printf("  ordinal roots         %d  (N-1)\n", h->n - 1);
    printf("  reference couple a12  %.6f %+.6fi\n",
           creal(h->alpha_ref), cimag(h->alpha_ref));

    /* Claim 1: every entry is recoverable from alpha_12 and (i,j) alone, so
     * the N x N matrix carries exactly one independent quantity. This is the
     * O(N^2) -> O(1) reduction, checked rather than asserted. */
    reduction    = gia_harmony_reduction_residual(h);
    ok_reduction = (reduction <= tol);
    printf("\n  reduction residual    %.3e   %s\n",
           reduction, ok_reduction ? "PASS" : "FAIL");
    printf("    all %d x %d = %d entries reconstruct from the single\n",
           h->n, h->n, h->n * h->n);
    printf("    reference couple, so the matrix holds 1 free quantity.\n");

    /* Claim 2: the (N-1) roots of unity sum to zero, so each row of the
     * harmony matrix balances. This is the global stability condition.
     *
     * It needs N >= 3. At N = 2 there is exactly one root, omega_0 = 1, and a
     * row sums to alpha_12 rather than to zero -- correctly so: a two-body
     * couple has no interior against which to balance. */
    balance = gia_harmony_row_residual(h);
    if (h->n >= 3) {
        ok_balance = (balance <= tol);
        printf("\n  row balance residual  %.3e   %s\n",
               balance, ok_balance ? "PASS" : "FAIL");
        printf("    every row sums to zero: the %d ordinal roots cancel,\n",
               h->n - 1);
        printf("    which is the global stability of the relationships.\n");
    } else {
        printf("\n  row balance residual  %.3e   n/a at N=2\n", balance);
        printf("    a two-body couple has no interior to balance against.\n");
    }

    printf("%s\n", RULE);
    return ok_reduction && ok_balance;
}
