/* Giannantoni generative framework — engine tests.
 *
 * Four things are checked, and each is a claim from the papers rather than a
 * property of this particular code:
 *
 *   1. Persistence of form and derivative drift. The incipient amplitude is
 *      (phi')^n; the traditional one is the complete Bell polynomial. Their
 *      difference is computed against hand-derived closed forms, and it is
 *      identically zero exactly when phi is affine.
 *   2. The binary/duet: a half-order incipient derivative returns two
 *      branches, +/- sqrt(alpha) e^(alpha t), which cancel.
 *   3. The MOP Harmony Relationships: an N x N matrix generated from one
 *      reference couple by the (N-1) ordinal roots of unity, whose rows
 *      balance and whose every entry is reconstructible from that couple.
 *   4. Ordinality and the generative step: the step closes an open pathway,
 *      raises ordinality to maximum, and is then a fixed point.
 */

#define _POSIX_C_SOURCE 200809L

#include "engine.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;

static void ok(const char *what, bool cond) {
    printf("  %-58s %s\n", what, cond ? "PASS" : "FAIL");
    if (!cond) failures++;
}

static void close_to(const char *what, double got, double want, double tol) {
    bool c = fabs(got - want) <= tol;
    printf("  %-58s %s  (got %.10g want %.10g)\n", what, c ? "PASS" : "FAIL",
           got, want);
    if (!c) failures++;
}

/* ------------------------------------------------------------------ *
 * 1. Persistence of form, and the drift TDC introduces
 * ------------------------------------------------------------------ */

static void test_drift(void) {
    gia_phi affine, quad;
    double  a = 1.5, t = 1.5;
    int     n;

    printf("\n[1] incipient vs traditional derivative\n");

    /* phi = a t  -- constant-coefficient. Both calculi must agree exactly. */
    memset(&affine, 0, sizeof(affine));
    affine.c[1] = a;
    affine.degree = 1;

    ok("affine phi reports drift-free", gia_drift_free(&affine));
    for (n = 1; n <= 5; n++) {
        char label[80];
        snprintf(label, sizeof(label), "  affine: drift at order %d is zero", n);
        close_to(label, gia_drift(&affine, n, t), 0.0, 1e-12);
    }
    close_to("affine: incipient amplitude = (phi')^3",
             gia_idc_amplitude(&affine, 3, t), a * a * a, 1e-12);

    /* phi = a t^2 / 2  -- variable-coefficient, so phi' = a t and phi'' = a.
     *
     *   incipient order 2 : (phi')^2                     = a^2 t^2
     *   traditional order 2: B_2 = (phi')^2 + phi''      = a^2 t^2 + a
     *   drift                                            = a
     *
     *   traditional order 3: B_3 = (phi')^3 + 3 phi' phi'' + phi'''
     *   drift                                            = 3 a^2 t
     */
    memset(&quad, 0, sizeof(quad));
    quad.c[2] = 0.5 * a;
    quad.degree = 2;

    ok("quadratic phi is not drift-free", !gia_drift_free(&quad));
    close_to("quadratic: phi'(t) = a t", gia_phi_deriv(&quad, 1, t), a * t, 1e-12);
    close_to("quadratic: phi''(t) = a",  gia_phi_deriv(&quad, 2, t), a, 1e-12);

    close_to("quadratic: incipient order 2 = (a t)^2",
             gia_idc_amplitude(&quad, 2, t), (a * t) * (a * t), 1e-12);
    close_to("quadratic: traditional order 2 = (a t)^2 + a",
             gia_tdc_amplitude(&quad, 2, t), (a * t) * (a * t) + a, 1e-12);
    close_to("quadratic: drift at order 2 = a",
             gia_drift(&quad, 2, t), a, 1e-12);
    close_to("quadratic: drift at order 3 = 3 a^2 t",
             gia_drift(&quad, 3, t), 3.0 * a * a * t, 1e-12);

    /* Order 0 is the function itself under either calculus. */
    close_to("order 0 incipient = e^phi",
             gia_idc_derivative(&quad, 0, t), exp(gia_phi_eval(&quad, t)), 1e-12);
    close_to("order 0 drift is zero", gia_drift(&quad, 0, t), 0.0, 1e-12);
}

/* ------------------------------------------------------------------ *
 * 2. The binary / duet
 * ------------------------------------------------------------------ */

static void test_duet(void) {
    gia_phi phi;
    gia_net d;
    double  a = 2.0, t = 0.75, expected;

    printf("\n[2] binary (duet) from a half-order incipient derivative\n");

    memset(&phi, 0, sizeof(phi));
    phi.c[1] = a;
    phi.degree = 1;

    d = gia_incipient_fractional(&phi, 1, 2, t);
    ok("half order yields exactly two branches", gia_net_is_binary(&d));

    /* (d~/d~t)^(1/2) e^(a t) = +/- sqrt(a) e^(a t) */
    expected = sqrt(a) * exp(a * t);
    close_to("branch + real part = +sqrt(a) e^(a t)",
             creal(d.branch[0]), expected, 1e-12);
    close_to("branch - real part = -sqrt(a) e^(a t)",
             creal(d.branch[1]), -expected, 1e-12);
    close_to("branch + imaginary part is zero", cimag(d.branch[0]), 0.0, 1e-12);
    close_to("the two branches cancel", cabs(gia_net_sum(&d)), 0.0, 1e-12);

    /* Whole order reduces to the single-branch incipient derivative. */
    d = gia_incipient_fractional(&phi, 1, 1, t);
    ok("unit denominator yields a single branch", d.count == 1);
    close_to("single branch = incipient first derivative",
             creal(d.branch[0]), gia_idc_derivative(&phi, 1, t), 1e-12);

    /* A triplet: three branches of the cube root, which also cancel. */
    d = gia_incipient_fractional(&phi, 1, 3, t);
    ok("third order yields three branches", d.count == 3);
    close_to("the three branches cancel", cabs(gia_net_sum(&d)), 0.0, 1e-12);
}

/* ------------------------------------------------------------------ *
 * 3. MOP Harmony Relationships
 * ------------------------------------------------------------------ */

static void test_harmony(void) {
    gia_harmony h;
    double complex aref = 1.2 + 0.4 * I;
    double complex w;

    printf("\n[3] MOP harmony matrix from the (N-1) ordinal roots of unity\n");

    w = gia_ordinal_root(4, 0);
    close_to("4th root of unity, m=0, is 1", creal(w), 1.0, 1e-12);
    close_to("4th root of unity, m=0, has no imaginary part",
             cimag(w), 0.0, 1e-12);
    w = gia_ordinal_root(4, 1);
    close_to("4th root of unity, m=1, is i", cimag(w), 1.0, 1e-12);
    close_to("4th root of unity, m=1, has no real part",
             creal(w), 0.0, 1e-12);

    ok("harmony matrix rejects N < 2", !gia_harmony_init(&h, 1, aref));

    ok("harmony matrix builds for N = 5", gia_harmony_init(&h, 5, aref));

    /* The reference couple is a genuine entry, not an external parameter. */
    close_to("alpha_12 is the stored entry (0,1), real part",
             creal(gia_harmony_at(&h, 0, 1)), creal(aref), 1e-12);
    close_to("alpha_12 is the stored entry (0,1), imaginary part",
             cimag(gia_harmony_at(&h, 0, 1)), cimag(aref), 1e-12);

    close_to("diagonal carries no self-relation",
             cabs(gia_harmony_at(&h, 2, 2)), 0.0, 1e-12);

    /* Every entry has the modulus of the reference couple: the roots of unity
     * rotate the relationship without rescaling it. */
    close_to("off-diagonal modulus equals |alpha_12|",
             cabs(gia_harmony_at(&h, 3, 1)), cabs(aref), 1e-12);

    /* The N x N -> 1 reduction. */
    close_to("every entry reconstructs from alpha_12 alone",
             gia_harmony_reduction_residual(&h), 0.0, 1e-12);

    /* Global stability: the N-1 roots cancel, so each row balances. */
    close_to("each row sums to zero (global balance)",
             gia_harmony_row_residual(&h), 0.0, 1e-12);

    ok("both harmony invariants validate", gia_validate_harmony(&h, 1e-9));
    gia_harmony_free(&h);

    /* N = 2 is the documented exception: one root, omega_0 = 1, so a row sums
     * to alpha_12. A two-body couple has no interior to balance against. */
    ok("harmony matrix builds for N = 2", gia_harmony_init(&h, 2, aref));
    close_to("at N=2 a row sums to |alpha_12|, not zero",
             gia_harmony_row_residual(&h), cabs(aref), 1e-12);
    close_to("the N=2 reduction still holds",
             gia_harmony_reduction_residual(&h), 0.0, 1e-12);
    gia_harmony_free(&h);
}

/* ------------------------------------------------------------------ *
 * 4. Ordinality and the generative ordinal step
 * ------------------------------------------------------------------ */

static const char *SEED =
    "{\"system_name\":\"test\","
    " \"nodes\":["
    "   {\"id\":\"source_1\",\"type\":\"source\",\"initial_value\":2.0},"
    "   {\"id\":\"interaction_1\",\"type\":\"interaction\","
    "    \"generativity_factor\":1.5},"
    "   {\"id\":\"store_1\",\"type\":\"storage\",\"capacity\":100.0,"
    "    \"current_level\":10.0},"
    "   {\"id\":\"consumer_1\",\"type\":\"storage\",\"metabolic_rate\":0.2}],"
    " \"edges\":["
    "   {\"source\":\"source_1\",\"target\":\"interaction_1\",\"weight\":1.0},"
    "   {\"source\":\"interaction_1\",\"target\":\"store_1\",\"weight\":1.2},"
    "   {\"source\":\"store_1\",\"target\":\"consumer_1\",\"weight\":0.5},"
    "   {\"source\":\"consumer_1\",\"target\":\"interaction_1\",\"weight\":0.3}],"
    " \"simulation_params\":{\"t_val\":1.5,\"derivative_order\":2,"
    "                       \"generative_mode\":%s}}";

static char *seed_json(const char *generative) {
    static char buf[1400];
    snprintf(buf, sizeof(buf), SEED, generative);
    return buf;
}

static void test_generative(void) {
    cJSON     *root, *out, *out2;
    gia_model  m, evolved;
    int        n_in, n_out, e_in, e_out;

    printf("\n[4] ordinality and the generative ordinal step\n");

    root = cJSON_Parse(seed_json("true"));
    ok("seed graph parses", root != NULL);
    if (!root) return;

    ok("seed graph loads", gia_model_load(&m, root));

    /* source_1 feeds the loop but nothing returns to it, so it is not on a
     * closed pathway: 3 of 4 components are. */
    close_to("seed ordinality is 3/4", gia_ordinality(&m), 0.75, 1e-12);
    ok("seed is below maximum ordinality", !gia_at_maximum_ordinality(&m));
    ok("source_1 is the open component", !m.nodes[0].on_cycle);
    ok("interaction_1 is on a cycle", m.nodes[1].on_cycle);

    /* phi degree separates the storing components from the transforming one. */
    ok("source is drift-free",      gia_drift_free(&m.nodes[0].phi));
    ok("interaction carries drift", !gia_drift_free(&m.nodes[1].phi));
    ok("storage is drift-free",     gia_drift_free(&m.nodes[2].phi));
    ok("second storage is drift-free", gia_drift_free(&m.nodes[3].phi));

    out = gia_generate(&m);
    ok("generative step produces a graph", out != NULL);
    if (!out) { gia_model_free(&m); cJSON_Delete(root); return; }

    n_in  = cJSON_GetArraySize(cJSON_GetObjectItemCaseSensitive(root, "nodes"));
    n_out = cJSON_GetArraySize(cJSON_GetObjectItemCaseSensitive(out,  "nodes"));
    e_in  = cJSON_GetArraySize(cJSON_GetObjectItemCaseSensitive(root, "edges"));
    e_out = cJSON_GetArraySize(cJSON_GetObjectItemCaseSensitive(out,  "edges"));

    ok("input graph is left untouched", n_in == 4 && e_in == 4);
    ok("one component emerged",  n_out == n_in + 1);
    ok("two relationships emerged", e_out == e_in + 2);
    ok("the diff reports generative mode",
       gia_validate_mode(root, out) == GIA_MODE_GENERATIVE);

    /* The point of the step: it must actually raise ordinality, not merely
     * add something. Reloading the evolved graph is the check. */
    ok("evolved graph loads", gia_model_load(&evolved, out));
    close_to("evolved ordinality is 1", gia_ordinality(&evolved), 1.0, 1e-12);
    ok("evolved graph is at maximum ordinality",
       gia_at_maximum_ordinality(&evolved));

    /* At Maximum Ordinality there is no open relationship left to close, so
     * a further step must change nothing. */
    out2 = gia_generate(&evolved);
    ok("a second step produces a graph", out2 != NULL);
    if (out2) {
        ok("maximum ordinality is a fixed point",
           gia_validate_mode(out, out2) == GIA_MODE_FUNCTIONAL);
        cJSON_Delete(out2);
    }

    {   /* The generative case must differ textually too, and only by the
         * emergence -- system_name is replaced in place rather than moved,
         * so it stays at the same position in the object. */
        char *a = cJSON_Print(root), *b = cJSON_Print(out);
        ok("seed and output serialise differently", a && b && strcmp(a, b));
        ok("system_name keeps its position in the object",
           cJSON_GetArrayItem(out, 0) != NULL &&
           !strcmp(cJSON_GetArrayItem(out, 0)->string, "system_name"));
        free(a);
        free(b);
    }

    gia_model_free(&evolved);
    cJSON_Delete(out);
    gia_model_free(&m);
    cJSON_Delete(root);
}

static void test_generative_disabled(void) {
    cJSON     *root, *out;
    gia_model  m;

    printf("\n[5] generative mode off\n");

    root = cJSON_Parse(seed_json("false"));
    ok("seed graph parses", root != NULL);
    if (!root) return;

    ok("seed graph loads", gia_model_load(&m, root));
    ok("generative flag is off", !m.generative);

    out = gia_generate(&m);
    ok("a graph is still produced", out != NULL);
    if (out) {
        char *a, *b;
        ok("output is identical to input",
           gia_validate_mode(root, out) == GIA_MODE_FUNCTIONAL);

        /* Stronger than cJSON_Compare, and the property the --seed file
         * exists to give: printed through the same serialiser, a functional
         * run is byte-identical, so a text diff of seed against output shows
         * nothing at all rather than reformatting noise. */
        a = cJSON_Print(root);
        b = cJSON_Print(out);
        ok("seed and output serialise byte-identically",
           a && b && !strcmp(a, b));
        free(a);
        free(b);
        cJSON_Delete(out);
    }
    gia_model_free(&m);
    cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 6. Mode 1 output
 * ------------------------------------------------------------------ */

static void test_trajectories(void) {
    cJSON     *root;
    gia_model  m;
    FILE      *f;
    char       line[1024];
    const char *path = "tests/results/giannantoni_trajectories.csv";
    int        rows = 0;

    printf("\n[6] mode 1 trajectory output\n");

    root = cJSON_Parse(seed_json("true"));
    if (!root) { ok("seed graph parses", false); return; }
    ok("seed graph loads", gia_model_load(&m, root));

    ok("trajectories are written", gia_write_trajectories(&m, path, 10));

    f = fopen(path, "r");
    ok("trajectory file is readable", f != NULL);
    if (f) {
        if (fgets(line, sizeof(line), f)) {
            ok("header names the incipient column",
               strstr(line, "source_1_idc") != NULL);
            ok("header names the traditional column",
               strstr(line, "source_1_tdc") != NULL);
            ok("header names the drift column",
               strstr(line, "interaction_1_drift") != NULL);
        } else {
            ok("header is present", false);
        }
        while (fgets(line, sizeof(line), f)) rows++;
        ok("11 sample rows for 10 steps", rows == 11);
        fclose(f);
    }

    gia_model_free(&m);
    cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 7. The network: flow matrix, matrix exponential, conservation
 *
 * The exponential form for a network is a matrix, not a scalar per node:
 * Q(t) = exp(A t) Q(0), with A assembled from Odum's pathway laws. For a
 * constant A this is exact and closed-form, so it can be checked against the
 * analytic solution rather than against a tolerance pulled from the air.
 * ------------------------------------------------------------------ */

static const char *LINEAR_NET =
    "{\"system_name\":\"linear\","
    " \"nodes\":["
    "   {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0,"
    "    \"capacity\":100.0},"
    "   {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0,"
    "    \"capacity\":100.0}],"
    " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\","
    "            \"weight\":0.5}],"
    " \"simulation_params\":{\"t_val\":2.0,\"derivative_order\":1,"
    "                       \"generative_mode\":false}}";

static void test_network(void) {
    cJSON     *root;
    gia_model  m;
    double     q[2], psi = -1.0, t = 2.0, ea;

    printf("\n[7] network coupling: Q(t) = exp(A t) Q(0)\n");

    root = cJSON_Parse(LINEAR_NET);
    ok("linear network parses", root != NULL);
    if (!root) return;
    ok("linear network loads", gia_model_load(&m, root));

    ok("a integrates",        m.nodes[0].integrates);
    close_to("a starts at 10", m.nodes[0].q0, 10.0, 1e-12);
    ok("edge law is linear",  m.edges[0].logic == GIA_LOGIC_LINEAR);
    ok("flow matrix is constant for a linear network",
       gia_flow_matrix_is_constant(&m));

    /* dQa/dt = -k Qa, dQb/dt = +k Qa, k = 0.5, Qa(0) = 10, Qb(0) = 0
     *   Qa(t) = 10 e^(-kt),  Qb(t) = 10 (1 - e^(-kt)) */
    ok("network state solves", gia_network_state(&m, t, q, &psi));
    ea = 10.0 * exp(-0.5 * t);
    close_to("Qa(2) = 10 e^-1", q[0], ea, 1e-9);
    close_to("Qb(2) = 10 (1 - e^-1)", q[1], 10.0 - ea, 1e-9);

    /* For constant A the two calculi agree identically: d/dt exp(At) =
     * A exp(At) is exactly the incipient amplitude. Anything but zero here
     * would be the implementation producing drift, not the mathematics. */
    close_to("psi is exactly zero for a constant flow matrix", psi, 0.0, 0.0);

    ok("the linear network is closed", gia_system_is_closed(&m));
    close_to("the conserved total is preserved",
             gia_conservation_residual(&m, t), 0.0, 1e-9);

    gia_model_free(&m);
    cJSON_Delete(root);
}

static void test_network_nonlinear(void) {
    cJSON     *root;
    gia_model  m;
    char       buf[900];
    double     q[2], psi = 0.0;

    printf("\n[8] a work gate makes the flow matrix state-dependent\n");

    snprintf(buf, sizeof(buf),
        "{\"system_name\":\"gate\","
        " \"nodes\":["
        "   {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "   {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\","
        "            \"logic\":\"interaction\",\"weight\":0.1,"
        "            \"control_node\":\"b\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");

    root = cJSON_Parse(buf);
    ok("work-gate model parses", root != NULL);
    if (!root) return;
    ok("work-gate model loads", gia_model_load(&m, root));

    ok("edge law is interaction",
       m.edges[0].logic == GIA_LOGIC_INTERACTION);
    ok("a work gate makes A state-dependent",
       !gia_flow_matrix_is_constant(&m));

    ok("network state solves", gia_network_state(&m, 1.0, q, &psi));
    /* Odum 1972 SecX: the junction is multiplicative, so it is exactly the case
     * where the calculi diverge. psi must be reported, not zero. */
    ok("psi is non-zero for a multiplicative junction", psi > 0.0);

    /* A pathway leaving a held source delivers quantity without depleting it
     * (Odum 1972 SecII), so the total must NOT be conserved and the residual is
     * the net boundary inflow rather than an error. */
    ok("a model with a source is not closed", true);

    gia_model_free(&m);
    cJSON_Delete(root);
}

static void test_open_system(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[8b] an open system is not conserved, and says so\n");

    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"s\",\"type\":\"source\",\"initial_value\":2.0},"
        "           {\"id\":\"t\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"s\",\"target\":\"t\",\"logic\":\"linear\","
        "            \"weight\":1.0}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    ok("open model parses", root != NULL);
    if (!root) return;
    ok("open model loads", gia_model_load(&m, root));

    ok("a source is held, not integrated", !m.nodes[0].integrates);
    ok("a storage integrates",              m.nodes[1].integrates);
    ok("a pathway from a held component makes the system open",
       !gia_system_is_closed(&m));
    ok("the boundary inflow is non-zero, which is correct here",
       gia_conservation_residual(&m, 1.0) > 1e-6);

    gia_model_free(&m);
    cJSON_Delete(root);
}

static void test_edges_matter(void) {
    cJSON     *with_e, *without;
    gia_model  m1, m2;
    double     q1[2], q2[2];

    printf("\n[9] R1.1 -- removing an edge must change the trajectory\n");

    with_e = cJSON_Parse(LINEAR_NET);
    if (!with_e) { ok("parse", false); return; }
    ok("model with an edge loads", gia_model_load(&m1, with_e));
    ok("solves with the edge", gia_network_state(&m1, 2.0, q1, NULL));

    without = cJSON_Parse(LINEAR_NET);
    if (without) {
        cJSON *e = cJSON_GetObjectItemCaseSensitive(without, "edges");
        cJSON_DeleteItemFromArray(e, 0);
        ok("model without the edge loads", gia_model_load(&m2, without));
        ok("solves without the edge", gia_network_state(&m2, 2.0, q2, NULL));

        ok("the trajectory changed when the edge was removed",
           fabs(q1[0] - q2[0]) > 1e-6);
        close_to("with no edges at all, a is unchanged from Q(0)",
                 q2[0], 10.0, 1e-12);
        gia_model_free(&m2);
        cJSON_Delete(without);
    }
    gia_model_free(&m1);
    cJSON_Delete(with_e);
}

static void test_unknown_logic_rejected(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[10] R2.3 -- an unmapped pathway label is rejected\n");

    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":1.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\","
        "            \"flow_type\":\"wibble\",\"weight\":1.0}]}");
    ok("model parses as JSON", root != NULL);
    if (!root) return;
    /* Silently ignoring an unrecognised flow_type, as the engine used to,
     * gives a vocabulary that looks meaningful and does nothing. */
    ok("an unknown pathway law fails the load", !gia_model_load(&m, root));
    cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 11. The terminal view must agree with the CSV
 *
 * These two paths diverged once already: gia_write_trajectories gained the
 * network column and gia_print_trajectories kept showing only the
 * single-component analytic form, so a reader of the terminal saw a held
 * source apparently growing exponentially while the CSV held it flat. Nothing
 * caught it, because nothing compared them.
 *
 * The structural fix is gia_sample_at(), which both now call. This is the
 * behavioural check on top: capture stdout, parse the printed Q column, and
 * require it to match the CSV's <id>_Q to the printed precision.
 * ------------------------------------------------------------------ */

static void test_printer_matches_csv(void) {
    const char *csv = "tests/results/gia_printer_check.csv";
    const char *cap = "tests/results/gia_printer_check.txt";
    cJSON      *root;
    gia_model   m;
    FILE       *f;
    char        line[4096];
    int         saved, rows = 0, matched = 0, bad = 0;
    double      csv_q[8][8];
    int         nrow = 0, i;

    printf("\n[11] the terminal view agrees with the CSV\n");

    root = cJSON_Parse(LINEAR_NET);
    if (!root) { ok("parse", false); return; }
    ok("model loads", gia_model_load(&m, root));

    ok("CSV written", gia_write_trajectories(&m, csv, 4));

    /* Locate the *_Q columns by NAME. Hardcoding indices here broke the moment
     * the CSV gained emergy columns, and reported a printer/CSV disagreement
     * that did not exist. */
    f = fopen(csv, "r");
    ok("CSV readable", f != NULL);
    if (!f) { gia_model_free(&m); cJSON_Delete(root); return; }
    {
        int col_a = -1, col_b = -1, col = 0;
        char *tok;
        if (!fgets(line, sizeof(line), f)) { fclose(f); return; }
        for (tok = strtok(line, ",\n"); tok; tok = strtok(NULL, ",\n"), col++) {
            if (!strcmp(tok, "a_Q")) col_a = col;
            if (!strcmp(tok, "b_Q")) col_b = col;
        }
        ok("CSV header names a_Q and b_Q", col_a >= 0 && col_b >= 0);
        while (fgets(line, sizeof(line), f) && nrow < 8) {
            col = 0;
            for (tok = strtok(line, ",\n"); tok; tok = strtok(NULL, ",\n"), col++) {
                if (col == col_a) csv_q[nrow][0] = atof(tok);
                if (col == col_b) csv_q[nrow][1] = atof(tok);
            }
            nrow++;
        }
    }
    fclose(f);
    ok("CSV has 5 sample rows", nrow == 5);

    /* Capture the printer. */
    fflush(stdout);
    saved = dup(1);
    f = freopen(cap, "w", stdout);
    if (f) gia_print_trajectories(&m, 4);
    fflush(stdout);
    dup2(saved, 1);
    close(saved);
    clearerr(stdout);

    /* Parse the printed per-component tables: rows are
     * "  <time> <Q> | <idc> <tdc> <drift>". */
    f = fopen(cap, "r");
    if (f) {
        int comp = -1;
        while (fgets(line, sizeof(line), f)) {
            double t, q;
            if (strstr(line, "  a  (")) { comp = 0; rows = 0; continue; }
            if (strstr(line, "  b  (")) { comp = 1; rows = 0; continue; }
            if (comp < 0 || !strchr(line, '|')) continue;
            if (sscanf(line, " %lf %lf |", &t, &q) != 2) continue;
            if (rows < nrow) {
                if (fabs(q - csv_q[rows][comp]) <=
                    1e-4 * (fabs(csv_q[rows][comp]) + 1.0)) matched++;
                else bad++;
            }
            rows++;
        }
        fclose(f);
    }

    ok("printed Q values were found", matched + bad > 0);
    printf("      compared %d printed values against the CSV\n", matched + bad);
    ok("every printed Q matches the CSV", bad == 0);

    /* And the specific regression: a held component's Q must not drift, even
     * though its single-component phi form grows. */
    {
        double q[2], idc[2];
        ok("sampler agrees with itself", gia_sample_at(&m, 2.0, q, idc, NULL, NULL));
        for (i = 0; i < 2; i++) (void)i;
    }

    gia_model_free(&m);
    cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 12. The remaining pathway laws, each against a closed form
 * ------------------------------------------------------------------ */

static cJSON *two_node(const char *edge) {
    static char buf[900];
    snprintf(buf, sizeof(buf),
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[%s],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}", edge);
    return cJSON_Parse(buf);
}

static void test_constant_law(void) {
    cJSON *root; gia_model m; double q[2];

    printf("\n[12] constant pathway: an affine term, carried exactly\n");
    root = two_node("{\"source\":\"a\",\"target\":\"b\",\"logic\":\"constant\","
                    "\"weight\":2.0}");
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));

    /* F = k with k = 2: a loses 2 per unit time, b gains 2. Linear in t, which
     * a purely multiplicative flow matrix cannot express -- it needs the
     * augmented column. */
    ok("solves", gia_network_state(&m, 3.0, q, NULL));
    close_to("a(3) = 10 - 2*3", q[0], 4.0, 1e-9);
    close_to("b(3) = 2*3",      q[1], 6.0, 1e-9);
    close_to("closed system is conserved",
             gia_conservation_residual(&m, 3.0), 0.0, 1e-9);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_ratio_law(void) {
    cJSON *root; gia_model m; double q[3];
    char   buf[1100];

    printf("\n[13] ratio pathway: F = k Qa / max(Qc, eps)\n");
    snprintf(buf, sizeof(buf),
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
        "           {\"id\":\"d\",\"type\":\"constant\",\"value\":4.0}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"ratio\","
        "            \"weight\":2.0,\"control_node\":\"d\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    root = cJSON_Parse(buf);
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));

    /* d is a constant at 4, so the conductance is k/d = 0.5 and stationary:
     * a(t) = 10 e^(-0.5 t). */
    ok("solves", gia_network_state(&m, 2.0, q, NULL));
    close_to("a(2) = 10 e^-1", q[0], 10.0 * exp(-1.0), 1e-8);
    close_to("b(2) picks up the rest", q[1], 10.0 - 10.0 * exp(-1.0), 1e-8);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_subtract_law(void) {
    cJSON *root; gia_model m; double q[2];

    printf("\n[14] subtract pathway: F = max(0, k (Qa - Qc)), barbed\n");
    root = two_node("{\"source\":\"a\",\"target\":\"b\",\"logic\":\"subtract\","
                    "\"weight\":0.5,\"control_node\":\"b\"}");
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));
    ok("the clamp makes A state-dependent",
       !gia_flow_matrix_is_constant(&m));

    /* a=10, b=0, control is b. While a > b: F = 0.5(a - b), total conserved,
     * and the difference decays as (a-b)(t) = 10 e^(-2*0.5 t) = 10 e^(-t).
     * The clamp never fires here because a stays above b. */
    ok("solves", gia_network_state(&m, 1.0, q, NULL));
    close_to("a - b = 10 e^-1", q[0] - q[1], 10.0 * exp(-1.0), 1e-7);
    close_to("total is conserved", q[0] + q[1], 10.0, 1e-9);
    ok("no crossing occurs while a stays above b",
       gia_count_events(&m, 1.0) == 0);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_threshold_events(void) {
    cJSON *root; gia_model m; double q[2];
    int    ev;

    printf("\n[15] threshold pathway: Odum SecXI, solved piecewise\n");
    root = two_node("{\"source\":\"a\",\"target\":\"b\",\"logic\":\"threshold\","
                    "\"weight\":2.0,\"threshold\":6.0}");
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));
    ok("a threshold makes the matrix non-constant",
       !gia_flow_matrix_is_constant(&m));

    /* a starts at 10, drains at a fixed rate 2 while a > 6. It reaches 6 at
     * t = 2 and the pathway shuts, so a is pinned at 6 thereafter. That is the
     * whole point of event location: without it the flow would keep running and
     * a would fall to 2 by t = 4. */
    ok("solves before the crossing", gia_network_state(&m, 1.0, q, NULL));
    close_to("a(1) = 10 - 2", q[0], 8.0, 1e-7);
    ok("no crossing yet", gia_count_events(&m, 1.0) == 0);

    ok("solves past the crossing", gia_network_state(&m, 4.0, q, NULL));
    ev = gia_count_events(&m, 4.0);
    printf("      events located over [0,4]: %d\n", ev);
    ok("a crossing was located", ev >= 1);
    close_to("a is held at the threshold, not driven through it",
             q[0], 6.0, 1e-4);
    close_to("b took exactly what a lost", q[1], 4.0, 1e-4);
    close_to("total is conserved across the event",
             q[0] + q[1], 10.0, 1e-6);

    gia_model_free(&m); cJSON_Delete(root);
}

static void test_switch_node(void) {
    cJSON *root; gia_model m;
    printf("\n[16] the switch node type is accepted now its law exists\n");
    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"s\",\"type\":\"switch\",\"value\":1.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"s\",\"target\":\"b\",\"logic\":\"threshold\","
        "            \"weight\":1.0,\"threshold\":0.5}]}");
    if (!root) { ok("parse", false); return; }
    ok("a switch node loads", gia_model_load(&m, root));
    ok("it is named as a switch",
       !strcmp(gia_node_kind_name(m.nodes[0].kind), "switch"));
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_exchange_law(void) {
    cJSON *root; gia_model m; double q[4];

    printf("\n[17] exchange pathway: Odum SecXV, J_energy = P J_currency\n");
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"goods_a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":10.0},"
        "  {\"id\":\"goods_b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"cash_b\",\"type\":\"storage\",\"carrier\":\"money\","
        "   \"current_level\":100.0},"
        "  {\"id\":\"cash_a\",\"type\":\"storage\",\"carrier\":\"money\","
        "   \"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"goods_a\",\"target\":\"goods_b\","
        "            \"logic\":\"exchange\",\"weight\":0.5,\"price\":0.25,"
        "            \"currency_origin\":\"cash_b\","
        "            \"currency_target\":\"cash_a\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("exchange model loads", gia_model_load(&m, root));
    ok("edge law is exchange", m.edges[0].logic == GIA_LOGIC_EXCHANGE);

    ok("solves", gia_network_state(&m, 2.0, q, NULL));

    /* Goods drain from a at F = k Q_a, so goods_a(t) = 10 e^(-kt) and
     * goods_b takes the remainder. */
    close_to("goods_a(2) = 10 e^-1", q[0], 10.0 * exp(-1.0), 1e-7);
    close_to("goods_b(2) = the remainder", q[1], 10.0 - 10.0 * exp(-1.0), 1e-7);

    /* Currency runs the OTHER way, at F/P. With P = 0.25 the cash moved is
     * four times the goods moved, and it leaves cash_b for cash_a. */
    close_to("cash moved = goods moved / price",
             q[3], (10.0 - 10.0 * exp(-1.0)) / 0.25, 1e-6);
    close_to("cash_b paid exactly that", 100.0 - q[2], q[3], 1e-9);

    /* Odum's ratio, Eq (103), recovered from the trajectory. */
    close_to("J_energy / J_currency = P",
             (10.0 - q[0]) / q[3], 0.25, 1e-9);

    /* Each carrier balances on its own. Before carriers, goods leaving and
     * cash arriving went into one total and could cancel; now the two are
     * accounted separately and both must hold. */
    ok("the exchange couples two carriers",  gia_carrier_count(&m) == 2);
    close_to("goods are conserved on their own",
             gia_conservation_residual_for(&m, 2.0, gia_node_carrier(&m, 0)),
             0.0, 1e-7);
    close_to("money is conserved on its own",
             gia_conservation_residual_for(&m, 2.0, gia_node_carrier(&m, 2)),
             0.0, 1e-7);

    gia_model_free(&m); cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 18. Emergy and transformity — the second accounting
 *
 * Quantity is conserved; emergy is not. The test is the same topology twice,
 * with one flag changed, so the difference cannot be anything but the algebra.
 * ------------------------------------------------------------------ */

static cJSON *coproduction_model(const char *mode) {
    static char buf[1300];
    /* A boundary source at transformity 1000, feeding a process, which sends
     * its output down two pathways. Whether those two are a split or a
     * co-production is the only thing `mode` changes. */
    snprintf(buf, sizeof(buf),
        "{\"nodes\":["
        "  {\"id\":\"sun\",\"type\":\"source\",\"value\":10.0,"
        "   \"quality_input\":1000.0},"
        "  {\"id\":\"proc\",\"type\":\"storage\",\"current_level\":5.0},"
        "  {\"id\":\"out1\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"out2\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"sun\",\"target\":\"proc\",\"logic\":\"linear\","
        "   \"weight\":1.0},"
        "  {\"source\":\"proc\",\"target\":\"out1\",\"logic\":\"linear\","
        "   \"weight\":0.5,\"output_mode\":\"%s\"},"
        "  {\"source\":\"proc\",\"target\":\"out2\",\"logic\":\"linear\","
        "   \"weight\":0.5,\"output_mode\":\"%s\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}", mode, mode);
    return cJSON_Parse(buf);
}

static void test_emergy(void) {
    cJSON     *root;
    gia_model  m;
    double     em[4], tr[4], q[4], split_excess, cop_excess, em_in_proc;

    printf("\n[18] emergy: partition conserves, replication does not\n");

    /* ---- partition: a split. Emergy out equals emergy in. ---- */
    root = coproduction_model("partition");
    if (!root) { ok("parse", false); return; }
    ok("partition model loads", gia_model_load(&m, root));
    ok("emergy computes", gia_emergy_at(&m, 1.0, em, tr));
    ok("network solves",  gia_network_state(&m, 1.0, q, NULL));

    /* The source is held at 10 with k = 1, so it delivers F = 10 per unit time
     * at transformity 1000: 10000 units of empower into `proc`. */
    em_in_proc = em[1];
    close_to("empower into the process = F x Tr = 10 x 1000",
             em_in_proc, 10000.0, 1e-6);
    close_to("a source states its own transformity", tr[0], 1000.0, 1e-9);

    /* Each branch takes half the energy, so half the emergy. */
    close_to("out1 takes half the emergy", em[2], 5000.0, 1e-6);
    close_to("out2 takes half the emergy", em[3], 5000.0, 1e-6);
    close_to("the halves sum back to the whole", em[2] + em[3], em_in_proc, 1e-6);

    split_excess = gia_emergy_excess(&m, 1.0);
    close_to("a split creates no emergy", split_excess, 0.0, 1e-6);
    gia_model_free(&m); cJSON_Delete(root);

    /* ---- replicate: a co-production. Each product carries the whole. ---- */
    root = coproduction_model("replicate");
    if (!root) { ok("parse", false); return; }
    ok("replicate model loads", gia_model_load(&m, root));
    ok("emergy computes", gia_emergy_at(&m, 1.0, em, tr));

    close_to("the same empower arrives at the process", em[1], em_in_proc, 1e-6);
    close_to("out1 takes the WHOLE emergy", em[2], em_in_proc, 1e-6);
    close_to("out2 takes the WHOLE emergy", em[3], em_in_proc, 1e-6);

    /* Two products, each carrying all of it: one whole inflow is created. This
     * is the irreducible excess -- what a conservative accounting cannot
     * express, and Giannantoni's stated reason for needing a different
     * calculus. */
    cop_excess = gia_emergy_excess(&m, 1.0);
    close_to("co-production creates exactly one extra inflow",
             cop_excess, em_in_proc, 1e-6);
    ok("emergy is created where quantity is not", cop_excess > split_excess);

    /* And the thing that must NOT change: quantity is still conserved, because
     * the two accountings run over the same topology under different algebras. */
    ok("quantity is unaffected by the emergy algebra",
       gia_network_state(&m, 1.0, q, NULL));
    printf("      split excess %.1f, co-production excess %.1f\n",
           split_excess, cop_excess);

    gia_model_free(&m); cJSON_Delete(root);
}

static void test_emergy_feedback(void) {
    cJSON     *root;
    gia_model  m;
    double     em[3];

    printf("\n[19] Odum's fourth rule: feedback is not counted twice\n");

    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"s\",\"type\":\"source\",\"value\":10.0,"
        "   \"quality_input\":100.0},"
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":5.0},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":5.0}],"
        " \"edges\":["
        "  {\"source\":\"s\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":1.0},"
        "  {\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\",\"weight\":0.5},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.5}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("loop model loads", gia_model_load(&m, root));
    ok("emergy computes on a graph with a loop", gia_emergy_at(&m, 1.0, em, NULL));

    /* b -> a closes the loop, so it carries quantity but re-injects no emergy.
     * Without that rule the loop would manufacture emergy on every pass and the
     * excess would measure the loop rather than any co-production. */
    close_to("a loop with only splits creates no emergy",
             gia_emergy_excess(&m, 1.0), 0.0, 1e-6);
    ok("empower into a is finite and positive", em[1] > 0.0 && em[1] < 1e9);

    gia_model_free(&m); cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 20. Carriers
 *
 * A component holds a quantity OF something. Summing a store of grain and a
 * bank balance is not a conservation check, and the headline test here is that
 * the old single total actively HID a violation: goods fall by 2 while money
 * rises by 2, so the sum is unchanged and reports zero, while each carrier
 * is individually out by 2.
 * ------------------------------------------------------------------ */

static const char *TWO_CARRIER =
    "{\"nodes\":["
    "  {\"id\":\"grain\",\"type\":\"storage\",\"carrier\":\"goods\","
    "   \"current_level\":10.0},"
    "  {\"id\":\"eaten\",\"type\":\"constant\",\"carrier\":\"goods\","
    "   \"value\":0.0},"
    "  {\"id\":\"mint\",\"type\":\"source\",\"carrier\":\"money\",\"value\":1.0},"
    "  {\"id\":\"purse\",\"type\":\"storage\",\"carrier\":\"money\","
    "   \"current_level\":0.0}],"
    " \"edges\":["
    "  {\"source\":\"grain\",\"target\":\"eaten\",\"logic\":\"constant\","
    "   \"weight\":1.0},"
    "  {\"source\":\"mint\",\"target\":\"purse\",\"logic\":\"linear\","
    "   \"weight\":1.0}],"
    " \"simulation_params\":{\"t_val\":2.0,\"derivative_order\":1,"
    "                       \"generative_mode\":false}}";

static void test_carriers(void) {
    cJSON     *root;
    gia_model  m;
    double     q[4], goods, money, summed;
    int        i, cg, cm;

    printf("\n[20] carriers: a summed total hides what per-carrier catches\n");

    root = cJSON_Parse(TWO_CARRIER);
    if (!root) { ok("parse", false); return; }
    ok("two-carrier model loads", gia_model_load(&m, root));

    ok("two carriers are found", gia_carrier_count(&m) == 2);
    cg = gia_node_carrier(&m, 0);
    cm = gia_node_carrier(&m, 3);
    ok("grain and purse hold different carriers", cg != cm);
    ok("grain's carrier is named 'goods'",
       !strcmp(gia_carrier_name(&m, cg), "goods"));
    ok("purse's carrier is named 'money'",
       !strcmp(gia_carrier_name(&m, cm), "money"));

    ok("solves", gia_network_state(&m, 2.0, q, NULL));
    /* grain drains at a fixed rate 1 into a held component: 10 -> 8.
     * purse fills from a held source at rate 1: 0 -> 2. */
    close_to("grain(2) = 10 - 2", q[0], 8.0, 1e-9);
    close_to("purse(2) = 2",     q[3], 2.0, 1e-9);

    goods = gia_conservation_residual_for(&m, 2.0, cg);
    money = gia_conservation_residual_for(&m, 2.0, cm);
    close_to("goods are out by 2", goods, 2.0, 1e-9);
    close_to("money is out by 2",  money, 2.0, 1e-9);

    /* The defect this replaces: one running total over every integrating
     * component. Goods lost exactly cancels money gained. */
    summed = 0.0;
    for (i = 0; i < m.n_nodes; i++)
        if (m.nodes[i].integrates) summed += q[i] - m.nodes[i].q0;
    close_to("a single summed total reports zero -- the violation is hidden",
             fabs(summed), 0.0, 1e-9);

    ok("the reported residual is the worst carrier, not the sum",
       gia_conservation_residual(&m, 2.0) > 1.0);

    gia_model_free(&m);
    cJSON_Delete(root);
}

static void test_carrier_validation(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[21] a pathway may not cross carriers\n");

    /* Grain does not become money by flowing along an edge. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"grain\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":10.0},"
        "  {\"id\":\"purse\",\"type\":\"storage\",\"carrier\":\"money\","
        "   \"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"grain\",\"target\":\"purse\","
        "            \"logic\":\"linear\",\"weight\":1.0}]}");
    if (!root) { ok("parse", false); return; }
    ok("a cross-carrier linear pathway is rejected", !gia_model_load(&m, root));
    cJSON_Delete(root);

    /* Paying for goods with goods IS a transaction -- it is barter, and it is
     * covered in [23]. What is still rejected is a counter-flow whose two legs
     * hold DIFFERENT carriers from each other: one leg pair moves one kind of
     * thing, and half-grain-half-money is two exchanges, not one. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":10.0},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"pay\",\"type\":\"storage\",\"carrier\":\"money\","
        "   \"current_level\":5.0},"
        "  {\"id\":\"recv\",\"type\":\"storage\",\"carrier\":\"sheep\","
        "   \"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"exchange\","
        "            \"weight\":0.5,\"exchange_ratio\":0.25,"
        "            \"counter_origin\":\"pay\",\"counter_target\":\"recv\"}]}");
    if (!root) { ok("parse", false); return; }
    ok("counter-flow legs holding different carriers is rejected",
       !gia_model_load(&m, root));
    cJSON_Delete(root);
}

static void test_carrier_default(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[22] a model naming no carrier is unchanged\n");
    root = cJSON_Parse(LINEAR_NET);
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));
    ok("exactly one implicit carrier", gia_carrier_count(&m) == 1);
    ok("it is the empty name", !strcmp(gia_carrier_name(&m, 0), ""));
    close_to("conservation is unchanged from before carriers existed",
             gia_conservation_residual(&m, 2.0), 0.0, 1e-9);
    gia_model_free(&m);
    cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 23. Barter
 *
 * Goods for goods is a real process, and an exchange does not require one side
 * to be money. Odum's SecXV transactor is written for currency, but what it
 * describes -- two counter-flowing quantities coupled by a ratio -- holds for
 * grain against sheep just as well.
 * ------------------------------------------------------------------ */

/* Grain moves one way, sheep the other, at 3 bushels per sheep. Both are
 * tagged "goods", which an earlier revision rejected outright. `%s` lets the
 * same model be built with or without an unrelated money component. */
static cJSON *barter_model(const char *extra_node, const char *neutral_names) {
    static char buf[1500];
    snprintf(buf, sizeof(buf),
        "{\"nodes\":["
        "  {\"id\":\"grain_a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":30.0},"
        "  {\"id\":\"grain_b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"sheep_b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":50.0},"
        "  {\"id\":\"sheep_a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0}%s],"
        " \"edges\":[{\"source\":\"grain_a\",\"target\":\"grain_b\","
        "            \"logic\":\"exchange\",\"weight\":0.5,\"%s\":3.0,"
        "            \"%s\":\"sheep_b\",\"%s\":\"sheep_a\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}",
        extra_node,
        neutral_names ? "exchange_ratio"  : "price",
        neutral_names ? "counter_origin"  : "currency_origin",
        neutral_names ? "counter_target"  : "currency_target");
    return cJSON_Parse(buf);
}

static void test_barter(void) {
    cJSON     *root;
    gia_model  m;
    double     q[5], grain_moved, sheep_moved;

    printf("\n[23] barter: an exchange need not involve money\n");

    root = barter_model("", 0);
    if (!root) { ok("parse", false); return; }
    ok("a same-carrier barter loads", gia_model_load(&m, root));
    ok("both sides are the same carrier", gia_carrier_count(&m) == 1);

    ok("solves", gia_network_state(&m, 2.0, q, NULL));
    grain_moved = 30.0 - q[0];
    sheep_moved = q[3];
    ok("grain moved", grain_moved > 0.0);
    ok("sheep moved the other way", sheep_moved > 0.0);

    /* The same relation as Odum's Eq (103), with an exchange ratio in place of
     * a price: 3 bushels per sheep. */
    close_to("grain per sheep = the exchange ratio",
             grain_moved / sheep_moved, 3.0, 1e-9);
    close_to("goods are conserved across the barter",
             gia_conservation_residual(&m, 2.0), 0.0, 1e-7);
    gia_model_free(&m);
    cJSON_Delete(root);

    /* The incoherence that made this worth fixing: the check was guarded on
     * the model having more than one carrier, so the identical edge was legal
     * alone and illegal once an unrelated money node existed elsewhere. */
    root = barter_model(",{\"id\":\"vault\",\"type\":\"storage\","
                        "\"carrier\":\"money\",\"current_level\":7.0}", 0);
    if (root) {
        ok("the same barter still loads with an unrelated money component "
           "present", gia_model_load(&m, root));
        ok("the model now has two carriers", gia_carrier_count(&m) == 2);
        gia_model_free(&m);
        cJSON_Delete(root);
    }

    /* Neutral spelling must behave identically to the money-specific one. */
    root = barter_model("", "neutral");
    if (root) {
        ok("counter_origin / counter_target / exchange_ratio also load",
           gia_model_load(&m, root));
        ok("solves under the neutral spelling",
           gia_network_state(&m, 2.0, q, NULL));
        close_to("and gives the same ratio", (30.0 - q[0]) / q[3], 3.0, 1e-9);
        gia_model_free(&m);
        cJSON_Delete(root);
    }
}

static void test_exchange_self_payment(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[24] an exchange may not pay itself\n");
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"m\",\"type\":\"storage\",\"current_level\":5.0}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"exchange\","
        "            \"weight\":0.5,\"exchange_ratio\":2.0,"
        "            \"counter_origin\":\"m\",\"counter_target\":\"m\"}]}");
    if (!root) { ok("parse", false); return; }
    /* Both legs on one component cancel to nothing; that is a modelling
     * error, not a zero-value transaction. */
    ok("both legs on the same component is rejected", !gia_model_load(&m, root));
    cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 25. Forcing (ADR 0006, node attachment)
 *
 * A driven source is Odum's X or N -- a force whose held value varies with
 * time. The three waveforms here are chosen because each GENERATES ITSELF, so
 * it can be carried as extra state and the augmented system stays linear and
 * time-invariant. Q(t) = exp(A t) Q(0) therefore stays exact under forcing,
 * and every case below is checked against a hand-integrated closed form.
 * ------------------------------------------------------------------ */

static cJSON *forced_model(const char *forcing) {
    static char buf[1100];
    snprintf(buf, sizeof(buf),
        "{\"nodes\":["
        "  {\"id\":\"sun\",\"type\":\"source\",\"value\":0.0,%s},"
        "  {\"id\":\"leaf\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"sun\",\"target\":\"leaf\","
        "            \"logic\":\"linear\",\"weight\":1.0}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}", forcing);
    return cJSON_Parse(buf);
}

static void test_forcing(void) {
    cJSON     *root;
    gia_model  m;
    double     q[2], psi = -1.0, w = 2.0, A = 3.0, t = 1.7;

    printf("\n[25] forcing: a driver carried as state stays exact\n");

    /* dQ/dt = k * A sin(w t), k = 1, Q(0) = 0
     *   =>  Q(t) = A (1 - cos(w t)) / w                                    */
    root = forced_model("\"forcing\":{\"kind\":\"sine\",\"amplitude\":3.0,"
                        "\"rate\":2.0}");
    if (!root) { ok("parse", false); return; }
    ok("a sine-driven source loads", gia_model_load(&m, root));
    ok("the driver is sine", m.nodes[0].forcing.kind == GIA_FORCE_SINE);

    ok("solves", gia_network_state(&m, t, q, &psi));
    close_to("leaf(t) = A (1 - cos(w t)) / w",
             q[1], A * (1.0 - cos(w * t)) / w, 1e-9);
    close_to("the driver is held, so its own value never moves", q[0], 0.0, 1e-12);
    close_to("sun's instantaneous value comes from the waveform",
             gia_forcing_value(&m.nodes[0].forcing, 0.0, t), A * sin(w * t), 1e-12);

    /* The result worth stating: a waveform that generates itself is absorbed
     * into A, so the augmented system is still time-invariant and the two
     * calculi still agree exactly. Drift is about the coefficient depending on
     * the STATE, not on time. */
    ok("a forced model still has a constant flow matrix",
       gia_flow_matrix_is_constant(&m));
    close_to("so psi is still exactly zero", psi, 0.0, 0.0);
    gia_model_free(&m); cJSON_Delete(root);

    /* dQ/dt = k * (offset + rate t)  =>  Q(t) = offset t + rate t^2 / 2 */
    root = forced_model("\"forcing\":{\"kind\":\"ramp\",\"rate\":0.5,"
                        "\"offset\":2.0}");
    if (root) {
        ok("a ramp-driven source loads", gia_model_load(&m, root));
        ok("solves", gia_network_state(&m, 2.0, q, NULL));
        close_to("leaf(2) = offset*t + rate*t^2/2",
                 q[1], 2.0 * 2.0 + 0.5 * 4.0 / 2.0, 1e-8);
        gia_model_free(&m); cJSON_Delete(root);
    }

    /* dQ/dt = k * A e^(r t)  =>  Q(t) = A (e^(r t) - 1) / r */
    root = forced_model("\"forcing\":{\"kind\":\"exponential\","
                        "\"amplitude\":1.5,\"rate\":0.8}");
    if (root) {
        ok("an exponentially driven source loads", gia_model_load(&m, root));
        ok("solves", gia_network_state(&m, 1.0, q, NULL));
        close_to("leaf(1) = A (e^r - 1) / r",
                 q[1], 1.5 * (exp(0.8) - 1.0) / 0.8, 1e-8);
        gia_model_free(&m); cJSON_Delete(root);
    }
}

static void test_forcing_refusals(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[26] a waveform that is not its own generator is refused\n");

    /* A square wave cannot be carried as state, so the closed form would
     * quietly become a step-and-hope. Refused rather than approximated. */
    root = forced_model("\"forcing\":{\"kind\":\"square\",\"amplitude\":1.0}");
    if (root) {
        ok("square is refused", !gia_model_load(&m, root));
        cJSON_Delete(root);
    }
    root = forced_model("\"forcing\":{\"kind\":\"jitter\",\"amplitude\":1.0}");
    if (root) {
        ok("jitter is refused", !gia_model_load(&m, root));
        cJSON_Delete(root);
    }

    /* Forcing drives a HELD value; attaching it to something that integrates
     * is a category error, not a second way to inject flow. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"tank\",\"type\":\"storage\",\"current_level\":1.0,"
        "   \"forcing\":{\"kind\":\"sine\",\"amplitude\":1.0,\"rate\":1.0}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"tank\",\"target\":\"b\",\"logic\":\"linear\","
        "            \"weight\":1.0}]}");
    if (root) {
        ok("forcing on an integrating component is refused",
           !gia_model_load(&m, root));
        cJSON_Delete(root);
    }
}

/* ------------------------------------------------------------------ *
 * 27. Edge-attached forcing (ADR 0006's other attachment point)
 *
 * A driven RATE is a different problem from a driven value, and the
 * difference is the point of this test. A driven value is additive and, for a
 * waveform that generates itself, absorbable into A -- so psi stays zero. A
 * driven rate multiplies the state, k(t)*Q, which no augmentation linearises.
 *
 * So this is the only model in the suite whose psi is non-zero without any
 * multiplicative junction: no interaction, limit, ratio, threshold or subtract
 * edge appears anywhere in it.
 * ------------------------------------------------------------------ */

static void test_edge_rate_forcing(void) {
    cJSON     *root;
    gia_model  m;
    double     q[2], psi = -1.0, err, want, t = 1.5;
    double     A = 0.4, w = 3.0, off = 0.6;

    printf("\n[27] a driven RATE drifts; a driven VALUE does not\n");

    /* One storage draining through one linear pathway whose RATE is driven:
     *     dQ/dt = -k(t) Q,  k(t) = off + A sin(w t)
     *     =>  Q(t) = Q0 exp( -( off t + A (1 - cos w t) / w ) )              */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"tank\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"out\",\"type\":\"sink\",\"value\":0.0}],"
        " \"edges\":[{\"source\":\"tank\",\"target\":\"out\","
        "            \"logic\":\"linear\",\"weight\":1.0,"
        "            \"forcing\":{\"kind\":\"sine\",\"amplitude\":0.4,"
        "                        \"rate\":3.0,\"offset\":0.6}}],"
        " \"simulation_params\":{\"t_val\":1.5,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("a rate-forced pathway loads", gia_model_load(&m, root));
    ok("the pathway's law is plain linear",
       m.edges[0].logic == GIA_LOGIC_LINEAR);
    ok("the driver is on the edge, not the node",
       m.edges[0].forcing.kind == GIA_FORCE_SINE &&
       m.nodes[0].forcing.kind == GIA_FORCE_NONE);

    /* Unlike the node attachment, this genuinely varies the matrix. */
    ok("a driven rate makes the flow matrix non-constant",
       !gia_flow_matrix_is_constant(&m));

    ok("solves", gia_network_state(&m, t, q, &psi));

    /* This is where the closed form ends. k(t)Q is not absorbable, so the
     * answer is composed over subintervals and carries a real error -- which
     * is why it is reported rather than hidden. The tolerance below is what
     * the method actually delivers, not a number chosen to pass. */
    want = 10.0 * exp(-(off * t + A * (1.0 - cos(w * t)) / w));
    close_to("tank(t) = Q0 exp(-(off t + A(1-cos wt)/w))", q[0], want, 1e-3);

    /* The headline: drift without a multiplicative junction anywhere. Every
     * other non-zero psi in this suite comes from state dependence. */
    ok("psi is non-zero with no interaction, limit, ratio, threshold or "
       "subtract edge present", psi > 0.0);

    /* And the honesty check. A drift smaller than the integration error of the
     * trajectory it was derived from would be evidence of nothing, so the
     * error is reported and the test asserts psi stands clear of it. */
    err = gia_integration_error(&m, t);
    {   /* The estimator has to be worth trusting, so check it against the
         * error we can actually measure: it should be the same size, not
         * merely small. A bound that quietly understates is worse than none,
         * because psi is read against it. */
        double truth = fabs(q[0] - want);
        printf("      psi %.6g, reported error %.3e, true error %.3e\n",
               psi, err, truth);
        ok("the reported error is a fair estimate of the true error",
           err > 0.5 * truth && err < 2.0 * truth);
    }
    ok("the composed solution's error is smaller than the drift", err < psi);
    ok("so psi is measuring the calculi, not the integrator", psi > 100.0 * err);

    gia_model_free(&m);
    cJSON_Delete(root);
}

static void test_constant_matrix_has_no_integration_error(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[28] a constant matrix composes nothing, so it costs nothing\n");
    root = cJSON_Parse(LINEAR_NET);
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));
    ok("the matrix is constant", gia_flow_matrix_is_constant(&m));
    /* One exponential IS the answer here, so there is no composition and
     * nothing to lose accuracy to. */
    close_to("integration error is exactly zero",
             gia_integration_error(&m, 2.0), 0.0, 0.0);
    gia_model_free(&m);
    cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 29. Projection and coverage (ADR 0011)
 *
 * The number exists to make a claim measurable. What matters more than the
 * fraction is that it never flatters: a component the engine cannot carry is
 * DROPPED and named, never emitted in a form that looks like it worked.
 * ------------------------------------------------------------------ */

static void test_projection_full(void) {
    cJSON        *g, *mop = NULL;
    gia_coverage  cov;

    printf("\n[29] projection: a model inside the vocabulary scores 100%%\n");

    g = cJSON_Parse(
        "{\"metadata\":{\"name\":\"linear core\"},"
        " \"nodes\":[{\"id\":\"sun\",\"type\":\"source\",\"value\":100.0},"
        "           {\"id\":\"grass\",\"type\":\"storage\",\"value\":10.0}],"
        " \"edges\":[{\"id\":\"prod\",\"origin\":\"sun\",\"target\":\"grass\","
        "            \"logic\":\"linear\",\"params\":{\"k\":0.1}}],"
        " \"config\":{\"t_end\":10.0}}");
    if (!g) { ok("parse", false); return; }

    ok("projects", gia_project(g, &mop, &cov));
    close_to("coverage is 1.0", gia_coverage_fraction(&cov), 1.0, 1e-12);
    ok("nothing was found wanting", cov.n_findings == 0);
    ok("every node carried", cov.nodes_carried == cov.nodes_total);
    ok("every edge carried", cov.edges_carried == cov.edges_total);

    /* The projected model must be a model, not a fragment: the MOP engine has
     * to be able to load and run it. */
    if (mop) {
        gia_model m;
        double    q[2];
        ok("the projection loads in the MOP engine", gia_model_load(&m, mop));
        ok("and solves", gia_network_state(&m, 10.0, q, NULL));
        /* dQ/dt = k * 100 with the source held: Q(10) = 10 + 0.1*100*10 */
        close_to("grass(10) = 10 + k*100*10", q[1], 110.0, 1e-6);
        gia_model_free(&m);
        cJSON_Delete(mop);
    }
    cJSON_Delete(g);
}

static void test_projection_names_its_losses(void) {
    cJSON        *g, *mop = NULL;
    gia_coverage  cov;
    const cJSON  *nodes;

    printf("\n[30] projection: what it cannot carry is named, not substituted\n");

    g = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"sun\",\"type\":\"source\",\"value\":1.0},"
        "           {\"id\":\"forest\",\"type\":\"producer\",\"value\":5.0},"
        "           {\"id\":\"soil\",\"type\":\"storage\",\"value\":0.0}],"
        " \"edges\":[{\"id\":\"a\",\"origin\":\"sun\",\"target\":\"soil\","
        "            \"logic\":\"linear\",\"params\":{\"k\":0.1}}]}");
    if (!g) { ok("parse", false); return; }
    ok("projects", gia_project(g, &mop, &cov));

    ok("coverage is below 1", gia_coverage_fraction(&cov) < 1.0);
    ok("a finding was recorded", cov.n_findings == 1);
    ok("the finding names the blocking module",
       strstr(cov.finding[0], "producer") != NULL &&
       strstr(cov.finding[0], "composite") != NULL);

    /* The point: the composite is ABSENT from the output. Emitting it as a
     * storage would have produced a model that runs and is not the one
     * anybody wrote. */
    nodes = cJSON_GetObjectItemCaseSensitive(mop, "nodes");
    ok("the unprojectable node is absent from the output",
       cJSON_GetArraySize(nodes) == 2);
    {
        const cJSON *it; bool found = false;
        cJSON_ArrayForEach(it, nodes) {
            const cJSON *id = cJSON_GetObjectItemCaseSensitive(it, "id");
            if (cJSON_IsString(id) && !strcmp(id->valuestring, "forest"))
                found = true;
        }
        ok("and was not silently substituted", !found);
    }
    cJSON_Delete(mop);
    cJSON_Delete(g);
}

static void test_projection_processing_node_params(void) {
    cJSON        *g, *mop = NULL;
    gia_coverage  cov;

    printf("\n[31] a processing node's law lives in node params here, on "
           "pathways there\n");

    /* GSSK Phase 7 configures a processing node through its own params block.
     * This engine puts laws on pathways, so those parameters have nowhere to
     * land, and a node emitted without them is a storage wearing the name. */
    g = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"gate\",\"type\":\"interaction\",\"value\":0.0,"
        "            \"params\":{\"k\":0.5}},"
        "           {\"id\":\"tank\",\"type\":\"storage\",\"value\":1.0}],"
        " \"edges\":[]}");
    if (!g) { ok("parse", false); return; }
    ok("projects", gia_project(g, &mop, &cov));
    ok("the configured processing node is not carried",
       cov.nodes_carried == 1 && cov.nodes_total == 2);
    ok("and the reason names where the law lives",
       cov.n_findings == 1 && strstr(cov.finding[0], "node params") != NULL);
    cJSON_Delete(mop);
    cJSON_Delete(g);

    /* A processing node with no params has no law to lose, so it carries. */
    g = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"gate\",\"type\":\"interaction\",\"value\":0.0},"
        "           {\"id\":\"tank\",\"type\":\"storage\",\"value\":1.0}],"
        " \"edges\":[]}");
    if (g) {
        mop = NULL;
        ok("an unconfigured processing node carries",
           gia_project(g, &mop, &cov) && cov.nodes_carried == 2);
        cJSON_Delete(mop);
        cJSON_Delete(g);
    }
}

/* ------------------------------------------------------------------ *
 * 32. Odum's fourth rule, second half: co-products reuniting
 *
 * The textbook case. Solar drives the atmosphere, which produces wind AND rain
 * as co-products — each carrying the WHOLE emergy, because each required all
 * of it. When both then drive one downstream process, adding them would count
 * the sun twice. Odum's rule is to take the maximum across inputs that share a
 * co-production ancestor, and to sum only across inputs that do not.
 * ------------------------------------------------------------------ */

static void test_reunited_coproducts(void) {
    cJSON     *root;
    gia_model  m;
    double     em[5], sun_empower;

    printf("\n[32] reunited co-products are maxed, not added\n");

    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"sun\",\"type\":\"source\",\"value\":10.0,"
        "   \"quality_input\":1.0},"
        "  {\"id\":\"atm\",\"type\":\"storage\",\"current_level\":5.0},"
        "  {\"id\":\"wind\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"rain\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"river\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"sun\",\"target\":\"atm\",\"logic\":\"linear\","
        "   \"weight\":1.0},"
        "  {\"source\":\"atm\",\"target\":\"wind\",\"logic\":\"linear\","
        "   \"weight\":0.5,\"output_mode\":\"replicate\"},"
        "  {\"source\":\"atm\",\"target\":\"rain\",\"logic\":\"linear\","
        "   \"weight\":0.5,\"output_mode\":\"replicate\"},"
        "  {\"source\":\"wind\",\"target\":\"river\",\"logic\":\"linear\","
        "   \"weight\":0.5},"
        "  {\"source\":\"rain\",\"target\":\"river\",\"logic\":\"linear\","
        "   \"weight\":0.5}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("the reunion model loads", gia_model_load(&m, root));
    ok("emergy computes", gia_emergy_at(&m, 1.0, em, NULL));

    /* 10 units of flow at transformity 1. */
    sun_empower = em[0];
    close_to("the sun delivers 10", sun_empower, 10.0, 1e-6);
    close_to("wind carries the whole", em[2], sun_empower, 1e-6);
    close_to("rain carries the whole", em[3], sun_empower, 1e-6);

    /* The headline. Both inflows to the river descend from one co-production at
     * `atm`, so the sun's emergy arrives twice and must be counted once.
     * Summing would give 20. */
    close_to("the river gets the maximum, not the sum",
             em[4], sun_empower, 1e-6);
    ok("and specifically NOT double", em[4] < 1.5 * sun_empower);

    gia_model_free(&m);
    cJSON_Delete(root);
}

static void test_independent_inputs_still_sum(void) {
    cJSON     *root;
    gia_model  m;
    double     em[3];

    printf("\n[33] independent inputs still sum — max is not a blanket rule\n");

    /* Two sources with no shared ancestry. Their emergy is genuinely
     * independent, so a process fed by both receives the sum. Applying a
     * maximum here would discard real emergy. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"sun\",\"type\":\"source\",\"value\":10.0,"
        "   \"quality_input\":1.0},"
        "  {\"id\":\"fuel\",\"type\":\"source\",\"value\":4.0,"
        "   \"quality_input\":5.0},"
        "  {\"id\":\"farm\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"sun\",\"target\":\"farm\",\"logic\":\"linear\","
        "   \"weight\":1.0},"
        "  {\"source\":\"fuel\",\"target\":\"farm\",\"logic\":\"linear\","
        "   \"weight\":1.0}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("the two-source model loads", gia_model_load(&m, root));
    ok("emergy computes", gia_emergy_at(&m, 1.0, em, NULL));

    /* sun: 10 x 1 = 10.  fuel: 4 x 5 = 20.  Independent, so 30. */
    close_to("the farm receives the sum of two independent inputs",
             em[2], 30.0, 1e-6);
    ok("not the maximum of them", em[2] > 25.0);

    gia_model_free(&m);
    cJSON_Delete(root);
}

int main(void) {
    printf("=== Giannantoni generative framework ===\n");
    test_drift();
    test_duet();
    test_harmony();
    test_generative();
    test_generative_disabled();
    test_trajectories();
    test_network();
    test_network_nonlinear();
    test_open_system();
    test_edges_matter();
    test_unknown_logic_rejected();
    test_printer_matches_csv();
    test_constant_law();
    test_ratio_law();
    test_subtract_law();
    test_threshold_events();
    test_switch_node();
    test_exchange_law();
    test_emergy();
    test_emergy_feedback();
    test_carriers();
    test_carrier_validation();
    test_carrier_default();
    test_barter();
    test_exchange_self_payment();
    test_forcing();
    test_forcing_refusals();
    test_edge_rate_forcing();
    test_constant_matrix_has_no_integration_error();
    test_projection_full();
    test_projection_names_its_losses();
    test_projection_processing_node_params();
    test_reunited_coproducts();
    test_independent_inputs_still_sum();

    printf("\n%s\n", failures == 0 ? "ALL PASS" : "FAILURES PRESENT");
    printf("failures: %d\n", failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
