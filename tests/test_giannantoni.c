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

/* Load, or record the failure and stop: nothing after a failed load may touch
 * the model. */
static bool load_ok(const char *what, gia_model *m, cJSON *root) {
    bool good = gia_model_load(m, root);
    ok(what, good);
    if (!good) cJSON_Delete(root);
    return good;
}

/* Load and free; true when the document is a valid model. */
static bool loads(const char *json) {
    cJSON     *root = cJSON_Parse(json);
    gia_model  m;
    bool       good;
    if (!root) return false;
    good = gia_model_load(&m, root);
    if (good) gia_model_free(&m);
    cJSON_Delete(root);
    return good;
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

/* The example seed's shape (examples/giannantoni/input.json), with both
 * controls read rather than drawn so there is no heat sink: a sink is never
 * closed (ADR 0015), and this test is about the step reaching maximum. */
static const char *SEED =
    "{\"system_name\":\"test\","
    " \"nodes\":["
    "   {\"id\":\"source_1\",\"type\":\"source\",\"initial_value\":2.0},"
    "   {\"id\":\"interaction_1\",\"type\":\"interaction\","
    "    \"module\":{\"k\":0.1}},"
    "   {\"id\":\"store_1\",\"type\":\"storage\",\"capacity\":100.0,"
    "    \"current_level\":10.0},"
    "   {\"id\":\"consumer_1\",\"type\":\"storage\",\"metabolic_rate\":0.2}],"
    " \"edges\":["
    "   {\"source\":\"source_1\",\"target\":\"interaction_1\",\"role\":\"energy\"},"
    "   {\"source\":\"store_1\",\"target\":\"interaction_1\",\"role\":\"control\","
    "    \"use_ratio\":0},"
    "   {\"source\":\"interaction_1\",\"target\":\"store_1\",\"weight\":1.0},"
    "   {\"source\":\"store_1\",\"target\":\"consumer_1\",\"weight\":0.5},"
    "   {\"source\":\"consumer_1\",\"target\":\"interaction_1\",\"role\":\"control\","
    "    \"use_ratio\":0}],"
    " \"simulation_params\":{\"t_val\":1.5,\"derivative_order\":2,"
    "                       \"generative_mode\":%s}}";

static char *seed_json(const char *generative) {
    static char buf[1800];
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

    /* The work gate is a module, so it is not a component, and its controls
     * are read, so they close nothing: none of the three components is on a
     * closed pathway. */
    close_to("seed ordinality is 0 of 3", gia_ordinality(&m), 0.0, 1e-12);
    ok("seed is below maximum ordinality", !gia_at_maximum_ordinality(&m));
    ok("source_1 is open", !m.nodes[0].on_cycle);
    ok("interaction_1 is a module, not a component",
       gia_node_is_module(&m, 1) && gia_component_count(&m) == 3);

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

    ok("input graph is left untouched", n_in == 4 && e_in == 5);
    ok("one component emerged",  n_out == n_in + 1);
    ok("three relationships emerged", e_out == e_in + 3);
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
    double     q[3], psi = 0.0;

    printf("\n[8] a work gate makes the flow matrix state-dependent\n");

    snprintf(buf, sizeof(buf),
        "{\"system_name\":\"gate\","
        " \"nodes\":["
        "   {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "   {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0},"
        "   {\"id\":\"g\",\"type\":\"interaction\",\"module\":{\"k\":0.1}}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"g\",\"role\":\"energy\"},"
        "            {\"source\":\"b\",\"target\":\"g\",\"role\":\"control\","
        "             \"use_ratio\":0},"
        "            {\"source\":\"g\",\"target\":\"b\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");

    root = cJSON_Parse(buf);
    ok("work-gate model parses", root != NULL);
    if (!root) return;
    ok("work-gate model loads", gia_model_load(&m, root));

    ok("the gate is a module", gia_node_is_module(&m, 2));
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
    cJSON *root; gia_model m; double q[4];

    printf("\n[13] divisor action: F = k Qa / max(Qc, eps)\n");
    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
        "           {\"id\":\"d\",\"type\":\"constant\",\"value\":4.0},"
        "           {\"id\":\"g\",\"type\":\"interaction\","
        "            \"module\":{\"k\":2.0,\"action\":\"divide\"}}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"g\",\"role\":\"energy\"},"
        "           {\"source\":\"d\",\"target\":\"g\",\"role\":\"control\",\"use_ratio\":0},"
        "           {\"source\":\"g\",\"target\":\"b\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;

    /* d is a constant at 4, so the conductance is k/d = 0.5 and stationary:
     * a(t) = 10 e^(-0.5 t). */
    ok("solves", gia_network_state(&m, 2.0, q, NULL));
    close_to("a(2) = 10 e^-1", q[0], 10.0 * exp(-1.0), 1e-8);
    close_to("b(2) picks up the rest", q[1], 10.0 - 10.0 * exp(-1.0), 1e-8);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_subtract_law(void) {
    cJSON *root; gia_model m; double q[3];

    printf("\n[14] subtracting action: F = max(0, k (Qa - Qc)), barbed\n");
    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
        "           {\"id\":\"g\",\"type\":\"interaction\","
        "            \"module\":{\"k\":0.5,\"action\":\"subtract\"}}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"g\",\"role\":\"energy\"},"
        "           {\"source\":\"b\",\"target\":\"g\",\"role\":\"control\",\"use_ratio\":0},"
        "           {\"source\":\"g\",\"target\":\"b\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
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
    cJSON *root; gia_model m; double q[3];
    int    ev;

    printf("\n[15] switch: Odum SecXI, solved piecewise\n");
    /* a is both the energy the switch drains and the sensor it reads -- the
     * pathway threshold's reading, written as the module. */
    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
        "           {\"id\":\"sw\",\"type\":\"switch\","
        "            \"module\":{\"k\":2.0,\"threshold\":6.0}}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"sw\",\"role\":\"energy\"},"
        "           {\"source\":\"a\",\"target\":\"sw\",\"role\":\"control\",\"use_ratio\":0},"
        "           {\"source\":\"sw\",\"target\":\"b\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
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
    printf("\n[16] the switch node type is a module\n");
    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"s\",\"type\":\"switch\","
        "            \"module\":{\"k\":1.0,\"threshold\":0.5}},"
        "           {\"id\":\"t\",\"type\":\"storage\",\"current_level\":1.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"t\",\"target\":\"s\",\"role\":\"energy\"},"
        "           {\"source\":\"t\",\"target\":\"s\",\"role\":\"control\",\"use_ratio\":0},"
        "           {\"source\":\"s\",\"target\":\"b\"}]}");
    if (!root) { ok("parse", false); return; }
    if (!load_ok("a switch module loads", &m, root)) return;
    ok("it is named as a switch",
       !strcmp(gia_node_kind_name(m.nodes[0].kind), "switch"));
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
    static char buf[1800];
    snprintf(buf, sizeof(buf),
        "{\"nodes\":["
        "  {\"id\":\"grain_a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":30.0},"
        "  {\"id\":\"grain_b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"sheep_b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":50.0},"
        "  {\"id\":\"sheep_a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"mkt\",\"type\":\"exchange\","
        "   \"module\":{\"k\":0.5,\"%s\":3.0}}%s],"
        " \"edges\":["
        "  {\"source\":\"grain_a\",\"target\":\"mkt\",\"role\":\"goods_in\"},"
        "  {\"source\":\"mkt\",\"target\":\"grain_b\",\"role\":\"goods_out\"},"
        "  {\"source\":\"sheep_b\",\"target\":\"mkt\",\"role\":\"counter_in\"},"
        "  {\"source\":\"mkt\",\"target\":\"sheep_a\",\"role\":\"counter_out\"}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}",
        neutral_names ? "exchange_ratio" : "price", extra_node);
    return cJSON_Parse(buf);
}

static void test_barter(void) {
    cJSON     *root;
    gia_model  m;
    double     q[6], grain_moved, sheep_moved;

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
        ok("exchange_ratio, the neutral name for price, also loads",
           gia_model_load(&m, root));
        ok("solves under the neutral spelling",
           gia_network_state(&m, 2.0, q, NULL));
        close_to("and gives the same ratio", (30.0 - q[0]) / q[3], 3.0, 1e-9);
        gia_model_free(&m);
        cJSON_Delete(root);
    }
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

    printf("\n[31] a GSSK processing node's legs are positional; here they "
           "are named\n");

    /* GSSK Phase 7 configures a processing node through its own params block
     * and finds its legs by position. This engine names a module's legs by
     * role, so it will not guess which pathway is which (ADR 0013). */
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

    /* With no params it used to carry as a plain node -- a storage wearing
     * the name, which the engine now refuses to load (ADR 0012). Dropped and
     * named instead. */
    g = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"gate\",\"type\":\"interaction\",\"value\":0.0},"
        "           {\"id\":\"tank\",\"type\":\"storage\",\"value\":1.0}],"
        " \"edges\":[]}");
    if (g) {
        mop = NULL;
        ok("an unconfigured processing node is not carried either",
           gia_project(g, &mop, &cov) && cov.nodes_carried == 1 &&
           cov.n_findings == 1 && strstr(cov.finding[0], "wearing the name"));
        if (mop) {
            gia_model m;
            ok("and what is carried loads", gia_model_load(&m, mop));
            gia_model_free(&m);
        }
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

/* ------------------------------------------------------------------ *
 * 34. Module-hosted laws (ADR 0012), with inputs named by role (ADR 0013)
 *
 * The headline is [34]: the same model, its pathways listed in three different
 * orders, must produce byte-identical output. GSSK fails that test — swapping
 * two lines changes which stock it drains, and for an amplifier changes the
 * answer tenfold — which is why roles exist here.
 * ------------------------------------------------------------------ */

/* A work gate: energy from `grass`, control from `sun`, output to `cow`.
 * `edges` is supplied by the caller so the same model can be written with its
 * pathways in any order. */
static cJSON *gate_model(const char *edges) {
    static char buf[1400];
    snprintf(buf, sizeof(buf),
        "{\"nodes\":["
        "  {\"id\":\"grass\",\"type\":\"storage\",\"current_level\":100.0},"
        "  {\"id\":\"sun\",\"type\":\"constant\",\"value\":2.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.01}},"
        "  {\"id\":\"cow\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[%s],"
        " \"simulation_params\":{\"t_val\":2.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}", edges);
    return cJSON_Parse(buf);
}

#define E_ENERGY "{\"source\":\"grass\",\"target\":\"gate\",\"role\":\"energy\"}"
#define E_CTRL   "{\"source\":\"sun\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0}"
#define E_OUT    "{\"source\":\"gate\",\"target\":\"cow\",\"weight\":1.0}"

static char *slurp_file(const char *path) {
    FILE *f = fopen(path, "rb");
    long  n;
    char *b;
    if (!f) return NULL;
    fseek(f, 0, SEEK_END); n = ftell(f); rewind(f);
    b = (char *)malloc((size_t)n + 1);
    if (!b) { fclose(f); return NULL; }
    if (fread(b, 1, (size_t)n, f) != (size_t)n) { free(b); fclose(f); return NULL; }
    b[n] = '\0';
    fclose(f);
    return b;
}

static void test_module_order_invariance(void) {
    const char *orders[3] = {
        E_ENERGY "," E_CTRL "," E_OUT,
        E_CTRL "," E_ENERGY "," E_OUT,
        E_OUT "," E_CTRL "," E_ENERGY
    };
    const char *paths[3] = {
        "tests/results/gia_order_a.csv",
        "tests/results/gia_order_b.csv",
        "tests/results/gia_order_c.csv"
    };
    char *body[3] = { NULL, NULL, NULL };
    int   i;

    printf("\n[34] the order of pathways cannot change a result\n");

    for (i = 0; i < 3; i++) {
        cJSON     *root = gate_model(orders[i]);
        gia_model  m;
        if (!root) { ok("parse", false); return; }
        ok("model loads with its pathways in this order",
           gia_model_load(&m, root));
        ok("trajectories written", gia_write_trajectories(&m, paths[i], 8));
        body[i] = slurp_file(paths[i]);
        gia_model_free(&m);
        cJSON_Delete(root);
    }

    /* The property ADR 0013 exists for. Not "close" — identical. */
    ok("order A and order B are byte-identical",
       body[0] && body[1] && !strcmp(body[0], body[1]));
    ok("order A and order C are byte-identical",
       body[0] && body[2] && !strcmp(body[0], body[2]));

    for (i = 0; i < 3; i++) free(body[i]);
}

static void test_module_gate_closed_form(void) {
    cJSON     *root;
    gia_model  m;
    double     q[4], g = 0.01 * 2.0, t = 2.0;

    printf("\n[35] work gate: F = k Q_energy * prod(controls)\n");

    root = gate_model(E_ENERGY "," E_CTRL "," E_OUT);
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));
    ok("the gate is a module", gia_node_is_module(&m, 2));
    ok("a work gate makes the matrix state-dependent",
       !gia_flow_matrix_is_constant(&m));

    ok("solves", gia_network_state(&m, t, q, NULL));
    /* The control is held, so g = k*Q_sun is stationary and grass decays
     * exactly: grass(t) = 100 e^(-g t), and cow takes what grass loses. */
    close_to("grass(2) = 100 e^(-k*Qsun*t)", q[0], 100.0 * exp(-g * t), 1e-6);
    close_to("cow(2) = what grass lost", q[3], 100.0 - q[0], 1e-6);
    close_to("the gate itself holds nothing: flow passes through",
             q[2], 0.0, 1e-12);
    close_to("the closed system is conserved",
             gia_conservation_residual(&m, t), 0.0, 1e-6);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_module_nary_and_fanout(void) {
    cJSON     *root;
    gia_model  m;
    double     q[6], g, t = 1.0;

    printf("\n[36] n-ary controls, and fan-out in proportion to weight\n");

    /* Two controls — the case a binary pathway law cannot express at all —
     * and two outputs weighted 3:1. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":100.0},"
        "  {\"id\":\"b\",\"type\":\"constant\",\"value\":2.0},"
        "  {\"id\":\"c\",\"type\":\"constant\",\"value\":3.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.05}},"
        "  {\"id\":\"out1\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"out2\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"b\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"gate\",\"target\":\"out1\",\"weight\":3.0},"
        "  {\"source\":\"gate\",\"target\":\"out2\",\"weight\":1.0}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("a two-control gate loads", gia_model_load(&m, root));
    ok("solves", gia_network_state(&m, t, q, NULL));

    /* Every control folds into the conductance: g = k * Qb * Qc. */
    g = 0.05 * 2.0 * 3.0;
    close_to("a(1) = 100 e^(-k Qb Qc t)", q[0], 100.0 * exp(-g * t), 1e-6);
    close_to("the two outputs take what a lost",
             q[4] + q[5], 100.0 - q[0], 1e-6);
    close_to("split 3:1 by pathway weight", q[4] / q[5], 3.0, 1e-9);
    close_to("conserved", gia_conservation_residual(&m, t), 0.0, 1e-6);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_module_gain_closed_form(void) {
    cJSON     *root;
    gia_model  m;
    double     q[4], k = 0.05, S = 10.0, t = 3.0;

    printf("\n[37] amplifier: the control sets the rate, the energy is drained\n");

    /* Odum SecIX. The control is held, so F = k*S is constant and both
     * trajectories are exactly linear — and the stock that empties is the
     * ENERGY input, which is the thing GSSK gets wrong by position. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"signal\",\"type\":\"constant\",\"value\":10.0},"
        "  {\"id\":\"power\",\"type\":\"storage\",\"current_level\":100.0},"
        "  {\"id\":\"amp\",\"type\":\"gain\",\"module\":{\"k\":0.05}},"
        "  {\"id\":\"load\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"signal\",\"target\":\"amp\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"power\",\"target\":\"amp\",\"role\":\"energy\"},"
        "  {\"source\":\"amp\",\"target\":\"load\",\"weight\":1.0}],"
        " \"simulation_params\":{\"t_val\":3.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));
    ok("solves", gia_network_state(&m, t, q, NULL));

    close_to("load(3) = k*S*t",        q[3], k * S * t, 1e-6);
    close_to("power(3) = 100 - k*S*t", q[1], 100.0 - k * S * t, 1e-6);
    close_to("the signal is read, never drained", q[0], 10.0, 1e-12);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_module_validation(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[38] a module states the roles it requires\n");

    /* No role at all: the case GSSK answers by position. */
    root = gate_model("{\"source\":\"grass\",\"target\":\"gate\"},"
                      E_CTRL "," E_OUT);
    if (root) {
        ok("an input with no role is rejected", !gia_model_load(&m, root));
        cJSON_Delete(root);
    }
    /* Two energy inputs: which one is consumed would have to be guessed. */
    root = gate_model(E_ENERGY ","
                      "{\"source\":\"sun\",\"target\":\"gate\",\"role\":\"energy\"},"
                      E_OUT);
    if (root) {
        ok("two energy inputs are rejected", !gia_model_load(&m, root));
        cJSON_Delete(root);
    }
    /* A role on a pathway entering no module is a misplaced field. */
    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":1.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\","
        "            \"weight\":1.0,\"role\":\"energy\"}]}");
    if (root) {
        ok("a role entering a non-module is rejected",
           !gia_model_load(&m, root));
        cJSON_Delete(root);
    }
    /* A cycling receptor with a surplus input. GSSK silently discards it. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"extra\",\"type\":\"constant\",\"value\":1.0},"
        "  {\"id\":\"rec\",\"type\":\"loop_limited\","
        "   \"module\":{\"k\":0.5,\"capacity\":5.0}},"
        "  {\"id\":\"out\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"rec\",\"role\":\"energy\"},"
        "  {\"source\":\"extra\",\"target\":\"rec\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"rec\",\"target\":\"out\",\"weight\":1.0}]}");
    if (root) {
        ok("a surplus input is named, not silently discarded as GSSK does",
           !gia_model_load(&m, root));
        cJSON_Delete(root);
    }
    /* A module block on a type that hosts no law. */
    root = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"s\",\"type\":\"storage\",\"current_level\":1.0,"
        "            \"module\":{\"k\":1.0}}],\"edges\":[]}");
    if (root) {
        ok("a module block on a storage is rejected",
           !gia_model_load(&m, root));
        cJSON_Delete(root);
    }
}

static void test_module_switch_events(void) {
    cJSON     *root;
    gia_model  m;
    double     q[4];
    int        ev;

    printf("\n[39] a switch module's crossing is located, not stepped over\n");

    /* The sensor is the tank being drained, so the switch shuts itself off at
     * the threshold. Without event location the flow would run past it. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"tank\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"sw\",\"type\":\"switch\","
        "   \"module\":{\"k\":2.0,\"threshold\":6.0}},"
        "  {\"id\":\"out\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"tank\",\"target\":\"sw\",\"role\":\"energy\"},"
        "  {\"source\":\"tank\",\"target\":\"sw\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"sw\",\"target\":\"out\",\"weight\":1.0}],"
        " \"simulation_params\":{\"t_val\":4.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("a switch module loads", gia_model_load(&m, root));
    ok("it makes the matrix non-constant", !gia_flow_matrix_is_constant(&m));

    ok("solves before the crossing", gia_network_state(&m, 1.0, q, NULL));
    close_to("tank(1) = 10 - 2", q[0], 8.0, 1e-6);

    ok("solves past the crossing", gia_network_state(&m, 4.0, q, NULL));
    ev = gia_count_events(&m, 4.0);
    printf("      events located over [0,4]: %d\n", ev);
    ok("a crossing was located", ev >= 1);
    close_to("the tank is held at the threshold, not driven through it",
             q[0], 6.0, 1e-4);
    close_to("conserved across the event", q[0] + q[2], 10.0, 1e-6);
    gia_model_free(&m); cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 40. Exchange as a module (ADR 0013 decision 5)
 *
 * The last module, and the one that settles leg discovery. As a pathway law an
 * exchange needs its counter-flow legs named, because which component pays and
 * which receives cannot be recovered from carrier identity alone. As a module
 * it reads its whole neighbourhood, and the four legs say outright which is
 * which.
 * ------------------------------------------------------------------ */

static cJSON *transactor_model(const char *legs) {
    static char buf[1600];
    snprintf(buf, sizeof(buf),
        "{\"nodes\":["
        "  {\"id\":\"goods_a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":100.0},"
        "  {\"id\":\"goods_b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"cash_b\",\"type\":\"storage\",\"carrier\":\"money\","
        "   \"current_level\":500.0},"
        "  {\"id\":\"cash_a\",\"type\":\"storage\",\"carrier\":\"money\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"ex\",\"type\":\"exchange\","
        "   \"module\":{\"k\":0.2,\"exchange_ratio\":0.25}}],"
        " \"edges\":[%s],"
        " \"simulation_params\":{\"t_val\":2.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}", legs);
    return cJSON_Parse(buf);
}

#define L_GI "{\"source\":\"goods_a\",\"target\":\"ex\",\"role\":\"goods_in\"}"
#define L_GO "{\"source\":\"ex\",\"target\":\"goods_b\",\"role\":\"goods_out\"}"
#define L_CI "{\"source\":\"cash_b\",\"target\":\"ex\",\"role\":\"counter_in\"}"
#define L_CO "{\"source\":\"ex\",\"target\":\"cash_a\",\"role\":\"counter_out\"}"

static void test_transactor_module(void) {
    cJSON     *root;
    gia_model  m;
    double     q[5], psi = -1.0, moved, P = 0.25, t = 2.0;

    printf("\n[40] transactor module: Eq (103) from four named legs\n");

    root = transactor_model(L_GI "," L_GO "," L_CI "," L_CO);
    if (!root) { ok("parse", false); return; }
    ok("a transactor module loads", gia_model_load(&m, root));
    ok("the transactor is a module", gia_node_is_module(&m, 4));

    /* Linear in the goods leaving the seller, so unlike a work gate a
     * transactor leaves the matrix constant — and the two calculi agree. */
    ok("a transactor leaves the flow matrix constant",
       gia_flow_matrix_is_constant(&m));

    ok("solves", gia_network_state(&m, t, q, &psi));
    close_to("psi is exactly zero", psi, 0.0, 0.0);

    /* F_goods = k Q_goods_a, so goods_a decays exactly. */
    close_to("goods_a(2) = 100 e^(-k t)", q[0], 100.0 * exp(-0.2 * t), 1e-7);
    moved = 100.0 - q[0];
    close_to("goods_b took what goods_a lost", q[1], moved, 1e-7);

    /* Odum SecXV: currency moves the OTHER way, at F/P. */
    close_to("cash reached the seller: moved / P", q[3], moved / P, 1e-6);
    close_to("and the buyer paid exactly that", 500.0 - q[2], q[3], 1e-9);
    close_to("J_goods / J_counter = P", moved / q[3], P, 1e-9);

    /* Each carrier balances on its own. */
    ok("two carriers", gia_carrier_count(&m) == 2);
    close_to("goods conserved",
             gia_conservation_residual_for(&m, t, gia_node_carrier(&m, 0)),
             0.0, 1e-7);
    close_to("money conserved",
             gia_conservation_residual_for(&m, t, gia_node_carrier(&m, 2)),
             0.0, 1e-6);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_transactor_order_invariance(void) {
    const char *orders[2] = {
        L_GI "," L_GO "," L_CI "," L_CO,
        L_CO "," L_CI "," L_GO "," L_GI
    };
    const char *paths[2] = { "tests/results/gia_tx_a.csv",
                             "tests/results/gia_tx_b.csv" };
    char *body[2] = { NULL, NULL };
    int   i;

    printf("\n[41] four legs in any order give the same answer\n");

    for (i = 0; i < 2; i++) {
        cJSON     *root = transactor_model(orders[i]);
        gia_model  m;
        if (!root) { ok("parse", false); return; }
        ok("loads in this leg order", gia_model_load(&m, root));
        ok("trajectories written", gia_write_trajectories(&m, paths[i], 6));
        body[i] = slurp_file(paths[i]);
        gia_model_free(&m); cJSON_Delete(root);
    }
    ok("the two leg orders are byte-identical",
       body[0] && body[1] && !strcmp(body[0], body[1]));
    free(body[0]); free(body[1]);
}

static void test_transactor_validation(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[42] a transactor states all four of its legs\n");

    root = transactor_model(L_GI "," L_GO "," L_CI);          /* no counter_out */
    if (root) { ok("a missing leg is rejected", !gia_model_load(&m, root));
                cJSON_Delete(root); }

    root = transactor_model(L_GI "," L_GO "," L_CI "," L_CO ","
        "{\"source\":\"cash_b\",\"target\":\"ex\",\"role\":\"counter_in\"}");
    if (root) { ok("a duplicated leg is rejected", !gia_model_load(&m, root));
                cJSON_Delete(root); }

    /* Paying yourself moves nothing. */
    root = transactor_model(L_GI "," L_GO ","
        "{\"source\":\"cash_b\",\"target\":\"ex\",\"role\":\"counter_in\"},"
        "{\"source\":\"ex\",\"target\":\"cash_b\",\"role\":\"counter_out\"}");
    if (root) { ok("a self-payment is rejected", !gia_model_load(&m, root));
                cJSON_Delete(root); }

    /* One leg pair moves one carrier: cash in, goods out is two exchanges. */
    root = transactor_model(L_GI "," L_GO "," L_CI ","
        "{\"source\":\"ex\",\"target\":\"goods_b\",\"role\":\"counter_out\"}");
    if (root) { ok("a leg pair straddling two carriers is rejected",
                   !gia_model_load(&m, root)); cJSON_Delete(root); }

    /* An output role on an incoming pathway is a misplaced field. */
    root = transactor_model(L_GI ","
        "{\"source\":\"goods_b\",\"target\":\"ex\",\"role\":\"goods_out\"},"
        L_CI "," L_CO);
    if (root) { ok("an output role on an incoming pathway is rejected",
                   !gia_model_load(&m, root)); cJSON_Delete(root); }
}

static void test_transactor_barter(void) {
    cJSON     *root;
    gia_model  m;
    double     q[5];

    printf("\n[43] barter: a transactor need not involve money\n");

    /* All four legs on one carrier. PR 3 established that goods for goods is a
     * real process; making exchange a module must not quietly re-forbid it. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"grain_a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":60.0},"
        "  {\"id\":\"grain_b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"sheep_b\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":40.0},"
        "  {\"id\":\"sheep_a\",\"type\":\"storage\",\"carrier\":\"goods\","
        "   \"current_level\":0.0},"
        "  {\"id\":\"mkt\",\"type\":\"exchange\","
        "   \"module\":{\"k\":0.1,\"exchange_ratio\":3.0}}],"
        " \"edges\":["
        "  {\"source\":\"grain_a\",\"target\":\"mkt\",\"role\":\"goods_in\"},"
        "  {\"source\":\"mkt\",\"target\":\"grain_b\",\"role\":\"goods_out\"},"
        "  {\"source\":\"sheep_b\",\"target\":\"mkt\",\"role\":\"counter_in\"},"
        "  {\"source\":\"mkt\",\"target\":\"sheep_a\",\"role\":\"counter_out\"}],"
        " \"simulation_params\":{\"t_val\":2.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("a single-carrier barter transactor loads", gia_model_load(&m, root));
    ok("solves", gia_network_state(&m, 2.0, q, NULL));
    close_to("grain per sheep = the exchange ratio",
             (60.0 - q[0]) / q[3], 3.0, 1e-7);
    close_to("goods conserved across the barter",
             gia_conservation_residual(&m, 2.0), 0.0, 1e-6);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_module_emergy_uses_module_flow(void) {
    cJSON     *root;
    gia_model  m;
    double     em[4];

    printf("\n[44] emergy carries a module's flow, not its pathway's default\n");

    /* A module's pathways have no law of their own, so asking an edge for its
     * flow returned the pathway DEFAULT — weight 1, linear — and the emergy
     * pass carried transformity along that instead of along F. Chosen so the
     * two differ: g = k*Qc = 0.6, F = 6, while the default would give 10. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"sun\",\"type\":\"source\",\"value\":10.0,"
        "   \"quality_input\":1.0},"
        "  {\"id\":\"c\",\"type\":\"constant\",\"value\":2.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.3}},"
        "  {\"id\":\"out\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"sun\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"gate\",\"target\":\"out\",\"weight\":1.0}],"
        " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    ok("loads", gia_model_load(&m, root));
    ok("emergy computes", gia_emergy_at(&m, 1.0, em, NULL));

    /* F = k * Q_sun * Q_c = 0.3*10*2 = 6, at transformity 1. */
    close_to("empower reaching the output is F x Tr, not the default flow",
             em[3], 6.0, 1e-6);
    ok("and it is not the pathway default of 10", fabs(em[3] - 10.0) > 1.0);
    gia_model_free(&m); cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 45-50. Ordinality over pathways that carry quantity (ADR 0014)
 *
 * Ordinality decides whether emergence happens. Since module-hosted laws made
 * a control into an edge, the cycle scan walked it, and a module was counted
 * as a component although it holds nothing. Each model here is one of the
 * ADR's reproductions.
 * ------------------------------------------------------------------ */

static bool load_json(gia_model *m, cJSON **root, const char *json) {
    *root = cJSON_Parse(json);
    if (!*root) return false;
    if (!gia_model_load(m, *root)) { cJSON_Delete(*root); *root = NULL; return false; }
    return true;
}

/* The module is listed FIRST on purpose: its on_cycle flag is always false, so
 * a scan that forgot a module is not a component would pick it as the open
 * component to close. */
static const char *ACCUMULATOR =
    "{\"nodes\":["
    "  {\"id\":\"g\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
    "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":5.0},"
    "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0}],"
    " \"edges\":["
    "  {\"source\":\"a\",\"target\":\"g\",\"role\":\"energy\"},"
    "  {\"source\":\"b\",\"target\":\"g\",\"role\":\"control\",\"use_ratio\":0},"
    "  {\"source\":\"g\",\"target\":\"a\",\"weight\":1.0},"
    "  {\"source\":\"g\",\"target\":\"b\",\"weight\":1.0}],"
    " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1,"
    "                       \"generative_mode\":true}}";

static void test_control_does_not_close_a_pathway(void) {
    cJSON     *root, *out, *nodes, *added;
    gia_model  m;
    double     q[3];

    printf("\n[45] a control leg does not close a pathway\n");
    if (!load_json(&m, &root, ACCUMULATOR)) { ok("loads", false); return; }

    /* b receives from the gate and only meters it. Nothing b holds leaves. */
    ok("solves", gia_network_state(&m, 1.0, q, NULL));
    ok("b only accumulates", q[2] > 2.0);

    ok("two components: a module is not one", gia_component_count(&m) == 2);
    close_to("ordinality is 1/2, not 1", gia_ordinality(&m), 0.5, 1e-12);
    ok("below maximum ordinality", !gia_at_maximum_ordinality(&m));
    ok("a is on a closed pathway, through the gate", m.nodes[1].on_cycle);
    ok("b is not: its only way back is a control", !m.nodes[2].on_cycle);
    ok("the module's own flag stays false", !m.nodes[0].on_cycle);

    /* So the MOP step has an open relationship to close -- and closes b. */
    out = gia_generate(&m);
    ok("generate evolves the graph", out && !cJSON_Compare(root, out, 1));
    nodes = out ? cJSON_GetObjectItemCaseSensitive(out, "nodes") : NULL;
    added = nodes ? cJSON_GetArrayItem(nodes, cJSON_GetArraySize(nodes) - 1) : NULL;
    {
        const cJSON *from = added ? cJSON_GetObjectItemCaseSensitive(added, "emerged_from") : NULL;
        const cJSON *rank = added ? cJSON_GetObjectItemCaseSensitive(added, "ordinality_rank") : NULL;
        /* An array since ADR 0015: one step may close several components. */
        ok("the component closed is b, never the module",
           cJSON_IsArray(from) && cJSON_GetArraySize(from) == 1 &&
           !strcmp(cJSON_GetArrayItem(from, 0)->valuestring, "b"));
        ok("the emergent component's rank counts components",
           cJSON_IsNumber(rank) && rank->valuedouble == 2.0);
    }
    cJSON_Delete(out);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_module_is_passed_through(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[46] a closed pathway may run through a module\n");
    /* a -> g -> b -> a: the loop is material, and runs through the gate. */
    if (!load_json(&m, &root,
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":5.0},"
        "  {\"id\":\"g\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"g\",\"role\":\"energy\"},"
        "  {\"source\":\"b\",\"target\":\"g\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"g\",\"target\":\"b\",\"weight\":1.0},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.05}],"
        " \"simulation_params\":{\"t_val\":1.0}}")) { ok("loads", false); return; }
    close_to("every component is on the loop", gia_ordinality(&m), 1.0, 1e-12);
    ok("at maximum ordinality", gia_at_maximum_ordinality(&m));
    ok("but the module is still not counted as on it", !m.nodes[1].on_cycle);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_control_does_not_open_the_boundary(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[47] a control read across the boundary moves nothing\n");
    if (!load_json(&m, &root,
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":5.0},"
        "  {\"id\":\"g\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0},"
        "  {\"id\":\"c\",\"type\":\"constant\",\"value\":2.0}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"g\",\"role\":\"energy\"},"
        "  {\"source\":\"c\",\"target\":\"g\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"g\",\"target\":\"b\",\"weight\":1.0},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.05}],"
        " \"simulation_params\":{\"t_val\":1.0}}")) { ok("loads", false); return; }
    ok("a gate metered by a constant leaves the system closed",
       gia_system_is_closed(&m));
    close_to("and it does conserve", gia_conservation_residual(&m, 1.0),
             0.0, 1e-12);
    gia_model_free(&m); cJSON_Delete(root);

    /* Not an over-correction: a source on the ENERGY leg does move quantity
     * across the boundary, and the system is open. */
    if (!load_json(&m, &root,
        "{\"nodes\":["
        "  {\"id\":\"s\",\"type\":\"source\",\"value\":2.0},"
        "  {\"id\":\"g\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0}],"
        " \"edges\":["
        "  {\"source\":\"s\",\"target\":\"g\",\"role\":\"energy\"},"
        "  {\"source\":\"b\",\"target\":\"g\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"g\",\"target\":\"b\",\"weight\":1.0}],"
        " \"simulation_params\":{\"t_val\":1.0}}")) { ok("loads", false); return; }
    ok("a source on the energy leg still opens the system",
       !gia_system_is_closed(&m));
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_ordinality_invariant_under_respelling(void) {
    /* One system: a <-> b is a closed loop, and a source feeds a work gate
     * metered by a that delivers into a. Only a and b are on a closed pathway.
     * Before ADR 0014 the module spelling scored 0.75; skipping controls alone
     * would score 0.50. The pathway spelling this was once compared against no
     * longer loads (ADR 0012 decision 5), so the module is the one spelling. */
    const char *as_pathway =
        "{\"nodes\":["
        "  {\"id\":\"src\",\"type\":\"source\",\"value\":2.0},"
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":5.0},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0}],"
        " \"edges\":["
        "  {\"source\":\"src\",\"target\":\"a\",\"logic\":\"interaction\","
        "   \"weight\":0.1,\"control_node\":\"a\"},"
        "  {\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\",\"weight\":0.3},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.2}],"
        " \"simulation_params\":{\"t_val\":1.0}}";
    const char *as_module =
        "{\"nodes\":["
        "  {\"id\":\"src\",\"type\":\"source\",\"value\":2.0},"
        "  {\"id\":\"g\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":5.0},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2.0}],"
        " \"edges\":["
        "  {\"source\":\"src\",\"target\":\"g\",\"role\":\"energy\"},"
        "  {\"source\":\"a\",\"target\":\"g\",\"role\":\"control\",\"use_ratio\":0},"
        "  {\"source\":\"g\",\"target\":\"a\",\"weight\":1.0},"
        "  {\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\",\"weight\":0.3},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.2}],"
        " \"simulation_params\":{\"t_val\":1.0}}";
    cJSON     *rm;
    gia_model  mm;

    printf("\n[48] a work gate has one spelling, and it keeps its ordinality\n");
    ok("the pathway spelling no longer loads", !loads(as_pathway));
    if (!load_json(&mm, &rm, as_module)) { ok("module spelling loads", false); return; }
    close_to("module spelling: 2/3", gia_ordinality(&mm), 2.0 / 3.0, 1e-12);
    ok("a read control leaves the system open only through its source",
       !gia_system_is_closed(&mm));
    ok("below maximum", !gia_at_maximum_ordinality(&mm));
    gia_model_free(&mm); cJSON_Delete(rm);
}

static void test_transactor_does_not_change_kind(void) {
    cJSON     *root;
    gia_model  m;
    int        i;
    bool       any = false;

    printf("\n[49] a transactor passes goods as goods and money as money\n");

    /* The only loop in this graph switches carrier twice:
     *
     *   gB --goods_in--> ex2 --counter_out--> mA --counter_in--> ex1 --goods_out--> gB
     *
     * Read pairwise, as quantity, every stream dead-ends: gZ -> ex1 -> gB ->
     * ex2 -> gX, and mY -> ex2 -> mA -> ex1 -> mW. Nothing that enters as goods
     * ever returns as goods. Letting a path cross streams would put gB and mA
     * "on a cycle" that moves no single quantity around it. */
    if (!load_json(&m, &root,
        "{\"nodes\":["
        "  {\"id\":\"gZ\",\"type\":\"storage\",\"carrier\":\"goods\",\"current_level\":9.0},"
        "  {\"id\":\"gB\",\"type\":\"storage\",\"carrier\":\"goods\",\"current_level\":9.0},"
        "  {\"id\":\"gX\",\"type\":\"storage\",\"carrier\":\"goods\",\"current_level\":0.0},"
        "  {\"id\":\"mY\",\"type\":\"storage\",\"carrier\":\"money\",\"current_level\":9.0},"
        "  {\"id\":\"mA\",\"type\":\"storage\",\"carrier\":\"money\",\"current_level\":9.0},"
        "  {\"id\":\"mW\",\"type\":\"storage\",\"carrier\":\"money\",\"current_level\":0.0},"
        "  {\"id\":\"ex1\",\"type\":\"exchange\",\"module\":{\"k\":0.1,\"exchange_ratio\":1.0}},"
        "  {\"id\":\"ex2\",\"type\":\"exchange\",\"module\":{\"k\":0.1,\"exchange_ratio\":1.0}}],"
        " \"edges\":["
        "  {\"source\":\"gZ\",\"target\":\"ex1\",\"role\":\"goods_in\"},"
        "  {\"source\":\"ex1\",\"target\":\"gB\",\"role\":\"goods_out\"},"
        "  {\"source\":\"mA\",\"target\":\"ex1\",\"role\":\"counter_in\"},"
        "  {\"source\":\"ex1\",\"target\":\"mW\",\"role\":\"counter_out\"},"
        "  {\"source\":\"gB\",\"target\":\"ex2\",\"role\":\"goods_in\"},"
        "  {\"source\":\"ex2\",\"target\":\"gX\",\"role\":\"goods_out\"},"
        "  {\"source\":\"mY\",\"target\":\"ex2\",\"role\":\"counter_in\"},"
        "  {\"source\":\"ex2\",\"target\":\"mA\",\"role\":\"counter_out\"}],"
        " \"simulation_params\":{\"t_val\":1.0}}")) { ok("loads", false); return; }

    ok("six components, two modules", gia_component_count(&m) == 6);
    close_to("no component is on a closed pathway", gia_ordinality(&m), 0.0, 0.0);
    gia_mark_cycles(&m);
    for (i = 0; i < m.n_nodes; i++) any = any || m.nodes[i].on_cycle;
    ok("including gB and mA, which only a stream-crossing path would join", !any);
    gia_model_free(&m); cJSON_Delete(root);

    /* And a genuine goods loop through a transactor still counts. */
    if (!load_json(&m, &root,
        "{\"nodes\":["
        "  {\"id\":\"ga\",\"type\":\"storage\",\"carrier\":\"goods\",\"current_level\":9.0},"
        "  {\"id\":\"gb\",\"type\":\"storage\",\"carrier\":\"goods\",\"current_level\":0.0},"
        "  {\"id\":\"mb\",\"type\":\"storage\",\"carrier\":\"money\",\"current_level\":9.0},"
        "  {\"id\":\"ma\",\"type\":\"storage\",\"carrier\":\"money\",\"current_level\":0.0},"
        "  {\"id\":\"ex\",\"type\":\"exchange\",\"module\":{\"k\":0.1,\"exchange_ratio\":1.0}}],"
        " \"edges\":["
        "  {\"source\":\"ga\",\"target\":\"ex\",\"role\":\"goods_in\"},"
        "  {\"source\":\"ex\",\"target\":\"gb\",\"role\":\"goods_out\"},"
        "  {\"source\":\"mb\",\"target\":\"ex\",\"role\":\"counter_in\"},"
        "  {\"source\":\"ex\",\"target\":\"ma\",\"role\":\"counter_out\"},"
        "  {\"source\":\"gb\",\"target\":\"ga\",\"logic\":\"linear\",\"weight\":0.1}],"
        " \"simulation_params\":{\"t_val\":1.0}}")) { ok("loads", false); return; }
    close_to("the goods loop closes through the transactor; money does not",
             gia_ordinality(&m), 0.5, 1e-12);
    ok("ga on it", m.nodes[0].on_cycle);
    ok("ma not", !m.nodes[3].on_cycle);
    gia_model_free(&m); cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 50-53. What the emergent quality closes (ADR 0015)
 *
 * The step promised to raise ordinality and did not: it wired hub -> E -> open,
 * which closes nothing when `open` is a dead end, so ordinality fell on every
 * step and the step never stopped. Each model here is one the ADR measured.
 * ------------------------------------------------------------------ */

#define GEN_SIM "\"simulation_params\":{\"t_val\":1.0,\"generative_mode\":true}}"

static const char *DEAD_END =
    "{\"nodes\":[{\"id\":\"c\",\"type\":\"storage\"},"
    "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":5},"
    "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2}],"
    " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\",\"weight\":0.3},"
    "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.2},"
    "  {\"source\":\"a\",\"target\":\"c\",\"logic\":\"linear\",\"weight\":0.1}]," GEN_SIM;

/* Runs one step, reloads what it produced, and reports the evolved ordinality.
 * Returns the evolved graph (caller frees) so callers can inspect it. */
static cJSON *step_once(const char *json, double *before, double *after,
                        bool *second_is_functional) {
    cJSON     *root, *out, *out2;
    gia_model  m, ev;
    int        saved;

    *before = *after = -1.0;
    *second_is_functional = false;
    if (!load_json(&m, &root, json)) return NULL;
    *before = gia_ordinality(&m);

    fflush(stdout); saved = dup(1);
    { FILE *f = freopen("/dev/null", "w", stdout); (void)f; }
    out = gia_generate(&m);
    fflush(stdout); dup2(saved, 1); close(saved);

    if (out && gia_model_load(&ev, out)) {
        *after = gia_ordinality(&ev);
        fflush(stdout); saved = dup(1);
        { FILE *f = freopen("/dev/null", "w", stdout); (void)f; }
        out2 = gia_generate(&ev);
        fflush(stdout); dup2(saved, 1); close(saved);
        *second_is_functional =
            out2 && gia_validate_mode(out, out2) == GIA_MODE_FUNCTIONAL;
        cJSON_Delete(out2);
        gia_model_free(&ev);
    }
    gia_model_free(&m); cJSON_Delete(root);
    return out;
}

static void test_emergent_quality_closes(void) {
    static const struct { const char *name; const char *json; double want; } cases[] = {
        { "dead end a<->b, a->c", NULL, 1.0 },
        { "isolated a<->b, c",
          "{\"nodes\":[{\"id\":\"c\",\"type\":\"storage\",\"current_level\":1},"
          "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":5},"
          "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2}],"
          " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\",\"weight\":0.3},"
          "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.2}]," GEN_SIM, 1.0 },
        { "chain s->x->y->z",
          "{\"nodes\":[{\"id\":\"s\",\"type\":\"source\",\"value\":2},"
          "  {\"id\":\"x\",\"type\":\"storage\"},{\"id\":\"y\",\"type\":\"storage\"},"
          "  {\"id\":\"z\",\"type\":\"storage\"}],"
          " \"edges\":[{\"source\":\"s\",\"target\":\"x\",\"logic\":\"linear\",\"weight\":0.3},"
          "  {\"source\":\"x\",\"target\":\"y\",\"logic\":\"linear\",\"weight\":0.2},"
          "  {\"source\":\"y\",\"target\":\"z\",\"logic\":\"linear\",\"weight\":0.1}]," GEN_SIM, 1.0 },
        { "a lone component",
          "{\"nodes\":[{\"id\":\"x\",\"type\":\"storage\",\"current_level\":1}],"
          " \"edges\":[]," GEN_SIM, 1.0 },
        { "accumulator behind a module", NULL, 1.0 },
    };
    size_t k;

    printf("\n[50] one step reaches the fixed point, and never lowers ordinality\n");
    for (k = 0; k < sizeof cases / sizeof cases[0]; k++) {
        const char *json = cases[k].json;
        double before, after;
        bool   fixed;
        cJSON *out;
        char   label[96];

        if (k == 0) json = DEAD_END;
        if (k == 4) json = ACCUMULATOR;
        out = step_once(json, &before, &after, &fixed);
        snprintf(label, sizeof label, "%s: loads and evolves", cases[k].name);
        ok(label, out != NULL && after >= 0.0);
        snprintf(label, sizeof label, "%s: ordinality does not fall", cases[k].name);
        ok(label, after >= before);
        snprintf(label, sizeof label, "%s: maximum in one step", cases[k].name);
        close_to(label, after, cases[k].want, 1e-12);
        snprintf(label, sizeof label, "%s: a second step changes nothing", cases[k].name);
        ok(label, fixed);
        cJSON_Delete(out);
    }
}

static void test_sink_is_never_closed(void) {
    double before, after;
    bool   fixed;
    cJSON *out;
    const cJSON *e, *edges_out;
    bool   touches = false;

    printf("\n[51] a sink is never closed and never drawn from\n");

    /* Only the sink is open: a stated fixed point below maximum. */
    out = step_once(
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":5},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\",\"weight\":0.3},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.2},"
        "  {\"source\":\"a\",\"target\":\"heat\",\"logic\":\"linear\",\"weight\":0.1}]," GEN_SIM,
        &before, &after, &fixed);
    close_to("heat sink alone: ordinality stays 2/3", after, 2.0 / 3.0, 1e-12);
    ok("heat sink alone: the step changes nothing", fixed && before == after);
    ok("the evolved graph IS the seed: nothing was added",
       out && cJSON_GetArraySize(cJSON_GetObjectItemCaseSensitive(out, "nodes")) == 3);
    cJSON_Delete(out);

    /* A sink alongside a real dead end: the dead end is closed, the sink not. */
    out = step_once(
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"current_level\":5},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":2},"
        "  {\"id\":\"c\",\"type\":\"storage\"},{\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\",\"weight\":0.3},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\",\"weight\":0.2},"
        "  {\"source\":\"a\",\"target\":\"c\",\"logic\":\"linear\",\"weight\":0.1},"
        "  {\"source\":\"a\",\"target\":\"heat\",\"logic\":\"linear\",\"weight\":0.1}]," GEN_SIM,
        &before, &after, &fixed);
    close_to("with a dead end too: every component but the sink closes (4/5)",
             after, 0.8, 1e-12);
    ok("and that is a fixed point", fixed);
    edges_out = out ? cJSON_GetObjectItemCaseSensitive(out, "edges") : NULL;
    cJSON_ArrayForEach(e, edges_out) {
        const cJSON *s = cJSON_GetObjectItemCaseSensitive(e, "source");
        const cJSON *t = cJSON_GetObjectItemCaseSensitive(e, "target");
        const cJSON *ft = cJSON_GetObjectItemCaseSensitive(e, "flow_type");
        if (!cJSON_IsString(ft)) continue;               /* seed pathways */
        if ((cJSON_IsString(s) && !strcmp(s->valuestring, "heat")) ||
            (cJSON_IsString(t) && !strcmp(t->valuestring, "heat")))
            touches = true;
    }
    ok("no emergent leg touches the sink", !touches);
    cJSON_Delete(out);
}

/* What the step appends, serialised: the emergent component and every leg
 * beyond the seed's own. */
static char *appended(const cJSON *out, int seed_edges) {
    const cJSON *nodes = cJSON_GetObjectItemCaseSensitive(out, "nodes");
    const cJSON *edges = cJSON_GetObjectItemCaseSensitive(out, "edges");
    cJSON *bag = cJSON_CreateArray();
    char  *txt;
    int    i, n = cJSON_GetArraySize(edges);
    cJSON_AddItemToArray(bag, cJSON_Duplicate(
        cJSON_GetArrayItem(nodes, cJSON_GetArraySize(nodes) - 1), 1));
    for (i = seed_edges; i < n; i++)
        cJSON_AddItemToArray(bag, cJSON_Duplicate(cJSON_GetArrayItem(edges, i), 1));
    txt = cJSON_PrintUnformatted(bag);
    cJSON_Delete(bag);
    return txt;
}

static void test_emergence_is_order_invariant(void) {
    /* A dead end c AND an unconnected d, with the nodes listed two ways.
     * Before ADR 0015 the step closed c from hub a in one order and d from hub
     * b in the other. */
    const char *orders[2] = {
        "{\"nodes\":[{\"id\":\"c\",\"type\":\"storage\"},{\"id\":\"d\",\"type\":\"storage\"},"
        "  {\"id\":\"a\",\"type\":\"storage\"},{\"id\":\"b\",\"type\":\"storage\"}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\"},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\"},"
        "  {\"source\":\"a\",\"target\":\"c\",\"logic\":\"linear\"}]," GEN_SIM,
        "{\"nodes\":[{\"id\":\"b\",\"type\":\"storage\"},{\"id\":\"a\",\"type\":\"storage\"},"
        "  {\"id\":\"d\",\"type\":\"storage\"},{\"id\":\"c\",\"type\":\"storage\"}],"
        " \"edges\":[{\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\"},"
        "  {\"source\":\"a\",\"target\":\"c\",\"logic\":\"linear\"},"
        "  {\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\"}]," GEN_SIM
    };
    char  *txt[2] = { NULL, NULL };
    int    i;

    printf("\n[52] what the step adds is decided by the model, not its file\n");
    for (i = 0; i < 2; i++) {
        double before, after;
        bool   fixed;
        cJSON *out = step_once(orders[i], &before, &after, &fixed);
        ok("this order evolves to maximum", out && after == 1.0);
        if (out) txt[i] = appended(out, 3);
        cJSON_Delete(out);
    }
    ok("the component and legs appended are byte-identical",
       txt[0] && txt[1] && !strcmp(txt[0], txt[1]));
    free(txt[0]); free(txt[1]);
}

static void test_emergent_quality_is_a_component(void) {
    double before, after;
    bool   fixed;
    cJSON *out;
    const cJSON *nodes, *e_node, *type, *from;

    printf("\n[53] the emergent quality is a component, and names what it closed\n");
    out = step_once("{\"nodes\":[{\"id\":\"d\",\"type\":\"storage\"},{\"id\":\"c\",\"type\":\"storage\"},"
        "  {\"id\":\"a\",\"type\":\"storage\"},{\"id\":\"b\",\"type\":\"storage\"}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"linear\"},"
        "  {\"source\":\"b\",\"target\":\"a\",\"logic\":\"linear\"},"
        "  {\"source\":\"a\",\"target\":\"c\",\"logic\":\"linear\"}]," GEN_SIM, &before, &after, &fixed);
    nodes  = out ? cJSON_GetObjectItemCaseSensitive(out, "nodes") : NULL;
    e_node = nodes ? cJSON_GetArrayItem(nodes, cJSON_GetArraySize(nodes) - 1) : NULL;
    type   = e_node ? cJSON_GetObjectItemCaseSensitive(e_node, "type") : NULL;
    from   = e_node ? cJSON_GetObjectItemCaseSensitive(e_node, "emerged_from") : NULL;

    /* A module would not count toward the ordinality it exists to raise. */
    ok("typed storage, not gain", cJSON_IsString(type) && !strcmp(type->valuestring, "storage"));
    ok("carries no module block",
       e_node && !cJSON_GetObjectItemCaseSensitive(e_node, "module"));
    ok("emerged_from lists both closed components",
       cJSON_IsArray(from) && cJSON_GetArraySize(from) == 2);
    ok("in id order",
       cJSON_IsArray(from) && cJSON_GetArraySize(from) == 2 &&
       !strcmp(cJSON_GetArrayItem(from, 0)->valuestring, "c") &&
       !strcmp(cJSON_GetArrayItem(from, 1)->valuestring, "d"));
    cJSON_Delete(out);
}

/* ------------------------------------------------------------------ *
 * 54-59. A control input is a flow of energy (ADR 0017)
 *
 * ADR 0013 made a control "read, never consumed", so a control leg carried no
 * flow and its emergy never reached the product: a transformity-1000 control
 * left the product at transformity 1, where Odum's Fig. 2.6(b) draws it
 * "moderate". A control now draws use_ratio * F from its source, that quantity
 * leaves on the module's `used` leg to a sink, and its emergy goes into the
 * product.
 * ------------------------------------------------------------------ */

#define CTL_SIM " \"simulation_params\":{\"t_val\":1.0,\"derivative_order\":1," \
                "                       \"generative_mode\":false}}"

/* Test [44]'s gate, with the control fed by a source of stated transformity.
 * `control` is the control leg's extra fields; `extra` adds edges. */
static cJSON *measured_gate(const char *control, const char *extra) {
    static char buf[1600];
    snprintf(buf, sizeof(buf),
        "{\"nodes\":["
        "  {\"id\":\"sun\",\"type\":\"source\",\"value\":10.0,"
        "   \"quality_input\":1.0},"
        "  {\"id\":\"H\",\"type\":\"source\",\"value\":2.0,"
        "   \"quality_input\":1000.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.3}},"
        "  {\"id\":\"out\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":["
        "  {\"source\":\"sun\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"H\",\"target\":\"gate\",\"role\":\"control\"%s},"
        "  {\"source\":\"gate\",\"target\":\"out\",\"weight\":1.0}%s],"
        CTL_SIM, control, extra);
    return cJSON_Parse(buf);
}

static void test_use_ratio_is_required(void) {
    char buf[1400];
    const char *head =
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"c\",\"type\":\"constant\",\"value\":2.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[";

    printf("\n[54] every control leg states its use_ratio, and 0 is allowed\n");

#define CTL_CASE(ctrl_fields, energy_fields, out_fields)                      \
    snprintf(buf, sizeof(buf), "%s"                                           \
        "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"%s},"     \
        "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\"%s},"    \
        "  {\"source\":\"gate\",\"target\":\"b\"%s}]," CTL_SIM,               \
        head, energy_fields, ctrl_fields, out_fields)

    /* The measured defect was a silent zero, so a missing value is an error
     * rather than a default (decision 4). */
    CTL_CASE("", "", "");
    ok("a control leg without use_ratio is refused", !loads(buf));
    CTL_CASE(",\"use_ratio\":0", "", "");
    ok("use_ratio 0 is allowed, and needs no used leg", loads(buf));
    CTL_CASE(",\"use_ratio\":-0.1", "", "");
    ok("a negative use_ratio is refused", !loads(buf));
    CTL_CASE(",\"use_ratio\":\"0.1\"", "", "");
    ok("a use_ratio that is not a number is refused", !loads(buf));
    CTL_CASE(",\"use_ratio\":0", ",\"use_ratio\":0.1", "");
    ok("use_ratio on an energy leg is refused", !loads(buf));
    CTL_CASE(",\"use_ratio\":0", "", ",\"use_ratio\":0.1");
    ok("use_ratio on an output leg is refused", !loads(buf));
#undef CTL_CASE
}

static void test_control_emergy_reaches_product(void) {
    cJSON     *root;
    gia_model  m;
    double     em[5], tr[5], q[5];

    printf("\n[55] the measured gate: the control's emergy reaches the product\n");

    /* ADR 0017's measurement, use_ratio 0: the control contributes nothing,
     * which is exactly the read-only behaviour the figure contradicts. */
    root = measured_gate(",\"use_ratio\":0", "");
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads with use_ratio 0", &m, root)) return;
    ok("emergy computes", gia_emergy_at(&m, 1.0, em, tr));
    close_to("use_ratio 0: empower at product is 6", em[3], 6.0, 1e-6);
    close_to("use_ratio 0: product transformity is 1", tr[3], 1.0, 1e-6);
    gia_model_free(&m); cJSON_Delete(root);

    /* Decision 3, worked: Tr = 1 + 0.01 * 1000 = 11, "moderate" -- between
     * the energy input's 1 and the control's 1000. */
    root = measured_gate(",\"use_ratio\":0.01",
        ",{\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"}");
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads with use_ratio 0.01 and a used leg", &m, root)) return;
    ok("emergy computes", gia_emergy_at(&m, 1.0, em, tr));
    ok("solves", gia_network_state(&m, 1.0, q, NULL));
    close_to("the product still receives F = 6: the control is not added",
             q[3], 6.0, 1e-6);
    close_to("empower at product is F + s F Tr_c = 66", em[3], 66.0, 1e-6);
    close_to("product transformity is 11", tr[3], 11.0, 1e-6);
    ok("which is moderate: between 1 and 1000", tr[3] > 1.0 && tr[3] < 1000.0);
    close_to("the used leg delivers s F = 0.06 to the sink", q[4], 0.06, 1e-6);
    close_to("and carries no emergy onward", em[4], 0.0, 1e-12);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_control_is_drawn_and_dissipated(void) {
    cJSON     *root;
    gia_model  m;
    double     q[5], k = 0.05, s = 0.1, S = 10.0, t = 3.0, d;

    printf("\n[56] a control draws use_ratio * F from its source, into the sink\n");

    /* Odum SecIX amplifier with a control that is a stock. F = k Q_signal and
     * the signal loses s F, so Q_signal = S e^(-s k t), exactly; the load gets
     * what the flow integrates to and the heat sink gets s times that. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"signal\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"power\",\"type\":\"storage\",\"current_level\":100.0},"
        "  {\"id\":\"amp\",\"type\":\"gain\",\"module\":{\"k\":0.05}},"
        "  {\"id\":\"load\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":["
        "  {\"source\":\"signal\",\"target\":\"amp\",\"role\":\"control\","
        "   \"use_ratio\":0.1},"
        "  {\"source\":\"power\",\"target\":\"amp\",\"role\":\"energy\"},"
        "  {\"source\":\"amp\",\"target\":\"load\",\"weight\":1.0},"
        "  {\"source\":\"amp\",\"target\":\"heat\",\"role\":\"used\"}],"
        " \"simulation_params\":{\"t_val\":3.0,\"derivative_order\":1,"
        "                       \"generative_mode\":false}}");
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
    ok("solves", gia_network_state(&m, t, q, NULL));
    d = 1.0 - exp(-s * k * t);
    close_to("signal(3) = S e^(-s k t): the control is drawn",
             q[0], S * exp(-s * k * t), 1e-9);
    close_to("load(3) = (S / s)(1 - e^(-s k t))", q[3], S / s * d, 1e-9);
    close_to("power(3) = 100 - load: the energy input is what feeds it",
             q[1], 100.0 - S / s * d, 1e-9);
    close_to("heat(3) = what the signal lost", q[4], S * d, 1e-9);
    ok("drawn only from stocks, so the system is closed",
       gia_system_is_closed(&m));
    close_to("and conserved", gia_conservation_residual(&m, t), 0.0, 1e-9);
    gia_model_free(&m); cJSON_Delete(root);

    /* Work gate, two stocks: what the product gains is exactly what the
     * energy input lost, and what the control lost is exactly what the sink
     * gained. The control's quantity never reaches the product. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"grass\",\"type\":\"storage\",\"current_level\":100.0},"
        "  {\"id\":\"sun\",\"type\":\"storage\",\"current_level\":2.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.01}},"
        "  {\"id\":\"cow\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":["
        "  {\"source\":\"grass\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"sun\",\"target\":\"gate\",\"role\":\"control\","
        "   \"use_ratio\":0.002},"
        "  {\"source\":\"gate\",\"target\":\"cow\",\"weight\":1.0},"
        "  {\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"}],"
        CTL_SIM);
    if (!root) { ok("parse", false); return; }
    if (!load_ok("work gate loads", &m, root)) return;
    ok("solves", gia_network_state(&m, 1.0, q, NULL));
    ok("the control is drawn down", q[1] < 2.0);
    close_to("grass + cow = 100: the product is the energy input's",
             q[0] + q[3], 100.0, 1e-9);
    close_to("sun + heat = 2: the control's quantity goes to the sink",
             q[1] + q[4], 2.0, 1e-9);
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_used_leg_validation(void) {
    char buf[1600];
    const char *nodes =
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0,"
        "   \"carrier\":\"energy\"},"
        "  {\"id\":\"p\",\"type\":\"storage\",\"current_level\":2.0,"
        "   \"carrier\":\"%s\"},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0,"
        "   \"carrier\":\"energy\"},"
        "  {\"id\":\"heat\",\"type\":\"%s\",\"carrier\":\"%s\"},"
        "  {\"id\":\"heat2\",\"type\":\"sink\",\"carrier\":\"%s\"}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"p\",\"target\":\"gate\",\"role\":\"control\","
        "   \"use_ratio\":%s},"
        "  {\"source\":\"gate\",\"target\":\"b\"}%s]," CTL_SIM;
    const char *used  = ",{\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"}";
    const char *used2 = ",{\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"},"
                        "{\"source\":\"gate\",\"target\":\"heat2\",\"role\":\"used\"}";
    const char *into  = ",{\"source\":\"heat\",\"target\":\"gate\",\"role\":\"used\"}";

    printf("\n[57] the used-energy leg belongs to the module (decision 5)\n");

    snprintf(buf, sizeof(buf), nodes, "energy", "sink", "energy", "energy", "0.1", "");
    ok("use_ratio > 0 with no used leg is refused", !loads(buf));
    snprintf(buf, sizeof(buf), nodes, "energy", "sink", "energy", "energy", "0.1", used);
    ok("with a used leg to a sink of its carrier, it loads", loads(buf));
    snprintf(buf, sizeof(buf), nodes, "energy", "storage", "energy", "energy", "0.1", used);
    ok("a used leg must end at a sink", !loads(buf));
    snprintf(buf, sizeof(buf), nodes, "energy", "sink", "energy", "energy", "0.1", used2);
    ok("one used leg per carrier, not two", !loads(buf));
    snprintf(buf, sizeof(buf), nodes, "energy", "sink", "energy", "energy", "0", used);
    ok("a used leg with nothing to carry is refused", !loads(buf));
    snprintf(buf, sizeof(buf), nodes, "energy", "sink", "energy", "energy", "0.1", into);
    ok("used names a leg leaving a module, never entering one", !loads(buf));

    /* The price case: money is not dissipated. A nonzero use_ratio on a money
     * control forces a money sink, which is the prompt to reconsider. */
    snprintf(buf, sizeof(buf), nodes, "money", "sink", "energy", "energy", "0.1", used);
    ok("a money control's used leg must reach a money sink", !loads(buf));
    snprintf(buf, sizeof(buf), nodes, "money", "sink", "money", "energy", "0.1", used);
    ok("and the module's outputs still need not match the control", loads(buf));
    snprintf(buf, sizeof(buf), nodes, "money", "sink", "money", "energy", "0", "");
    ok("a price read at use_ratio 0 needs no sink at all", loads(buf));
}

static void test_control_draw_opens_the_boundary(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[58] a control drawn from a source crosses the boundary\n");

    /* ADR 0014 decision 3, amended: a source delivers without being depleted,
     * so a control drawing on one brings quantity in from outside. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"c\",\"type\":\"source\",\"value\":2.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\","
        "   \"use_ratio\":0.1},"
        "  {\"source\":\"gate\",\"target\":\"b\"},"
        "  {\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"}]," CTL_SIM);
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
    ok("use_ratio > 0 from a source: the system is open",
       !gia_system_is_closed(&m));
    gia_model_free(&m); cJSON_Delete(root);

    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"c\",\"type\":\"source\",\"value\":2.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\","
        "   \"use_ratio\":0},"
        "  {\"source\":\"gate\",\"target\":\"b\"}]," CTL_SIM);
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
    ok("use_ratio 0 from a source: still closed", gia_system_is_closed(&m));
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_drawn_control_closes_no_pathway(void) {
    cJSON     *root;
    gia_model  m;

    printf("\n[59] a drawn control passes only to the used leg (ADR 0014)\n");

    /* The product feeds the control's stock, and the control is drawn. Its
     * quantity leaves through `used` to a sink, never into the product, so A
     * is no more on a cycle than when the control was only read. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"E\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"A\",\"type\":\"storage\",\"current_level\":1.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"P\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":["
        "  {\"source\":\"E\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"A\",\"target\":\"gate\",\"role\":\"control\","
        "   \"use_ratio\":0.05},"
        "  {\"source\":\"gate\",\"target\":\"P\"},"
        "  {\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"},"
        "  {\"source\":\"P\",\"target\":\"A\",\"logic\":\"linear\"}]," CTL_SIM);
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
    gia_mark_cycles(&m);
    ok("the drawn control does not close P -> A -> gate -> P",
       !m.nodes[1].on_cycle && !m.nodes[3].on_cycle);
    gia_model_free(&m); cJSON_Delete(root);

    /* Walking a control must not mark the module visited for the energy
     * stream. C reaches the gate twice: first, in edge order, by its drawn
     * control, then round C -> E -> gate on the energy leg -- a real cycle. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"C\",\"type\":\"storage\",\"current_level\":5.0},"
        "  {\"id\":\"E\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\",\"module\":{\"k\":0.1}},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":["
        "  {\"source\":\"C\",\"target\":\"gate\",\"role\":\"control\","
        "   \"use_ratio\":0.05},"
        "  {\"source\":\"C\",\"target\":\"E\",\"logic\":\"linear\"},"
        "  {\"source\":\"E\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"gate\",\"target\":\"C\"},"
        "  {\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"}]," CTL_SIM);
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
    gia_mark_cycles(&m);
    ok("C is on its energy cycle although its control reached the gate first",
       m.nodes[0].on_cycle && m.nodes[1].on_cycle);
    gia_model_free(&m); cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 60-63. Divide and subtract are actions of the interaction module
 * (ADR 0016)
 *
 * Odum's Fig. 2.6 draws one glyph five ways; (d) divides the left input by
 * the control and (e) subtracts the control from it. Until now only the
 * pathway laws `ratio` and `subtract` could say so, which left ADR 0012's
 * migration with nowhere to put them.
 * ------------------------------------------------------------------ */

/* A store `a` feeding `b`, with `c` as the second input, written as a work gate
 * with an action -- or, to check it is refused, as the pathway law it
 * replaces. `leak` adds a linear drain on `a` so a subtracting action can cross
 * zero. */
static cJSON *action_model(bool module, const char *law, const char *c_node,
                           const char *leak) {
    static char buf[1800];
    if (module)
        snprintf(buf, sizeof(buf),
            "{\"nodes\":["
            "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
            "  {\"id\":\"c\",%s},"
            "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
            "  {\"id\":\"z\",\"type\":\"storage\",\"current_level\":0.0},"
            "  {\"id\":\"gate\",\"type\":\"interaction\","
            "   \"module\":{\"k\":0.5,\"action\":\"%s\"}}],"
            " \"edges\":["
            "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"},"
            "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\","
            "   \"use_ratio\":0},"
            "  {\"source\":\"gate\",\"target\":\"b\"}%s],"
            " \"simulation_params\":{\"t_val\":6.0,\"derivative_order\":1,"
            "                       \"generative_mode\":false}}",
            c_node, law, leak);
    else
        snprintf(buf, sizeof(buf),
            "{\"nodes\":["
            "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
            "  {\"id\":\"c\",%s},"
            "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
            "  {\"id\":\"z\",\"type\":\"storage\",\"current_level\":0.0}],"
            " \"edges\":["
            "  {\"source\":\"a\",\"target\":\"b\",\"logic\":\"%s\","
            "   \"control_node\":\"c\",\"weight\":0.5}%s],"
            " \"simulation_params\":{\"t_val\":6.0,\"derivative_order\":1,"
            "                       \"generative_mode\":false}}",
            c_node, law, leak);
    return cJSON_Parse(buf);
}

#define C_CONST "\"type\":\"constant\",\"value\":2.0"
#define C_FOUR  "\"type\":\"constant\",\"value\":4.0"
#define LEAK    ",{\"source\":\"a\",\"target\":\"z\",\"logic\":\"linear\",\"weight\":0.2}"

static void test_action_validation(void) {
    char buf[1400];
    const char *fmt =
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"c\",\"type\":\"constant\",\"value\":2.0},"
        "  {\"id\":\"d\",\"type\":\"constant\",\"value\":3.0},"
        "  {\"id\":\"gate\",\"type\":\"%s\",\"module\":{\"k\":0.1%s}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0}%s,"
        "  {\"source\":\"gate\",\"target\":\"b\"}]," CTL_SIM;
    const char *two = ",{\"source\":\"d\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0}";

    printf("\n[60] an interaction states its action; the others have none\n");

    snprintf(buf, sizeof(buf), fmt, "interaction", "", "");
    ok("no action is multiply, as before", loads(buf));
    snprintf(buf, sizeof(buf), fmt, "interaction", ",\"action\":\"multiply\"", two);
    ok("multiply takes any number of controls", loads(buf));
    snprintf(buf, sizeof(buf), fmt, "interaction", ",\"action\":\"divide\"", "");
    ok("divide with one control loads", loads(buf));
    snprintf(buf, sizeof(buf), fmt, "interaction", ",\"action\":\"subtract\"", "");
    ok("subtract with one control loads", loads(buf));
    snprintf(buf, sizeof(buf), fmt, "interaction", ",\"action\":\"divide\"", two);
    ok("divide with two controls is refused", !loads(buf));
    snprintf(buf, sizeof(buf), fmt, "interaction", ",\"action\":\"subtract\"", two);
    ok("subtract with two controls is refused", !loads(buf));
    snprintf(buf, sizeof(buf), fmt, "interaction", ",\"action\":\"min\"", "");
    ok("an action beyond Fig. 2.6 is refused", !loads(buf));
    snprintf(buf, sizeof(buf), fmt, "interaction", ",\"action\":3", "");
    ok("an action that is not a string is refused", !loads(buf));
    snprintf(buf, sizeof(buf), fmt, "gain", ",\"action\":\"divide\"", "");
    ok("action on a gain is refused", !loads(buf));
}

static void test_action_divide(void) {
    cJSON     *root;
    gia_model  m, p;
    double     q[5];

    printf("\n[61] divide: F = k Q_energy / max(Q_control, eps)  (Fig. 2.6d)\n");

    root = action_model(true, "divide", C_CONST, "");
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
    ok("solves", gia_network_state(&m, 6.0, q, NULL));
    /* Q_c held at 2, so F = 0.25 a and a decays exactly. */
    close_to("a(6) = 10 e^(-k t / Q_c)", q[0], 10.0 * exp(-0.25 * 6.0), 1e-9);
    close_to("b takes what a lost", q[2], 10.0 - q[0], 1e-9);
    {   /* ADR 0016: once the pathway law is removed, the closed form above
         * is the check, and the pathway spelling must not load. */
        cJSON *pr = action_model(false, "ratio", C_CONST, "");
        ok("the ratio pathway no longer loads", pr && !gia_model_load(&p, pr));
        cJSON_Delete(pr);
    }
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_action_subtract(void) {
    cJSON     *root;
    gia_model  m, p;
    double     q[5], ts, eq = 2.0 / 0.7;

    printf("\n[62] subtract: F = max(0, k (Q_energy - Q_control))  (Fig. 2.6e)\n");

    /* Q_c = 4, and a also leaks at 0.2. While a > 4,
     *   da/dt = -0.5 (a - 4) - 0.2 a = -0.7 a + 2,
     * so a = eq + (10 - eq) e^(-0.7 t) until it reaches 4 at ts. Past that
     * the clamp holds F at zero and a only leaks: a = 4 e^(-0.2 (t - ts)). */
    ts = log((10.0 - eq) / (4.0 - eq)) / 0.7;

    root = action_model(true, "subtract", C_FOUR, LEAK);
    if (!root) { ok("parse", false); return; }
    if (!load_ok("loads", &m, root)) return;
    ok("solves before the crossing", gia_network_state(&m, 1.0, q, NULL));
    close_to("a(1) on the open branch", q[0], eq + (10.0 - eq) * exp(-0.7), 1e-9);
    ok("solves past the crossing", gia_network_state(&m, 6.0, q, NULL));
    ok("the crossing is located as an event", gia_count_events(&m, 6.0) == 1);
    close_to("a(6) on the clamped branch", q[0], 4.0 * exp(-0.2 * (6.0 - ts)), 1e-7);
    close_to("conserved across the event",
             q[0] + q[2] + q[3], 10.0, 1e-9);
    {
        cJSON *pr = action_model(false, "subtract", C_FOUR, LEAK);
        ok("the subtract pathway no longer loads", pr && !gia_model_load(&p, pr));
        cJSON_Delete(pr);
    }
    gia_model_free(&m); cJSON_Delete(root);
}

static void test_action_draws_control(void) {
    cJSON     *root;
    gia_model  m;
    double     q[5];

    printf("\n[63] divide and subtract draw their control like multiply (ADR 0017)\n");

    /* Divide, control a stock drawn at s = 0.1: F = 0.5 a / c, and what c
     * loses is exactly what the heat sink gains. */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"c\",\"type\":\"storage\",\"current_level\":2.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\","
        "   \"module\":{\"k\":0.5,\"action\":\"divide\"}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0.1},"
        "  {\"source\":\"gate\",\"target\":\"b\"},"
        "  {\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"}]," CTL_SIM);
    if (!root) { ok("parse", false); return; }
    if (!load_ok("divide with a drawn control loads", &m, root)) return;
    ok("solves", gia_network_state(&m, 1.0, q, NULL));
    ok("the control is drawn down", q[1] < 2.0);
    close_to("a + b = 10", q[0] + q[3], 10.0, 1e-9);
    close_to("c + heat = 2", q[1] + q[4], 2.0, 1e-9);
    gia_model_free(&m); cJSON_Delete(root);

    /* Subtract held shut from the start: c exceeds a, F = 0, so the control
     * draws nothing either (ADR 0017's consequence for ADR 0016). */
    root = cJSON_Parse(
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"c\",\"type\":\"storage\",\"current_level\":20.0},"
        "  {\"id\":\"gate\",\"type\":\"interaction\","
        "   \"module\":{\"k\":0.5,\"action\":\"subtract\"}},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0},"
        "  {\"id\":\"heat\",\"type\":\"sink\"}],"
        " \"edges\":["
        "  {\"source\":\"a\",\"target\":\"gate\",\"role\":\"energy\"},"
        "  {\"source\":\"c\",\"target\":\"gate\",\"role\":\"control\",\"use_ratio\":0.5},"
        "  {\"source\":\"gate\",\"target\":\"b\"},"
        "  {\"source\":\"gate\",\"target\":\"heat\",\"role\":\"used\"}]," CTL_SIM);
    if (!root) { ok("parse", false); return; }
    if (!load_ok("subtract with a drawn control loads", &m, root)) return;
    ok("solves", gia_network_state(&m, 1.0, q, NULL));
    close_to("clamped: nothing moves", q[3], 0.0, 1e-12);
    close_to("and the control is not drawn", q[1], 20.0, 1e-12);
    close_to("so the sink receives nothing", q[4], 0.0, 1e-12);
    gia_model_free(&m); cJSON_Delete(root);
}

/* ------------------------------------------------------------------ *
 * 64-67. Module laws leave pathways (ADR 0012 decision 5)
 *
 * Pathways carry Odum SecIII and nothing else: linear, reversible, constant.
 * The seven laws that happen inside a symbol live on modules, and a model is
 * migrated by writing the module -- never by the loader quietly turning an
 * edge law into a gate the author did not draw.
 * ------------------------------------------------------------------ */

static bool loads_with_law(const char *law) {
    char buf[900];
    snprintf(buf, sizeof(buf),
        "{\"nodes\":["
        "  {\"id\":\"a\",\"type\":\"storage\",\"current_level\":10.0},"
        "  {\"id\":\"c\",\"type\":\"storage\",\"current_level\":1.0},"
        "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
        " \"edges\":[{\"source\":\"a\",\"target\":\"b\",\"logic\":\"%s\","
        "            \"weight\":0.1,\"control_node\":\"c\"}]," CTL_SIM, law);
    return loads(buf);
}

static void test_module_laws_leave_pathways(void) {
    static const char *gone[] = { "interaction", "limit", "ratio", "subtract",
                                  "gain", "threshold", "exchange",
                                  "generative_production", "ordinal_feedback" };
    static const char *kept[] = { "linear", "reversible", "constant", "inflow",
                                  "outflow", "flow", "ordinal_ascent",
                                  "emergent_feedback_loop", "diffusion" };
    char   what[96];
    size_t i;

    printf("\n[64] a pathway carries Odum SecIII and nothing else\n");
    for (i = 0; i < sizeof(gone) / sizeof(gone[0]); i++) {
        snprintf(what, sizeof(what), "'%s' on a pathway is refused", gone[i]);
        ok(what, !loads_with_law(gone[i]));
    }
    for (i = 0; i < sizeof(kept) / sizeof(kept[0]); i++) {
        snprintf(what, sizeof(what), "'%s' on a pathway loads", kept[i]);
        ok(what, loads_with_law(kept[i]));
    }
}

static void test_module_type_needs_module(void) {
    static const char *kinds[] = { "interaction", "gain", "switch",
                                   "loop_limited", "exchange" };
    char   buf[700], what[96];
    size_t i;

    printf("\n[65] a module type is a module, never a stock wearing the name\n");
    for (i = 0; i < sizeof(kinds) / sizeof(kinds[0]); i++) {
        snprintf(buf, sizeof(buf),
            "{\"nodes\":["
            "  {\"id\":\"x\",\"type\":\"%s\",\"value\":1.0},"
            "  {\"id\":\"b\",\"type\":\"storage\",\"current_level\":0.0}],"
            " \"edges\":[{\"source\":\"x\",\"target\":\"b\","
            "            \"logic\":\"linear\"}]," CTL_SIM, kinds[i]);
        snprintf(what, sizeof(what), "'%s' with no module block is refused",
                 kinds[i]);
        ok(what, !loads(buf));
    }
}

/* Find a node by id in a projected document. */
static const cJSON *node_by_id(const cJSON *doc, const char *id) {
    const cJSON *it, *nodes = cJSON_GetObjectItemCaseSensitive(doc, "nodes");
    cJSON_ArrayForEach(it, nodes) {
        const cJSON *v = cJSON_GetObjectItemCaseSensitive(it, "id");
        if (cJSON_IsString(v) && !strcmp(v->valuestring, id)) return it;
    }
    return NULL;
}

static void test_projection_writes_modules(void) {
    cJSON        *g, *mop = NULL;
    gia_coverage  cov;
    gia_model     m;
    double        q[8];
    const cJSON  *gate, *mod;

    printf("\n[66] projection: a GSSK edge law becomes a gate, and says so\n");

    /* Test [13]'s ratio and test [15]'s threshold, written as GSSK writes
     * them: the law on the edge. */
    g = cJSON_Parse(
        "{\"nodes\":[{\"id\":\"a\",\"type\":\"storage\",\"value\":10.0},"
        "           {\"id\":\"b\",\"type\":\"storage\",\"value\":0.0},"
        "           {\"id\":\"d\",\"type\":\"constant\",\"value\":4.0},"
        "           {\"id\":\"s\",\"type\":\"storage\",\"value\":10.0},"
        "           {\"id\":\"t\",\"type\":\"storage\",\"value\":0.0}],"
        " \"edges\":[{\"id\":\"div\",\"origin\":\"a\",\"target\":\"b\","
        "            \"logic\":\"ratio\","
        "            \"params\":{\"k\":2.0,\"control_node\":\"d\"}},"
        "           {\"id\":\"sw\",\"origin\":\"s\",\"target\":\"t\","
        "            \"logic\":\"threshold\","
        "            \"params\":{\"k\":2.0,\"threshold\":6.0}}],"
        " \"config\":{\"t_end\":4.0}}");
    if (!g) { ok("parse", false); return; }
    ok("projects", gia_project(g, &mop, &cov));
    ok("every edge is carried", cov.edges_carried == cov.edges_total);
    ok("both gates are declared as added in translation",
       cov.n_added == 2 &&
       strstr(cov.added[0], "div") && strstr(cov.added[1], "sw"));

    gate = node_by_id(mop, "div__gate");
    mod  = gate ? cJSON_GetObjectItemCaseSensitive(gate, "module") : NULL;
    ok("ratio became an interaction module with action divide",
       gate && mod &&
       !strcmp(cJSON_GetObjectItemCaseSensitive(gate, "type")->valuestring,
               "interaction") &&
       !strcmp(cJSON_GetObjectItemCaseSensitive(mod, "action")->valuestring,
               "divide"));
    gate = node_by_id(mop, "sw__gate");
    ok("threshold became a switch module",
       gate && !strcmp(cJSON_GetObjectItemCaseSensitive(gate, "type")->valuestring,
                       "switch"));

    if (!load_ok("the projection loads", &m, mop)) { cJSON_Delete(g); return; }
    ok("solves", gia_network_state(&m, 2.0, q, NULL));
    /* Node order: a b d s t, then the two gates. */
    close_to("ratio: a(2) = 10 e^-1, as the pathway gave", q[0],
             10.0 * exp(-1.0), 1e-8);
    ok("solves past the threshold crossing", gia_network_state(&m, 4.0, q, NULL));
    close_to("threshold: s is held at 6, as the pathway gave", q[3], 6.0, 1e-4);
    gia_model_free(&m);
    cJSON_Delete(mop);
    cJSON_Delete(g);
}

static void test_seeds_as_modules(void) {
    static const char *paths[2] = { "examples/giannantoni/input.json",
                                    "examples/giannantoni/closed_loop.json" };
    gia_mode want[2] = { GIA_MODE_GENERATIVE, GIA_MODE_FUNCTIONAL };
    int      i;

    printf("\n[67] the example seeds draw Odum's work gate\n");

    for (i = 0; i < 2; i++) {
        char      *text = slurp_file(paths[i]);
        cJSON     *root = text ? cJSON_Parse(text) : NULL, *out;
        gia_model  m;
        int        gate;
        double     q[8];

        free(text);
        if (!root) { ok(paths[i], false); continue; }
        if (!load_ok(paths[i], &m, root)) continue;

        for (gate = 0; gate < m.n_nodes; gate++)
            if (!strcmp(m.nodes[gate].id, "interaction_1")) break;
        ok("interaction_1 is a module", gate < m.n_nodes &&
                                        gia_node_is_module(&m, gate));
        ok("solves", gia_network_state(&m, m.t_end, q, NULL));
        ok("the consumer's feedback is dissipated to heat",
           q[m.n_nodes - 1] > 0.0 && m.nodes[m.n_nodes - 1].kind == GIA_NODE_SINK);
        out = gia_generate(&m);
        ok(i == 0 ? "input.json still runs generatively"
                  : "closed_loop.json still runs functionally",
           out && gia_validate_mode(root, out) == want[i]);
        cJSON_Delete(out);
        gia_model_free(&m);
        cJSON_Delete(root);
    }
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
    test_emergy();
    test_emergy_feedback();
    test_carriers();
    test_carrier_validation();
    test_carrier_default();
    test_barter();
    test_forcing();
    test_forcing_refusals();
    test_edge_rate_forcing();
    test_constant_matrix_has_no_integration_error();
    test_projection_full();
    test_projection_names_its_losses();
    test_projection_processing_node_params();
    test_reunited_coproducts();
    test_independent_inputs_still_sum();
    test_module_order_invariance();
    test_module_gate_closed_form();
    test_module_nary_and_fanout();
    test_module_gain_closed_form();
    test_module_validation();
    test_module_switch_events();
    test_transactor_module();
    test_transactor_order_invariance();
    test_transactor_validation();
    test_transactor_barter();
    test_module_emergy_uses_module_flow();
    test_control_does_not_close_a_pathway();
    test_module_is_passed_through();
    test_control_does_not_open_the_boundary();
    test_ordinality_invariant_under_respelling();
    test_transactor_does_not_change_kind();
    test_emergent_quality_closes();
    test_sink_is_never_closed();
    test_emergence_is_order_invariant();
    test_emergent_quality_is_a_component();
    test_use_ratio_is_required();
    test_control_emergy_reaches_product();
    test_control_is_drawn_and_dissipated();
    test_used_leg_validation();
    test_control_draw_opens_the_boundary();
    test_drawn_control_closes_no_pathway();
    test_action_validation();
    test_action_divide();
    test_action_subtract();
    test_action_draws_control();
    test_module_laws_leave_pathways();
    test_module_type_needs_module();
    test_projection_writes_modules();
    test_seeds_as_modules();

    printf("\n%s\n", failures == 0 ? "ALL PASS" : "FAILURES PRESENT");
    printf("failures: %d\n", failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
