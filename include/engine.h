/* engine.h — Giannantoni Generative Computational Framework.
 *
 * A self-contained implementation of the two constructions in
 *
 *   Giannantoni, C. (2006) "Mathematics for Generative Processes:
 *       Living and Non-Living Systems", J. Comp. Appl. Math. 189.
 *   Giannantoni, C. (2023) "Generativity of Self-Organizing Processes
 *       and Their Maximum Ordinality" (preprint).
 *
 * This translation unit deliberately depends on nothing but cJSON and the
 * C99 standard library. It does NOT include "gssk.h".
 *
 * The reason is not tidiness. The GSSK kernel carries system state as a flat
 * `double *` of scalar storages; Giannantoni's Relational Space carries
 * ordinal coordinates {r} = e^(sigma i (+) phi j (+) theta k), which are
 * complex and, at higher order, non-commutative. Threading the latter through
 * an API shaped for the former is the integration cost we are choosing not to
 * pay. The two engines share a build and a JSON parser, nothing else.
 *
 * See docs/giannantoni_assessment.md for how this relates to the kernel's
 * Phase 1 IDC solver, which solves a different (scalar, conservative) problem.
 */

#ifndef GIANNANTONI_ENGINE_H
#define GIANNANTONI_ENGINE_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

#include "cJSON.h"

/* ------------------------------------------------------------------ *
 * Limits
 * ------------------------------------------------------------------ */

/* phi(t) is carried as a polynomial, so every derivative it needs is exact.
 * Degree 8 is far past the order at which the drift term stops being the
 * interesting part of the comparison. */
#define GIA_PHI_MAX_DEGREE   8

/* Highest incipient/traditional derivative order the Bell recursion runs to. */
#define GIA_MAX_ORDER        8

/* Branch cap for the n-et. The binary (duet) case is n = 2. */
#define GIA_MAX_BRANCHES     16

/* ================================================================== *
 * 1. Exponential form:  f(t) = e^phi(t)
 *
 * Giannantoni's persistence-of-form theorem is stated for the exponential
 * form, and any positive f can be put in it via phi = ln f. Carrying phi as
 * a polynomial keeps every phi^(n) exact, which is what makes the drift below
 * a computation rather than an assertion.
 * ================================================================== */

typedef struct {
    double c[GIA_PHI_MAX_DEGREE + 1];  /* phi(t) = sum_k c[k] t^k */
    int    degree;
} gia_phi;

/* phi(t) itself. */
double gia_phi_eval(const gia_phi *p, double t);

/* d^order phi / dt^order at t, computed analytically from the coefficients.
 * order 0 returns phi(t). Orders above the degree return exactly 0. */
double gia_phi_deriv(const gia_phi *p, int order, double t);

/* ================================================================== *
 * 2. Incipient vs. traditional derivative, and the drift between them
 *
 * Incipient (2006, persistence of form):
 *
 *     (d~/d~t)^n e^phi(t)  =  [phi^o(t)]^n e^phi(t)
 *
 * Traditional (Faa di Bruno, via the complete Bell polynomial B_n):
 *
 *     (d/dt)^n  e^phi(t)  =  B_n(phi', phi'', ..., phi^(n)) e^phi(t)
 *
 * The amplitude of the incipient derivative is a bare power of phi'. The
 * traditional one accumulates every partition of n into cross-derivative
 * terms. Their difference is the derivative drift:
 *
 *     psi_n(t)  =  B_n(phi', ..., phi^(n))  -  (phi')^n
 *
 * This is the whole of Giannantoni's objection to TDC, and it is decidable:
 * psi_n vanishes identically when phi is affine in t (constant-coefficient
 * case, where the two calculi agree) and is non-zero the moment phi carries
 * a genuine t^2 or higher term (variable-coefficient case, where TDC drifts
 * away from the generative solution). gia_drift() reports it directly.
 * ================================================================== */

/* [phi^o]^n — the incipient amplitude factor. */
double gia_idc_amplitude(const gia_phi *p, int n, double t);

/* B_n(phi', ..., phi^(n)) — the traditional amplitude factor. */
double gia_tdc_amplitude(const gia_phi *p, int n, double t);

/* Full derivatives, i.e. amplitude * e^phi(t). */
double gia_idc_derivative(const gia_phi *p, int n, double t);
double gia_tdc_derivative(const gia_phi *p, int n, double t);

/* psi_n(t) = traditional amplitude - incipient amplitude. Exactly 0 for an
 * affine phi at every order; non-zero otherwise. */
double gia_drift(const gia_phi *p, int n, double t);

/* True when phi is affine in t, i.e. when the two calculi provably coincide
 * at every order. Cheap predicate over the coefficients, no evaluation. */
bool gia_drift_free(const gia_phi *p);

/* ================================================================== *
 * 3. Binary / duet / n-et functions
 *
 * A fractional incipient derivative of order p/q does not produce a single
 * function. Because
 *
 *     (d~/d~t)^(p/q) e^phi  =  (phi^o)^(p/q) e^phi
 *
 * and (phi^o)^(p/q) has q branches, the result is a q-tuple carried as one
 * object. For p/q = 1/2 and phi = alpha t this is the +/- sqrt(alpha) e^(alpha t)
 * duet of the 2006 paper. The branches are a single state, not q alternatives
 * to choose between.
 * ================================================================== */

typedef struct {
    double complex branch[GIA_MAX_BRANCHES];
    int            count;      /* = q, the denominator of the fractional order */
} gia_net;

/* (d~/d~t)^(num/den) e^phi(t), evaluated at t. `den` is clamped to
 * [1, GIA_MAX_BRANCHES]. den == 1 yields the ordinary single-branch case;
 * den == 2 yields the binary/duet. */
gia_net gia_incipient_fractional(const gia_phi *p, int num, int den, double t);

/* A duet is exactly two branches. */
bool gia_net_is_binary(const gia_net *n);

/* Sum over branches. For a complete set of q-th roots this cancels to ~0,
 * which is the balance condition of the n-et. */
double complex gia_net_sum(const gia_net *n);

/* ================================================================== *
 * 4. MOP: the Harmony Relationships
 *
 * At Maximum Ordinality the N x N matrix of ordinal relationships is not
 * free. Every off-diagonal entry is generated from a single reference
 * interaction couple alpha_12 by the (N-1) Ordinal Roots of Unity,
 *
 *     omega_m = exp(2 pi i m / (N-1)),   m = 0 .. N-2
 *
 * so that in row i, the N-1 partners j != i take alpha_ref * omega^m in
 * ascending order of j. Two consequences are worth stating because both are
 * checkable rather than decorative:
 *
 *   - Every row sums to alpha_ref * sum_m omega^m = 0 exactly, for N-1 >= 2.
 *     That is the global balance the Harmony Relationships assert.
 *   - Any entry is recoverable from alpha_ref and (i, j) alone. That is the
 *     O(N^2) -> O(1) reduction, and gia_harmony_reconstruct() lets a test
 *     confirm the stored matrix really does carry no extra information.
 *
 * Note alpha_{0,1} == alpha_ref by construction: row 0's first partner is 1,
 * taking m = 0. The reference couple is a genuine entry of the matrix.
 * ================================================================== */

typedef struct {
    int             n;          /* number of components N */
    double complex  alpha_ref;  /* the reference couple alpha_12 */
    double complex *a;          /* N*N, row-major; diagonal is 0 */
} gia_harmony;

/* The m-th of the `roots` roots of unity. */
double complex gia_ordinal_root(int roots, int m);

/* Build the N x N harmony matrix from the reference couple. Requires n >= 2.
 * Returns false on bad argument or allocation failure. */
bool gia_harmony_init(gia_harmony *h, int n, double complex alpha_ref);
void gia_harmony_free(gia_harmony *h);

/* Stored entry (i, j). Out-of-range indices return 0. */
double complex gia_harmony_at(const gia_harmony *h, int i, int j);

/* Entry (i, j) recomputed from alpha_ref and the roots of unity alone,
 * touching no stored state. Equality with gia_harmony_at() over all (i, j)
 * is the N x N -> 1 reduction claim. */
double complex gia_harmony_reconstruct(int n, double complex alpha_ref,
                                       int i, int j);

/* max_i |sum_j alpha_ij| — the balance residual. ~0 at Maximum Ordinality. */
double gia_harmony_row_residual(const gia_harmony *h);

/* max_ij |stored - reconstructed| — the reduction residual. ~0 by construction. */
double gia_harmony_reduction_residual(const gia_harmony *h);

/* ================================================================== *
 * 4b. The network: the matrix exponential form
 *
 * For a single component the exponential form is a scalar, f = e^phi. For a
 * network it is a matrix. Giannantoni 2023 gives Relational Space as
 *
 *     {r}_s = e^{alpha(t)},   alpha an N x N matrix of ordinal coordinates
 *
 * and the operational consequence is that a linear constant-coefficient network
 * is solved by ONE matrix exponentiation rather than by stepping:
 *
 *     Q(t) = exp(A t) Q(0)
 *
 * where A is the flow matrix assembled from Odum's edge laws. A component's
 * trajectory therefore depends on its neighbours', which is what makes this a
 * simulator rather than N independent scalar functions.
 *
 * The drift follows from the same object rather than from a stipulated
 * exponent. For constant A, d/dt exp(At) = A exp(At), so the traditional and
 * incipient derivatives agree identically and psi is exactly zero -- the
 * correct answer for a constant-coefficient network, and the matrix analogue of
 * "affine phi is drift-free". psi becomes non-zero exactly when A varies, with
 * t or with Q; a work gate is the latter case.
 * ================================================================== */

/* Odum's pathway laws. The names follow the kernel's GSSK_LogicType so the two
 * vocabularies do not drift apart. */
typedef enum {
    GIA_LOGIC_LINEAR,      /* Odum 1972 SecIII barbed:   F = k Q_origin        */
    GIA_LOGIC_INTERACTION, /* Odum 1972 SecX work gate:  F = k Q_origin Q_ctl  */
    GIA_LOGIC_REVERSIBLE,  /* Odum 1972 SecIII barb-less: F = k (Q_a - Q_b)    */
    GIA_LOGIC_CONSTANT,    /* F = k                                            */
    GIA_LOGIC_LIMIT,       /* Odum 1972 SecXIII receptor: F = k Q C / (C + Q)  */
    GIA_LOGIC_GAIN,        /* Odum 1972 SecIX amplifier: F = k Q_control.
                            * The control signal sets the rate; the origin
                            * supplies the power. Unlike `linear`, the flow does
                            * NOT scale with the origin's own quantity. */
    GIA_LOGIC_RATIO,       /* ADR 0002 divisor action (Odum 2000 Fig 2.6d):
                            * F = k Q_origin / max(Q_control, eps).           */
    GIA_LOGIC_SUBTRACT,    /* ADR 0008 subtracting action (Fig 2.6e):
                            * F = max(0, k (Q_origin - Q_control)). Barbed, so
                            * clamped at zero -- a negative flow would drain the
                            * target along a line whose barb says it cannot.   */
    GIA_LOGIC_THRESHOLD,   /* Odum 1972 SecXI switch: F = k when Q_origin
                            * exceeds `threshold`, else 0. Discontinuous, so it
                            * is solved piecewise between located crossings.   */
    GIA_LOGIC_EXCHANGE     /* Odum 1972 SecXV transactor: J_energy = P J_currency.
                            * A primary flow origin->target, plus a counter-flow
                            * of currency running the OTHER way between the two
                            * currency legs, of magnitude F / price.           */
} gia_logic;

const char *gia_logic_name(gia_logic l);

/* ---- forcing ----
 *
 * ADR 0006's vocabulary: one set of waveforms, attachable in two places. A
 * forcing on a NODE drives that component's held value (Odum's X or N, a
 * force); a forcing on an EDGE drives that pathway's rate (Odum's J, a flow).
 * This engine implements the node attachment; see the note on
 * gia_flow_matrix_is_constant for why the edge attachment is a different
 * mathematical problem and is not here.
 *
 * These three waveforms are chosen because each is its OWN GENERATOR, so the
 * forcing can be carried as extra state and the system stays linear and
 * time-invariant:
 *
 *     sine         d/dt [s; c] = [[0, w], [-w, 0]] [s; c]
 *     ramp         d/dt r = 1
 *     exponential  d/dt e = lambda e
 *
 * That is what keeps Q(t) = exp(A t) Q(0) exact under forcing: the driver is
 * absorbed into A rather than making A depend on t. A waveform that is not its
 * own generator -- a square wave, jitter -- is refused rather than approximated
 * silently. */
typedef enum {
    GIA_FORCE_NONE,        /* the declared value is used as-is (default) */
    GIA_FORCE_SINE,        /* offset + amplitude * sin(w t + phase)      */
    GIA_FORCE_RAMP,        /* offset + rate * t                          */
    GIA_FORCE_EXPONENTIAL  /* offset + amplitude * exp(rate * t)         */
} gia_forcing_kind;

typedef struct {
    gia_forcing_kind kind;
    double amplitude;
    double rate;       /* angular frequency for sine; growth rate otherwise */
    double phase;
    double offset;
} gia_forcing;

/* The driver's value at t, for reporting. The solver does not call this: it
 * carries the waveform as state so the solution stays closed-form. */
double gia_forcing_value(const gia_forcing *f, double base, double t);

/* How a component divides emergy among its outgoing pathways. Declared here
 * because gia_edge carries it; the accounting that uses it is in section 8b. */
typedef enum {
    GIA_OUT_PARTITION,   /* a split: proportional. Emergy is conserved.   */
    GIA_OUT_REPLICATE    /* a co-product: the whole, to each. It is not.  */
} gia_output_mode;

/* Dense N x N flow matrix, row-major, owned by the caller via gia_matrix_free. */
typedef struct {
    int     n;
    double *a;
} gia_matrix;

bool gia_matrix_init(gia_matrix *m, int n);
void gia_matrix_free(gia_matrix *m);
double gia_matrix_at(const gia_matrix *m, int i, int j);

/* B = exp(M), by Pade (3,3) with scaling and squaring.
 * Both matrices must be n x n and distinct. Returns false on bad input or
 * allocation failure. */
bool gia_matrix_exp(const gia_matrix *M, gia_matrix *B);

/* ================================================================== *
 * 5. The model graph
 *
 * Odum energy-systems graph read from the JSON seed. String fields are
 * borrowed from the cJSON tree, so the tree must outlive the model.
 * ================================================================== */

/* Odum's fundamentals, named as the kernel's schema names them so the two
 * vocabularies cannot drift. A type whose law this engine does not implement is
 * rejected at load rather than accepted as a label -- see gia_model_load. */
typedef enum {
    GIA_NODE_SOURCE,       /* Odum 1972 SecII   — held at its value, not solved */
    GIA_NODE_STORAGE,      /* Odum 1972 SecVI   — tank, dQ/dt = sum J          */
    GIA_NODE_SINK,         /* Odum 1972 SecV    — accumulates, never depleted  */
    GIA_NODE_CONSTANT,     /* no Odum symbol    — held, read, never consumed   */
    GIA_NODE_INTERACTION,  /* Odum 1972 SecX    — work gate                    */
    GIA_NODE_GAIN,         /* Odum 1972 SecIX   — constant gain amplifier      */
    GIA_NODE_LOOP_LIMITED, /* Odum 1972 SecXIII — cycling receptor             */
    GIA_NODE_SWITCH,       /* Odum 1972 SecXI   — on/off threshold process     */
    GIA_NODE_EXCHANGE,     /* Odum 1972 SecXV   — transaction diamond          */
    GIA_NODE_CONSUMER,     /* GSSK composite, not expanded here                */
    GIA_NODE_UNKNOWN
} gia_node_kind;

typedef struct {
    const char   *id;
    const char   *label;
    gia_node_kind kind;
    gia_phi       phi;       /* single-component analytic form (Sections 1-3) */
    double        q0;        /* initial quantity, for the network solution     */
    double        quality_input; /* Tr injected by a source; 0 if none        */
    gia_forcing   forcing;   /* drives the held value; ADR 0006 node attachment */
    const char   *carrier;   /* what this component holds; "" is the implicit
                              * single carrier, so a model that names none
                              * behaves exactly as before. Borrowed from the
                              * cJSON tree.                                    */
    bool          integrates;/* false for source/constant: Q is held, not solved */
    bool          on_cycle;  /* filled in by gia_mark_cycles() */
} gia_node;

typedef struct {
    int         from;        /* index into nodes, or -1 if unresolved */
    int         to;
    const char *flow_type;   /* the seed's label; maps onto `logic` */
    gia_logic   logic;
    double      weight;      /* the rate coefficient k */
    int         control;     /* control node index, or -1 */
    double      capacity;    /* C, for limit edges */
    double      threshold;   /* crossing level, for threshold edges */
    double      price;       /* P, for exchange edges */
    int         cur_from;    /* currency leg paying, for exchange edges */
    int         cur_to;      /* currency leg receiving, for exchange edges */
    gia_output_mode out_mode;/* partition (split) or replicate (co-product) */
    gia_forcing  forcing;    /* drives the RATE k; ADR 0006 edge attachment.
                              * Unlike the node attachment this is NOT
                              * absorbable: the flow becomes k(t)*Q, bilinear
                              * in driver and state, so the matrix genuinely
                              * depends on t and the closed form ends here. */
} gia_edge;

typedef struct {
    const char *system_name;
    gia_node   *nodes;
    int         n_nodes;
    gia_edge   *edges;
    int         n_edges;
    double      t_end;        /* simulation_params.t_val */
    int         order;        /* simulation_params.derivative_order */
    bool        generative;   /* simulation_params.generative_mode */
    cJSON      *root;         /* borrowed, not owned */
} gia_model;

const char *gia_node_kind_name(gia_node_kind k);

/* Populate `m` from a parsed seed graph. `root` is borrowed and must outlive
 * `m`. Returns false if the document is not a usable graph. */
bool gia_model_load(gia_model *m, cJSON *root);

/* Release only what gia_model_load allocated. Does not touch `root`. */
void gia_model_free(gia_model *m);

/* Assemble the flow matrix at the operating point `q` (length n_nodes).
 *
 * The matrix is AUGMENTED: it is (n+1) x (n+1), and the extra row and column
 * carry the constant part. Some of Odum's pathways deliver a rate that does not
 * scale with any quantity -- a `constant` pathway, and a `threshold` pathway
 * while it is open -- so the system is affine, dQ/dt = A Q + b, not linear. The
 * standard device is to solve
 *
 *     d/dt [Q; 1] = [[A, b], [0, 0]] [Q; 1]
 *
 * which is exactly linear again in n+1 dimensions, so one matrix exponential
 * still gives the closed-form solution. Column n holds b; row n is zero, which
 * pins the phantom component at 1.
 *
 * Returns false on allocation failure. */
bool gia_build_flow_matrix(const gia_model *m, const double *q, double t,
                           gia_matrix *out);

/* True when the flow matrix is the same at two different operating points, i.e.
 * when A is constant and the incipient solution is exact and closed-form. */
bool gia_flow_matrix_is_constant(const gia_model *m);

/* Q(t) = exp(A t) Q(0), written into `out` (length n_nodes). For a constant A
 * this is exact. For a state-dependent A it is the Pade linearisation about
 * Q(0), and `out_drift`, when non-NULL, receives the largest component-wise
 * difference between the incipient and traditional derivatives at t -- exactly
 * zero when A is constant.
 *
 * Discontinuous pathways (`threshold`, and the clamp on `subtract`) are solved
 * PIECEWISE. The persistence-of-form theorem needs a smooth alpha, and there is
 * none across a crossing, so the run is broken at each located crossing and an
 * exact exponential is taken over each smooth interval. Crossings are found by
 * the Illinois method, as the kernel does. */
bool gia_network_state(const gia_model *m, double t, double *out,
                       double *out_drift);

/* Richardson estimate of the integration error in Q(t), from the difference
 * between the composed solution and one computed with twice as many
 * subintervals, scaled by 4/3 because the step is second order and the raw
 * difference is only three quarters of the error it is estimating.
 *
 * Exactly zero when the flow matrix is constant, because there is nothing being
 * composed -- the single exponential IS the answer. Non-zero only where the
 * engine had to compose, and it exists so that a reported psi can be read
 * against it: a drift smaller than the integration error of the trajectory it
 * came from is not evidence of anything. */
double gia_integration_error(const gia_model *m, double t);

/* Number of switching crossings located over [0, t]. Zero for a smooth model.
 * A non-zero count is the reason the solution over that span is piecewise
 * rather than a single closed form. */
int gia_count_events(const gia_model *m, double t);

/* ---- carriers ----
 *
 * A component holds a quantity OF something. Summing a store of grain and a
 * bank balance into one total is not a conservation check, it is a category
 * error -- and it is what this engine did until carriers existed: goods and
 * money went into the same running total, so a violation in one could be
 * cancelled by a movement in the other and the residual would report zero.
 *
 * Carriers are named per component. A model that names none has a single
 * implicit carrier "" and behaves exactly as it did before.
 *
 * A pathway may not cross carriers. Grain does not turn into money by flowing
 * along an edge; it is exchanged for money, and that coupling is what Odum's
 * transaction diamond (1972 SecXV) is for. A cross-carrier edge whose law is
 * not `exchange` is therefore rejected at load rather than quietly moving
 * quantity between incommensurable stocks. */

/* Number of distinct carriers in the model; at least 1. */
int gia_carrier_count(const gia_model *m);

/* Name of carrier `idx`; "" for the implicit single carrier. */
const char *gia_carrier_name(const gia_model *m, int idx);

/* Index of the carrier a component holds. */
int gia_node_carrier(const gia_model *m, int node_idx);

/* Conservation residual for ONE carrier over [0, t]. */
double gia_conservation_residual_for(const gia_model *m, double t, int carrier);

/* Worst residual across all carriers over [0, t]. Reported, never silently
 * corrected.
 *
 * This is deliberately a maximum and not a sum: a sum lets a surplus in one
 * carrier mask a deficit in another, which is precisely the defect carriers
 * exist to remove.
 *
 * It is a conservation *check* only for a closed system. Odum holds a source at
 * its value rather than integrating it (1972 SecII), so a pathway leaving a
 * source delivers quantity without being depleted: in an open system this
 * residual is the net inflow across the boundary, and a non-zero value is the
 * correct answer rather than an error. Use gia_system_is_closed() to know which
 * reading applies.
 *
 * Emergy is deliberately outside this. It is non-conservative by construction,
 * and that non-conservativeness is what Giannantoni's argument for IDC rests
 * on -- folding the two accountings together would erase the thing being
 * demonstrated. */
double gia_conservation_residual(const gia_model *m, double t);

/* True when no pathway originates at a held (non-integrating) component, so
 * nothing enters or leaves across the boundary and the total must be conserved. */
bool gia_system_is_closed(const gia_model *m);

/* ================================================================== *
 * 8b. Emergy and transformity — the second accounting
 *
 * Quantity is conserved. Emergy is not, and that is the whole point.
 *
 * Giannantoni's argument for IDC (2006) begins from the observation that Odum's
 * Emergy Algebra could not be written as a dynamic differential equation using
 * classical derivatives, because TDC imposes conservative relationships that
 * emergy violates by construction: the output of a generative process retains
 * the genetic structure of its inputs but is irreducible to them. The
 * `d~/d~t` generator is offered as the operator that encodes exactly that.
 *
 * So an engine implementing his calculus without emergy is implementing the
 * machinery without the quantity that motivated it. This is that quantity.
 *
 * The rules follow Odum. Emergy is carried along a pathway as F x Tr_origin.
 * At a bifurcation the treatment depends on what kind it is, and the two are
 * NOT the same:
 *
 *   partition (a split, the default) — one kind of flow divided between
 *       branches. Each branch takes emergy in proportion to its share of the
 *       energy. Emergy out equals emergy in, so this is the conservative case.
 *
 *   replicate (a co-production) — genuinely different products of one process.
 *       EACH takes the whole emergy, not a share, because each required all of
 *       it to exist. Emergy out exceeds emergy in. This is the irreducible
 *       excess, and it is where the accounting stops being conservative.
 *
 * Feedback is not counted twice (Odum's fourth rule): emergy propagates from
 * the sources along the acyclic part of the graph, and an edge closing a cycle
 * carries quantity without re-injecting emergy that has already been counted.
 *
 * The vocabulary is the neutral one from docs/emergy_synthesis.md --
 * `quality_input` on a source, `output_mode` of "partition" or "replicate" on
 * an edge -- so the two engines cannot drift apart on terminology.
 * ================================================================== */

/* Flow along one pathway at the operating point `q` and time `t`, by that
 * pathway's law. This is the quantity F that the emergy pass carries Tr along.
 *
 * `t` is needed because a driven component's value is not in `q`: a held
 * component never moves, so its instantaneous value comes from its waveform. */
double gia_edge_flow(const gia_model *m, const gia_edge *e, const double *q,
                     double t);

/* Empower (emergy per unit time) and transformity per component at time t.
 * Either output array may be NULL; both are length n_nodes.
 *
 * `tr[i]` is emergy per unit quantity, so it is the *quality* of component i,
 * and it is why two components holding the same amount are not equivalent. */
bool gia_emergy_at(const gia_model *m, double t, double *em, double *tr);

/* Emergy created at co-productions over the whole network at time t:
 *
 *     sum over components of max(0, emergy out - emergy in)
 *
 * Exactly zero when every bifurcation is a partition, and positive as soon as
 * one is a replication. This number is the non-conservativeness itself -- the
 * thing Giannantoni says a conservative calculus cannot express -- reported
 * rather than argued about. */
double gia_emergy_excess(const gia_model *m, double t);

/* ================================================================== *
 * 6. Ordinality
 *
 * A component participates in the system's ordinal structure to the extent
 * that it sits on a closed pathway — an input that never returns is a flow
 * through the system, not a relationship within it. Maximum Ordinality is
 * reached when every component lies on some feedback loop, and that is the
 * condition the generative step below works to satisfy.
 * ================================================================== */

/* Set node.on_cycle for every node; returns how many are on a cycle. */
int gia_mark_cycles(gia_model *m);

/* Fraction of components on a cycle, in [0, 1]. Calls gia_mark_cycles(). */
double gia_ordinality(gia_model *m);

/* True when every component is on a cycle. */
bool gia_at_maximum_ordinality(gia_model *m);

/* ================================================================== *
 * 7. Mode 1 — functional: trajectories over a fixed topology
 *
 * Writes Time plus, per node, the incipient value, the traditional value and
 * the drift between them, at the requested derivative order. Returns false if
 * the file cannot be written.
 * ================================================================== */

/* One sampled row of the run, filling caller-provided arrays of length
 * n_nodes. Both the CSV writer and the terminal printer go through this, so
 * the two cannot disagree about what a row contains -- they did once, when the
 * CSV gained a network column and the printer kept showing only the
 * single-component form, and a reader of the terminal saw a held source
 * apparently growing.
 *
 * `q` receives the network quantity, `idc`/`tdc` the single-component analytic
 * form. Any of the out-parameters may be NULL. */
bool gia_sample_at(const gia_model *m, double t,
                   double *q, double *idc, double *tdc, double *psi);

bool gia_write_trajectories(const gia_model *m, const char *path, int steps);

/* The same values gia_write_trajectories() writes, rendered as aligned tables
 * on stdout: one section per component, then a summary carrying each
 * component's final values and its largest drift over the run. */
void gia_print_trajectories(const gia_model *m, int steps);

/* ================================================================== *
 * 8. Mode 2 — generative: an ordinal step that changes the graph
 *
 * Returns a NEW graph the caller owns and must cJSON_Delete(). When the seed
 * is already at Maximum Ordinality, or generative mode is off, the returned
 * graph compares equal to the input and the run is functional. Otherwise a
 * regulator is spawned to close an open pathway, raising ordinality.
 * ================================================================== */

cJSON *gia_generate(const gia_model *m);

/* ================================================================== *
 * 8c. Projection — a GSSK model read into MOP terms, declaring its losses
 *
 * ADR 0011 decisions 2 and 3. The direction is one-way and forced: GSSK is the
 * Odum conformance authority, and MOP has no typed relations to invent, so a
 * MOP-to-GSSK path would have to make up the typing it lacks. There is no such
 * path here and there should not be one.
 *
 * The projection exists to produce a NUMBER. docs/giannantoni_assessment.md
 * 5.3 makes a Level 1 claim -- that any system modelled through the MOP lens
 * has an explicit solution -- and flags it as unquantified. Coverage measures
 * it against real models: what fraction of a model's components and pathways
 * the MOP engine can actually carry, and, where it cannot, WHICH module
 * stopped it.
 *
 * A bare percentage would be the wrong deliverable. A model scoring 60% does
 * not mean MOP is 60% correct; it means that model uses modules this engine
 * cannot carry, and the report has to name them. Nothing is ever silently
 * substituted -- a component that cannot be projected is dropped and recorded,
 * never approximated into something that looks like it worked.
 * ================================================================== */

#define GIA_MAX_FINDINGS 64

typedef struct {
    int  nodes_total,  nodes_carried;
    int  edges_total,  edges_carried;
    int  n_findings;
    int  n_elided;                    /* findings past the cap */
    char finding[GIA_MAX_FINDINGS][192];
} gia_coverage;

/* Read a GSSK-schema model and produce a MOP input. The caller owns *out_mop
 * and must cJSON_Delete() it. `cov` is filled in either way. Returns false only
 * on a malformed document, not on incomplete coverage: a model the projection
 * can only partly carry is a result, not an error. */
bool gia_project(const cJSON *gssk, cJSON **out_mop, gia_coverage *cov);

/* Carried elements over total, in [0, 1]. 1.0 for a model wholly within the
 * MOP engine's vocabulary. */
double gia_coverage_fraction(const gia_coverage *cov);

/* Human-readable report: the fraction, and every module that blocked it. */
void gia_report_coverage(const gia_coverage *cov);

/* ================================================================== *
 * 9. Validation — which mode did the engine actually run in?
 *
 * Decided structurally, by comparing the graph that went in against the one
 * that came out, rather than by trusting a flag.
 * ================================================================== */

typedef enum {
    GIA_MODE_FUNCTIONAL,   /* Input_Graph === Output_Graph */
    GIA_MODE_GENERATIVE    /* Input_Graph !== Output_Graph */
} gia_mode;

gia_mode gia_validate_mode(const cJSON *in, const cJSON *out);

/* Human-readable report, including the node/edge deltas. */
void gia_report_mode(gia_mode mode, const cJSON *in, const cJSON *out);

/* Check the two harmony invariants against `tol`; prints and returns pass. */
bool gia_validate_harmony(const gia_harmony *h, double tol);

#endif /* GIANNANTONI_ENGINE_H */
