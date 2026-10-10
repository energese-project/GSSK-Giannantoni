# Interface control — Giannantoni kernel

The contracts an integrator programs against and a test can check. Signatures here are normative in
shape (names, parameters, ownership, error behaviour); a PR may add parameters only by revising this
document. Existing `engine.h` functions keep their signatures; the new units use `gia_status`.

## 1. Units and headers

| Header | Unit | Contents | Includes |
|---|---|---|---|
| `include/gia_status.h` | — | `gia_status`, `gia_status_str` | `<stddef.h>` |
| `include/idc.h` | `src/idc.c` | Single-variable IDC solvers (FR-IDC-002, 006–010) | `gia_status.h`, `<complex.h>` |
| `include/relational.h` | `src/relational.c` | Relational elements (FR-REL) | `gia_status.h` |
| `include/mop.h` | `src/mop.c` | First and Second Equation, EQS, harmony detector, network β (FR-MOP, FR-HAR-001–003) | `engine.h`, `relational.h` |
| `include/engine.h` | `src/engine.c` | Model, network, emergy, ordinality, generative step (existing; FR-EM, FR-ORD additions) | unchanged |
| `src/harmony.c` | — | `gia_harmony_assume_*` (FR-HAR-004), declared in `engine.h` | — |

None includes `gssk.h` (NFR-SEP-001). `bin/test_mop_emergence` links `mop.o` and `relational.o`
without `harmony.o` (FR-HAR-003).

## 2. API contracts

### IF-API-001 — Status codes and the `why` out-parameter
```c
typedef enum {
    GIA_OK = 0,
    GIA_E_ARG,          /* NULL where not allowed, size or index out of range        */
    GIA_E_DOMAIN,       /* mathematically outside the function's stated domain       */
    GIA_E_RANGE,        /* result would overflow or be non-finite (NFR-NUM-003)      */
    GIA_E_CONVERGENCE,  /* quadrature or root tracking failed its tolerance          */
    GIA_E_UNSUPPORTED,  /* defined by no available source (srs FR-IDC-013, PLAN §6)  */
    GIA_E_LIMIT,        /* a fixed limit exceeded (NFR-LIM-001)                      */
    GIA_E_NOMEM
} gia_status;
const char *gia_status_str(gia_status s);
```
Every new function returns `gia_status` and takes `const char **why` last. On any status other than
`GIA_OK`, if `why` is non-NULL it is set to a static string naming the cause and its source, e.g.
`"[06 Eq 3.22] is not derivable from Eq 3.19 (PLAN X3)"`; outputs are left unmodified. On `GIA_OK`,
`*why` is not written. No function prints, exits or aborts (NFR-ERR-001).
- **Source:** srs NFR-ERR-001, FR-OUT-002
- **Verification:** T-API-01, T-ERR-01 · **Priority:** Must · **Status:** implemented · **Task:** idc-general-f

### IF-API-002 — `idc.h`
```c
typedef double complex (*gia_cfn)(double t, void *ctx);      /* caller-owned coefficient */

gia_status gia_idc_of(double complex f, double complex df, int n,
                      double complex *out, const char **why);                 /* FR-IDC-002 */

typedef struct gia_lde2_sol gia_lde2_sol;                                      /* opaque     */
gia_status gia_lde2_solve(gia_cfn a1, gia_cfn a0, void *ctx,
                          double complex f0, double complex f1, double t_max,
                          gia_lde2_sol **sol, double *t_fail /* may be NULL */,
                          const char **why);                                    /* FR-IDC-006 */
gia_status gia_lde2_eval(const gia_lde2_sol *sol, double t,
                         double complex *f, double complex *trad_residual,
                         const char **why);
gia_status gia_lde2_terms(const gia_lde2_sol *sol, double t, double complex c[2],
                          double complex r[2], double complex E[2], const char **why);
void       gia_lde2_free(gia_lde2_sol *sol);

typedef struct { double complex u[2]; double complex c[2][2]; } gia_binary_sol;
gia_status gia_binary_solve(double complex A, double complex B,
                            const double complex f0[2], const double complex fhalf0[2],
                            gia_binary_sol *sol, const char **why);              /* FR-IDC-007 */
gia_status gia_binary_eval(const gia_binary_sol *sol, double t,
                           double complex f[2], double complex fhalf[2],
                           const char **why);

gia_status gia_riccati_solve(gia_cfn Q, gia_cfn R, gia_cfn dR, gia_cfn P, void *ctx,
                             double complex f0, double t_max, gia_lde2_sol **sol,
                             double *t_fail, const char **why);                  /* FR-IDC-008 */
gia_status gia_riccati_eval(const gia_lde2_sol *sol, double t,
                            double complex *f, double complex *trad_residual,
                            const char **why);

gia_status gia_nl1410_roots(double complex A, double complex B,
                            double complex u[2], const char **why);              /* FR-IDC-009 */
gia_status gia_idc_taylor(double f0, double df0, double dt, int n,
                          double *out, const char **why);                        /* FR-IDC-010 */
gia_status gia_idc_refuse(const char *feature, const char **why);                /* FR-IDC-013 */
```
`t_fail` is written only on `GIA_E_CONVERGENCE`, with the time at which the roots collide: it is
the diagnosis that a static `why` cannot carry, and the one exception to "outputs are left
unmodified". `gia_lde2_terms` returns the solution's terms (constants, roots at t, and
`E = e^{∫₀ᵗ r}`) so that a test can form the termwise residual itself (vv-plan.md §2 rule 2).
Ownership: `gia_lde2_sol` is allocated by `*_solve` and freed by `gia_lde2_free`; it holds the
caller's function pointers and `ctx` without copying, so `ctx` must outlive it. `t` must lie in
`[0, t_max]`, else `GIA_E_ARG`.
- **Source:** srs §2.1
- **Verification:** T-API-02 · **Priority:** Must · **Status:** planned · **Task:** idc-lde2

### IF-API-003 — `relational.h`
```c
typedef struct { double i, j, k; } rel_t;          /* x_i·i ⊕ x_j·j ⊕ x_k·k; i is the real unit */

rel_t      rel_mul(rel_t x, rel_t y);                                     /* FR-REL-001 */
rel_t      rel_mul3(rel_t x, rel_t y, rel_t z);    /* (x∘y)∘z, never x∘(y∘z)   FR-REL-004 */
gia_status rel_exp(rel_t x, rel_t *out, const char **why);                /* FR-REL-002 */
gia_status rel_root(int N, int l, rel_t *out, const char **why);          /* FR-REL-003 */
gia_status rel_root_pow(int N, int l, int m, rel_t *out, const char **why);  /* by angle */
rel_t      rel_mul_pow(rel_t x, int m);            /* by repeated ∘, left to right */
```
Pure functions; no allocation. `rel_exp` returns a status because `eᵃ` overflows (NFR-NUM-003).
- **Source:** srs §2.4
- **Verification:** T-API-02 · **Priority:** Must · **Status:** implemented · **Task:** mop-relational-algebra

### IF-API-004 — `mop.h`
```c
typedef struct { int num, den; } gia_rational;                    /* k; den ≥ 1, reduced */
typedef enum { GIA_BETA_NONE, GIA_BETA_AFFINE, GIA_BETA_SAMPLES } gia_beta_kind;
typedef struct {
    gia_beta_kind kind;
    double complex a, b; double p;                                /* AFFINE: (a + b t)^p */
    int n; const double *t; const double complex *v;              /* SAMPLES, caller-owned */
} gia_beta;

gia_status gia_mop_couple(const gia_beta *beta, gia_rational k, double t,
                          double complex *alpha, const char **why);         /* FR-MOP-001,2,4 */

typedef struct { int N; double complex *a; unsigned char *related; } gia_matrioska; /* row-major */
gia_status gia_mop_solve(int N, const gia_beta *beta /* N*N, diagonal ignored */,
                         gia_rational k, double t, gia_matrioska *out,
                         const char **why);                                 /* FR-MOP-003 */
void       gia_matrioska_free(gia_matrioska *m);
gia_status gia_mop_couple_rel(const gia_beta beta[3], gia_rational k, double t,
                              rel_t *alpha, const char **why);              /* FR-MOP-007 */

typedef double complex (*gia_quad_fn)(double t, void *ctx);
gia_status gia_quad_gk15(gia_quad_fn g, void *ctx, double lo, double hi, double tol,
                         int max_depth, double complex *integral, double *err,
                         const char **why);                                 /* NFR-NUM-005 */

typedef struct { double complex A; double complex B[2][2]; } gia_second;
gia_status gia_mop_second(double complex alpha12_0, double c1, double c2, int N, double t,
                          gia_second *out, gia_matrioska *r, const char **why);  /* FR-MOP-005 */

typedef struct { double psi1[3], psi2, eps[3], A; int N; } gia_eqs_params;
gia_status gia_eqs(const gia_eqs_params *p, rel_t ref /* Σ0, Φ0, Θ0 at t */, int l,
                   double out[3] /* ρ, φ, θ */, const char **why);          /* FR-MOP-006 */

gia_status gia_mop_network_beta(const gia_model *m, double t, gia_beta *beta /* N*N */,
                                const char **why);                           /* FR-MOP-008 */

typedef enum { GIA_H_IMPOSED, GIA_H_TRANSPORTED, GIA_H_PRESENT, GIA_H_ABSENT } gia_verdict;
const char *gia_verdict_str(gia_verdict v);                                 /* FR-OUT-003 */
gia_status gia_harmony_residual(const gia_matrioska *alpha, double *R_H,
                                const char **why);                          /* FR-HAR-001 */
typedef gia_status (*gia_construction)(const gia_beta *in, int N, gia_matrioska *out, void *ctx);
gia_status gia_harmony_verdict(gia_construction c, void *ctx, int N,
                               gia_verdict *v, const char **why);           /* FR-HAR-002 */
gia_status gia_harmony_observed(const gia_matrioska *alpha, gia_verdict *v,
                                const char **why);                          /* FR-HAR-002 */
```
`gia_matrioska` is allocated by the producing function and freed by `gia_matrioska_free`; `related[i*N+j]`
is 0 for unrelated couples, whose `a` entry is unspecified.
- **Source:** srs §2.3, §2.6
- **Verification:** T-API-02 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### IF-API-005 — `engine.h` additions
```c
typedef struct { int k, n22, n2, nhalf, nunrelated; } gia_ordinality_rec;
gia_status gia_ordinality_record(const gia_model *m, gia_ordinality_rec *r, const char **why); /* FR-ORD-001,2 */
bool       gia_at_maximum_ordinality(gia_model *m);    /* redefined: strong connectivity, FR-ORD-003 */
double     gia_closure(gia_model *m);                  /* renamed proxy, FR-ORD-004 */
gia_status gia_solution_drift(const gia_model *m, int order, double *out,
                              const char **why);                              /* FR-IDC-011 */
gia_status gia_drift_projection(const gia_model *m, double t, double dt,
                                double *out /* per node */, const char **why); /* FR-IDC-014 */
gia_status gia_emergy_source_term(const gia_model *m, double t, int node,
                                  double *phi, const char **why);               /* FR-EM-002, 004 */
gia_status gia_emergy_check_limits(const gia_model *m, const char **why);       /* NFR-LIM-001 */
gia_status gia_emergy_carried(const gia_model *m, double t, double *carried /* per edge */,
                              const char **why);                                /* FR-MOP-008 */
typedef enum { GIA_OF_SCALAR, GIA_OF_BINARY, GIA_OF_DUET, GIA_OF_DUET_BINARY } gia_oform_kind;
typedef struct { gia_oform_kind kind; int rows, cols; double v[2][2]; } gia_oform;
gia_oform  gia_oform_scalar(double a);
gia_oform  gia_oform_binary(double em_u);                                       /* FR-EM-005 */
gia_oform  gia_oform_duet(double em_u1, double em_u2);
gia_oform  gia_oform_duet_binary(double a1, double a2);
typedef struct { int rows, cols; double pair[2][2][2]; } gia_circle;
gia_status gia_circle_product(const gia_oform *a, const gia_oform *b, gia_circle *out,
                              const char **why);                                /* FR-EM-006 */
gia_oform  gia_circle_reduce(const gia_circle *c);
typedef struct { double value, weight; } gia_balance_term;
gia_status gia_emergy_global_balance(const gia_balance_term *in, int n_in,
                                     const gia_balance_term *out, int n_out,
                                     double *residual, const char **why);       /* FR-EM-007 */
gia_status gia_emergy_balance_solve(const gia_balance_term *in, int n_in,
                                    const gia_balance_term *out, int n_out,
                                    const double *phi_w, int n_phi, double *phi,
                                    const char **why);                          /* FR-EM-007 */
```
`gia_ordinality` keeps its name for one release, returns `gia_closure`, and is marked deprecated.
- **Source:** srs §2.5, FR-IDC-014
- **Verification:** T-API-02 · **Priority:** Must · **Status:** implemented · **Task:** mop-ordinality

## 3. Seed format

### IF-JSON-001 — The `mop` block
An optional top-level `"mop"` object in a Giannantoni seed. Unknown keys anywhere inside it are load
errors; so are type mismatches (`AGENTS.md`: strict parsing).

```jsonc
"mop": {
  "k": 1,                                  // integer ≥ 1, or {"num": 1, "den": 2}   (FR-MOP-002)
  "reference": ["a", "b"],                 // optional; default: first two components by id
  "beta": "network",                       // or an array of couples:              (FR-MOP-008)
  "beta": [
    {"from": "a", "to": "b", "a": 1.0, "b": 0.25, "p": 1.0},          // affine-power; a, b may be [re, im]
    {"from": "a", "to": "c", "samples": [[0.0, 1.0, 0.0], [1.0, 2.0, 0.5]]}   // [t, re, im], t increasing
  ],
  "second_equation": {"alpha12_0": [0.3, 0.1], "c1": 1.0, "c2": 0.5},   // optional (FR-MOP-005)
  "eqs": {"psi1": [1, 1, 1], "psi2": 1, "epsilon": [0, 0, 0], "A": 1}    // optional (FR-MOP-006)
}
```
`from`/`to` must name components (not modules or boundary nodes); a couple may appear once; the time
grid is the seed's existing one.
The parser is `gia_mop_seed_load` in `include/mop_seed.h`: a load error is `GIA_E_ARG`, with the JSON
path of the offending value in a caller-owned `detail` buffer.
- **Source:** srs §2.3
- **Verification:** T-IN-01, T-IN-02 · **Priority:** Must · **Status:** implemented · **Task:** mop-first-equation

## 4. Outputs

### IF-OUT-001 — The trajectory CSV
Columns, in order: `time`; per node in seed order `<id>_Q, <id>_Em, <id>_Tr, <id>_drift_proj`;
then `conservation, emergy_excess`. Changes from baseline: `_idc`, `_tdc`, `_drift` and `psi_network`
are removed (they applied the drift identity to invented per-node φ, PLAN E4); `_drift_proj` is added
(FR-IDC-014). Every value is finite (NFR-NUM-003).
- **Source:** srs FR-IDC-011, FR-IDC-014, FR-OUT-001
- **Verification:** T-OUT-01 · **Priority:** Must · **Status:** implemented · **Task:** idc-drift-coupled

### IF-OUT-002 — The MOP CSV
Written when the seed has a `mop` block and `--mop-out PATH` is given. Columns: `time`; per related
couple in `(from, to)` id order `<from>__<to>_re, <from>__<to>_im`; then `R_H`.
- **Source:** srs §2.3, FR-HAR-001
- **Verification:** T-OUT-01 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### IF-OUT-003 — The run report
Stdout, line-oriented `key: value`. It shall include `ordinality: {k, n22, n2, nhalf, nunrelated}`,
`maximum_ordinality: yes|no`, `closure (proxy): <x>`, one `label.<column>: <label>` per output
(FR-OUT-001), and one `harmony.<construction>: <verdict>` per evaluated construction (FR-OUT-003).
- **Source:** srs §2.5–2.7
- **Verification:** T-OUT-02 · **Priority:** Must · **Status:** planned · **Task:** mop-claims-remediation

## 5. Command line

### IF-CLI-001 — `giannantoni_sim`
Existing options unchanged (`--csv`, `--out`, `--steps`, `--seed`, `--project`, `--print`, `--help`);
added `--mop-out PATH`. Exit status 0 on success, 1 on a load or validation error, 2 on a refusal
(`GIA_E_UNSUPPORTED`/`GIA_E_DOMAIN`), with the reason on stderr (FR-OUT-002).
- **Source:** `src/sim_main.c`; srs FR-OUT-002
- **Verification:** T-OUT-03 · **Priority:** Must · **Status:** implemented · **Task:** mop-first-equation
