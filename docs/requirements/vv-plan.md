# Verification and validation plan — Giannantoni kernel

How every requirement in [stakeholder.md](stakeholder.md), [srs.md](srs.md) and [icd.md](icd.md) is
shown to be met. The catalogue in §7 is the single list of verification IDs; `make check-trace`
reads it.

## 1. Levels and methods

| Level | What it exercises | Binary |
|---|---|---|
| Unit | One library function against its oracle | `bin/test_mop`, `bin/test_giannantoni` |
| Integration | Seed JSON → load → solve → emergy → outputs | `bin/test_mop` |
| System | `giannantoni_sim` end to end: exit codes, CSV, report | `tests/mop_cli.sh` |
| Validation | Published results (`VAL-…`) | `bin/test_mop` |
| Structural | Symbols, includes, static data, coverage, trace | `scripts/*.sh`, `make` targets |

| Method | Prefix | Meaning |
|---|---|---|
| Test | `T-` | Executed; pass/fail against an oracle |
| Validation | `VAL-` | A test whose oracle is a number or form printed in a source |
| Inspection | `INS-` | A mechanical check of code or docs (grep, `nm`), automated where possible |
| Analysis | `ANA-` | A derivation or probe in `docs/sources/probes/`, reviewed |
| Demonstration | `DEM-` | Observed behaviour across environments (CI matrix) |

## 2. Oracle policy

1. **No golden files.** No Giannantoni test compares against output the code produced. Oracles are a
   source equation's residual, a hand-derived closed form, a number printed in a source, or an
   independent computation in the probe scripts. `make test-update` never touches Giannantoni tests.
2. **Residual first.** A solver test always checks the defining equation's residual on the solver's
   output; a closed form is a second check, never the only one.
3. **Every test can fail.** Each catalogue entry names a mutation that must make it fail; the PR that
   adds the test quotes that failure (PLAN G7).
4. **Errata are asserted, not skipped.** For every erratum X1–X12, a test asserts that the printed form
   fails the oracle (VAL-08).
5. **Red before green.** The PR that implements a requirement quotes the test's failing run from before
   the implementation.

## 3. Tolerances

| Kind | Tolerance | Applies to |
|---|---|---|
| Closed form, well-conditioned (numerics.md) | ≤ 1e-12 relative | N1, N4, N8, N9, ordinal forms |
| Defining-equation residual (central difference, `h = 1e-5·max(1, t)`) | ≤ 1e-6 relative | T-MOP-01, T-IDC-04, T-IDC-06, T-IDC-07 |
| Quadrature | ≤ 1e-9 relative (requested 1e-10) | N2, N3 |
| Published numbers | The printed precision | VAL-01, VAL-04, VAL-05 |
| "Must fail" assertions | ≥ 100× the matching tolerance | Errata, negative controls |

A tolerance is never widened without a PLAN §5 erratum and a PLANLOG entry.

## 4. Determinism and portability

- **Same binary:** byte-identical output on repeated runs (T-DET-01).
- **Across CI toolchains** (macOS clang, Linux clang, Linux GCC): every test at its tolerance
  (DEM-POR-01). Differences from libm and from FMA contraction on arm64 clang are within §3's
  tolerances by construction (numerics.md §3).

## 5. Existing tests: disposition

| Function in `tests/test_giannantoni.c` | Disposition | Catalogue |
|---|---|---|
| `test_drift` | Keep; tag | T-IDC-01 |
| `test_duet` | Keep; tag | T-IDC-02, VAL-06 |
| `test_harmony` | Kept as a test of the constructor only; renamed `test_harmony_constructor` (FR-HAR-004 landed) | T-HAR-08 |
| `test_emergy` | Keep; tag | T-EM-01 |
| `test_emergy_feedback`, `test_reunited_coproducts`, `test_independent_inputs_still_sum` | Keep; tag | T-EM-02 |
| `test_generative`, `test_generative_disabled` | Keep; tag. When FR-ORD-005 lands, keep the mode assertions and revise the ADR 0015 ones ("one component emerged") | T-ORD-05 |
| `test_module_order_invariance`, `test_emergence_is_order_invariant` | Keep; tag | T-ORD-06 |
| `test_trajectories`, `test_printer_matches_csv` | Revise when IF-OUT-001 removes `_idc/_tdc/_drift` | T-OUT-01 |
| `test_emergent_quality_closes`, `test_sink_is_never_closed`, `test_emergent_quality_is_a_component` | Retire when FR-ORD-005 replaces ADR 0015's emergent component | T-ORD-04 |
| All others (laws, modules, carriers, forcing, projection) | Keep; they verify GSSK-vocabulary behaviour outside this baseline | — |

## 6. Harmony verdict procedure (FR-HAR-002)

For a construction c and N:

1. Build a harmonic input `β_h`, with `β_{1,j+1} = e^{2πi(j−1)/(N−1)}·β₁₂` and `β₁₂ = (1, 0.25, 1)`.
2. Build 8 perturbations `β_δ⁽ᵐ⁾ = β_h + δ·ξ⁽ᵐ⁾` with `δ = 1e-3`. ξ is drawn from a fixed 64-bit LCG
   (`x ← 6364136223846793005·x + 1442695040888963407`, seed m) mapped to `[−1, 1]²` per couple. There
   is no global RNG.
3. Evaluate c on all 9 inputs at `t = 1`.
4. Classify, with tolerance 1e-9:
   - **`imposed`** if every output has `R_H ≤ tol`;
   - **`transported`** if every output's `R_H` equals its input's within `tol·max(1, R_H(in))`, and
     some input has `R_H > tol`;
   - **`absent`** otherwise.
5. An observed Matrioska is `present` if `R_H ≤ tol`, else `absent`.

## 7. Catalogue

Columns: ID · Verifies · Level · Oracle / procedure · Mutation that must fail · Where.

### 7.1 Incipient calculus

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| T-IDC-01 | FR-IDC-001, FR-IDC-003, FR-IDC-004 | Unit | `(φ')ⁿ e^φ`, `Bₙ` by hand: `ψ₂ = a`, `ψ₃ = 3a²t` for `φ = at²/2`; affine ⇒ 0 | Swap incipient and Bell | `test_drift` (exists) |
| T-IDC-02 | FR-IDC-005 | Unit | `(d̃/dt)^{½} e^{at} = ±√a e^{at}`; q branches cancel | Drop the `2πm` | `test_duet` (exists) |
| T-IDC-03 | FR-IDC-002 | Unit | `f = 1+t²`: `(2t/(1+t²))ⁿ(1+t²)`; `f = 0` refused | Drop `/f` | `test_mop.c` |
| T-IDC-04 | FR-IDC-006, NFR-NUM-001 | Unit | (i) `a₁=0, a₀=−(1+t)²` ⇒ `α̃ = ±(1+t)`; termwise residual 0; ICs; traditional residual `= c₁e^{φ} − c₂e^{−φ}`, `φ = t+t²/2` (`a₀ = −t²`, whose roots coincide at `t = 0`, is refused; PLANLOG); also `a₀ = −(1+sin t)²` and complex `a₀ = −(1+it)²` by hand. (ii) Double root `a₁=−2t, a₀=t²`: consistent ICs solved, inconsistent refused (X12). (iii) Roots colliding at `t = ½` (`a₁ = 0, a₀ = −(t−½)²`) refused with the time | Traditional characteristic `[10 Eq 10]` | `test_mop.c` |
| T-IDC-05 | FR-IDC-006 | Unit | [10 App. Eq 28–34]: `u² + ψ_f u = 0`, `ψ_f = 1+t` ⇒ `g = C₁ + C₂e^{−∫ψ_f}`; `F = f·g̃'` not constant | Force F constant | `test_mop.c` |
| T-IDC-06 | FR-IDC-007 | Unit | Termwise residual of `f' + A f^{(½)} + B f`; four ICs; `u₁ = u₂` refused | Swap `u₁ ↔ u₂` in one half-derivative | `test_mop.c` |
| T-IDC-07 | FR-IDC-008 | Unit | Constant Q, R, P: Riccati residual of `f = y'/(Ry)` ≤ tol; printed [06 Eq 3.17] gives Eq 3.18 residual > 100× tol (X2); variable coefficients: residual reported | Printed substitution | `test_mop.c` |
| T-IDC-08 | FR-IDC-009 | Unit | `4u²+Au+B = 0`; each `e^{ut}` has zero residual; exactly two roots (X8) | Coefficient 1 for 4 | `test_mop.c` |
| T-IDC-09 | FR-IDC-010 | Unit | Horner vs direct sum; `f₀ = 0` refused; `n > 170` refused | n off by one | `test_mop.c` |
| T-IDC-10 | FR-IDC-012 | Unit | `solve(ic₁+ic₂) = solve(ic₁) + solve(ic₂)` for N3, N4 | Clamp a coefficient | `test_mop.c` |
| T-IDC-11 | FR-IDC-011 | Integration | Every constant-matrix seed: solution drift exactly 0 for n = 1–4 | `φ = ln Q_i` (nonzero on a two-mode network) | `test_mop.c` |
| T-IDC-12 | FR-IDC-014 | Integration | Constant matrix: `((A²Q)_i − (AQ)_i²/Q_i)Δ²/2` exactly; interaction seed: N7 converges and is stable as Δ → 0 | Drop the `Q'²/Q` term | `test_mop.c` |
| T-IDC-13 | FR-IDC-011, FR-IDC-013 | Unit | Refusals: [06 Eq 3.22] duet, Abel n-et, solution drift on a non-constant network; each with `why` naming the source | Return GIA_OK | `test_mop.c` |

### 7.2 Emergy algebra

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| T-EM-01 | FR-EM-001, FR-EM-002 | Integration | Split: legs 5000/5000, excess 0; co-production: each product the whole, excess one inflow | Replicate ↔ partition | `test_emergy` (exists) |
| T-EM-02 | FR-EM-003 | Integration | Reunited co-products maxed; feedback not re-injected; independent inputs sum | Sum reunions | `test_emergy_feedback`, `test_reunited_coproducts`, `test_independent_inputs_still_sum` (exist) |
| T-EM-03 | FR-EM-002 | Integration | Excess `(n−1)·Em(u)` for n = 2, 3, 4 [02 Eq 3.8] | n for n−1 | `test_mop.c` |
| T-EM-04 | FR-EM-004 | Integration | Fully drawn interaction: `Em(y) = Em(u₁)+Em(u₂)`; source term 0 | Product `k·Em₁·Em₂` | `test_mop.c` |
| T-EM-05 | FR-EM-005 | Unit | Binary, duet, duet-binary shapes and values | Transpose duet-binary | `test_mop.c` |
| T-EM-06 | FR-EM-006 | Unit | [06b Eq 2] reduced form; `l∘l` stored `[l, l]`, reduces to `l²` | Store `l²` | `test_mop.c` |
| T-EM-07 | FR-EM-007 | Unit | [02 Eq 3.23–3.26] on S = 10,000, F = 20,000, Z = 30,000, Y = 7,500: case A balances; case B gives `Φ_D = 30,000`, `Φ_E = 15,000` | Factor 4 in case B | `test_mop.c` |

### 7.3 Maximum Ordinality Principle

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| T-MOP-01 | FR-MOP-001, NFR-NUM-001 | Unit | Residual of `(α̇/α)^k α − β` on `k ∈ {1,2,3} × p ∈ {0,½,1,2} × b ∈ {0,0.25,1}`, real and complex; `k = ½`, real positive | Printed 5.5.8 | `test_mop.c` |
| T-MOP-02 | FR-MOP-001 | Unit | Printed 5.5.7/5.5.8 at (1, 0.25, 1, 2), t = 1: residual 1.25 and −0.625 (X1) | — | `test_mop.c` |
| T-MOP-03 | FR-MOP-001, NFR-NUM-002 | Unit | `α(0) = 0`; continuity `b → 0` to the `b = 0` limit within 1e-12 | Separate b = 0 branch with naive form | `test_mop.c` |
| T-MOP-04 | FR-MOP-003 | Unit | Zero diagonal; `N(N−1)` couples; unrelated ≠ zero | Fill unrelated with 0 | `test_mop.c` |
| T-MOP-05 | FR-MOP-002, NFR-NUM-004 | Unit | Each refused domain (non-integer k with complex a; `a+bt = 0`; `p = −k`) returns `GIA_E_DOMAIN` with `why` | Accept | `test_mop.c` |
| T-MOP-06 | FR-MOP-004, NFR-NUM-004 | Unit | Samples of an affine β reproduce N1 to 1e-9; β whose path crosses the negative real axis integrates continuously (no jump at the crossing); t beyond samples refused | Principal root per point | `test_mop.c` |
| T-MOP-07 | FR-MOP-005 | Unit | `Ȧ` satisfies `u' + u² = 0`; `B = [[A,−A],[−A,A]]`; `{c₂,t} → c₂t`; `c₁+c₂t ≤ 0` refused | `ln(c₁+c₂+t)` | `test_mop.c` |
| T-MOP-08 | FR-MOP-006 | Unit | Hand-evaluated N = 4 EQS; `ε₂ ≠ ε₃` refused (X11) | Anticommuting product | `test_mop.c` |
| T-MOP-09 | FR-MOP-007 | Unit | k = 1 relational componentwise; k = 2 refused | Accept k = 2 | `test_mop.c` |
| T-MOP-10 | FR-MOP-008 | Integration | Network β equals the emergy pass's empower per couple; unrelated where no pathway | Use quantity flow | `test_mop.c` |

### 7.4 Relational algebra

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| T-REL-01 | FR-REL-001 | Unit | All nine products of the table | `j∘k = −k∘j` | `test_mop.c` |
| T-REL-02 | FR-REL-004 | Unit | `rel_mul3(j, j, k) = −k`, and `j∘(j∘k) = k` (left-to-right is observable) | Reassociate | `test_mop.c` |
| T-REL-03 | FR-REL-002 | Unit | `Exp` vs cos/sin by hand; `ρ → 0` series continuity | Power series | `test_mop.c` |
| T-REL-04 | FR-REL-003 | Unit | `rel_root_pow(N, l, N−1) = 1`, N = 3…7; `rel_mul_pow(rel_root(4,1), 3) = (0.540721, 0, −0.665721)` (X10) | Power by `rel_mul` | `test_mop.c` |

### 7.5 Ordinality and generative step

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| T-ORD-01 | FR-ORD-001, FR-ORD-002 | Unit | Six hand graphs (pure split, one co-production, one interaction, both, a strongly connected graph, a module-only accumulator): counts exact; boundary nodes excluded | Count sinks | `test_mop.c` |
| T-ORD-02 | FR-ORD-003 | Unit | Two disjoint 2-cycles: closure 1, **not** maximum; strongly connected ⇒ maximum | Cycle-coverage definition | `test_mop.c` |
| T-ORD-03 | FR-ORD-004 | System | Report shows `closure (proxy)`; no decision reads it | Gate on closure | `tests/mop_cli.sh` |
| T-ORD-04 | FR-ORD-005 | Integration | Reaches maximum in ≤ #sources + #sinks steps; each choice equals brute-force argmax of total empower on all graphs ≤ 5 components from a fixed enumeration; fixed point at maximum | First candidate | `test_mop.c` |
| T-ORD-05 | FR-ORD-006 | Integration | Mode from structural diff | Trust the flag | `test_generative`, `test_generative_disabled` (exist) |
| T-ORD-06 | FR-ORD-007 | Integration | Byte-identical output under reordering | Positional legs | `test_module_order_invariance`, `test_emergence_is_order_invariant` (exist) |

### 7.6 Harmony

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| T-HAR-01 | FR-HAR-003 | Structural | `test_mop_emergence` links without `harmony.o`; `nm mop.o` lists no `gia_harmony_*`, `gia_ordinal_root` | Call the constructor from `mop.c` | Makefile, `scripts/check_symbols.sh` |
| T-HAR-02 | FR-HAR-001 | Unit | Harmonic Matrioska, N = 4, 5, 7: `R_H < 1e-9` | — | `test_mop_emergence.c` |
| T-HAR-03 | FR-HAR-001 | Unit | Random complex Matrioska: `R_H > 1e-2` | Tolerance 1 | `test_mop_emergence.c` |
| T-HAR-04 | FR-HAR-001 | Unit | Any real Matrioska, N ≥ 4: `R_H ≥ \|sin(2π/(N−1))\|`; `α₁₂ ≈ 0` refused | — | `test_mop_emergence.c` |
| T-HAR-05 | FR-HAR-002 | Unit | First Equation construction ⇒ `transported` | Return the constructed Matrioska | `test_mop_emergence.c` |
| T-HAR-06 | FR-HAR-002 | Unit | Second Equation construction ⇒ `imposed` | — | `test_mop_emergence.c` |
| T-HAR-07 | FR-HAR-002 | Unit | EQS construction ⇒ `imposed` | — | `test_mop_emergence.c` |
| T-HAR-08 | FR-HAR-004 | Unit | Constructor reconstructs from `α₁₂`; rows cancel; output labelled `assumed` | — | `test_harmony` (exists; to rename) |

### 7.7 Interfaces and outputs

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| T-API-01 | IF-API-001 | Unit | Every status has a string; `why` set on failure, untouched on success; outputs untouched on failure | Write output before failing | `test_mop.c` |
| T-API-02 | IF-API-002, IF-API-003, IF-API-004, IF-API-005, NFR-API-001 | Structural | `scripts/check_api_called.sh`: every declared function referenced in a test; no not-implemented return | Declare an untested function | `scripts/` |
| T-IN-01 | IF-JSON-001 | Integration | Valid `mop` blocks of each form load; values reach the solver | — | `test_mop.c` |
| T-IN-02 | IF-JSON-001, NFR-ROB-001 | Integration | Unknown key, wrong type, unknown component, duplicate couple, non-increasing samples: each a load error | Ignore unknown keys | `test_mop.c` |
| T-OUT-01 | IF-OUT-001, IF-OUT-002 | System | Header columns exactly as specified; every value finite | Keep `_idc` | `tests/mop_cli.sh` |
| T-OUT-02 | FR-OUT-001, FR-OUT-003, IF-OUT-003 | System | One `label.` line per column; one `harmony.` line per construction | Omit a label | `tests/mop_cli.sh` |
| T-OUT-03 | FR-OUT-002, IF-CLI-001 | System | Each refusal exits 2 with the source on stderr; load error exits 1 | Exit 0 | `tests/mop_cli.sh` |
| T-OUT-04 | FR-OUT-004 | System | Committed goldens: exit 0 and a row per model, with the allowlisted model reported as skipped and decay's sparkline falling █→▁. A golden moved by 1e-3: exit 1, `decay_model` differs, deviation 1.0e-03. Moved by 5e-7: matches, deviation 5.0e-07. A golden deleted: exit 1, named | Perturb a copy of `tests/expected/` by a known amount | `tests/sim_report.sh` |
| T-OUT-05 | FR-OUT-005 | System | The page and plots equal a fresh run (`check-engine-comparison`). A plot exists for every model; every plot is valid SVG with a title and description; the page links each one. `decay_model`'s RK4 figure equals `gia_bridge`'s `max_scaled_difference` and its expm figure is below 1e-13. `--full-precision` rounds to the default output row for row | Delete a plot; edit a figure; drop the precision flag | `scripts/engine_compare.py --check`, `tests/engine_compare.sh` |
| T-BRG-01 | FR-BRG-001 | System | `decay_model` (h = 0.025, n = 40): the reported difference equals the closed form's, for Euler `Q0(1−h)ⁿ`, RK4 `Q0 Rⁿ` and expm `Q0 e^{−1}`; a partial projection is refused, exit 2 | Compare on a shifted grid; compare a partial projection | `tests/bridge_cli.sh` |
| T-BRG-02 | FR-BRG-002 | System | `docs/results/level1_survey.md` equals a fresh `scripts/level1_survey.py` run; a sine-forced source agrees to 1e-6 under RK4; square and step forcing are refused by name, exit 2 | Regenerate and diff; bridge the sine and forced-source fixtures | `scripts/level1_survey.py`, `tests/bridge_cli.sh` |
| T-KER-01 | FR-KER-001 | System | `"incipient"` output byte-identical to `"expm"`, one notice | — | kernel regression |
| T-KER-02 | FR-KER-002 | System | Hand-computed values on a four-point series for step, linear, hold and cycle. Tank integrals equal to the trapezoid and rectangle sums (rk4 and euler); an edge-rate table. A 1000-point bisection check; thirteen malformed tables rejected, on all three parsers. The round-trip reproduces values and trajectory; the projection refuses a table by name | ASan build, leak detection on | `tests/test_forcing.c`, `tests/test_giannantoni.c` [29b] |

### 7.8 Non-functional

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| T-NUM-01 | NFR-NUM-002, FR-MOP-001 | Unit | Relative error ≤ 1e-12 vs exact series for `x = bt/a ∈ {1e-2, 1e-5, 1e-8, 1e-11}` (probe `first_equation_numerics.py`) | Naive form (5.9e-5 at 1e-11) | `test_mop.c` |
| T-NUM-02 | NFR-NUM-003 | Unit | Overflowing α returns `GIA_E_RANGE`; no NaN or Inf in any output | Return inf | `test_mop.c` |
| T-NUM-03 | NFR-NUM-005, NFR-NUM-001 | Unit | Quadrature error estimate ≥ true error; an integrable singularity beyond the subdivision limit returns `GIA_E_CONVERGENCE` | Report success | `test_mop.c` |
| T-NUM-04 | NFR-NUM-006 | Unit | Negative real and complex coordinates survive solve and output unchanged in sign | Clamp | `test_mop.c` |
| T-DET-01 | NFR-DET-001 | System | Two runs, byte-identical CSV and report | — | `tests/mop_cli.sh` |
| T-REE-01 | NFR-REE-001 | Integration | Two models on two pthreads ≡ sequential results | Shared static | `test_mop_threads.c` |
| T-REE-02 | NFR-REE-001 | Structural | `nm` on the Giannantoni objects lists no writable data symbols (`b B d D S`) | Add a static | `scripts/check_symbols.sh` |
| T-MEM-01 | NFR-MEM-001 | Integration | `test_mop` under ASan + LSan: no error, no leak | Drop a free | `make test-mop-asan` |
| T-ERR-01 | NFR-ERR-001, IF-API-001 | Structural | `nm -u` on library objects lists no `exit`, `abort`, `__assert_fail` | Call `exit` | `scripts/check_symbols.sh` |
| T-LIM-01 | NFR-LIM-001 | Integration | 65 inflows into one component, and 65 components with co-production: `GIA_E_LIMIT`, not a silent wrong answer | Truncate | `test_mop.c` |
| T-ROB-01 | NFR-ROB-001 | Unit | Committed fuzz corpus for the `mop` block: no crash, every case loads or errors | — | `tests/mop_fuzz/` |
| T-PERF-01 | NFR-PERF-001 | System | N = 64, 1,000 times, affine β, < 2 s | — | `make bench-mop` |
| T-COV-01 | NFR-COV-001 | Structural | Coverage gate ≥ 90% on the Giannantoni units; garbage input fails | Remove a test | `make coverage-check` |
| T-TRC-01 | NFR-TRC-001 | Structural | `make check-trace` passes; planting an orphan tag or an untraced requirement fails it | — | `scripts/check_trace.sh` |
| INS-SEP-01 | NFR-SEP-001 | Structural | No `#include "gssk.h"` in the Giannantoni units | Add the include | `scripts/check_symbols.sh` |
| DEM-POR-01 | NFR-POR-001, NFR-DET-002, BR-011 | Demonstration | The whole catalogue passes on macOS clang, Linux clang, Linux GCC in CI with `-Werror` | Drop a suite's CI step | `deploy.yml` matrix; `scripts/check_ci_matrix.sh` (every `CI_TESTS` suite has a step) |

### 7.9 Validation

| ID | Verifies | Level | Oracle / procedure | Mutation | Where |
|---|---|---|---|---|---|
| VAL-01 | BR-001, BR-002, BR-008, FR-IDC-010 | Validation | [09] with n = 2: 16.4 (Eq 16), 3.01 (Eq 17), 178.0 (Eq 21), 172.06 (Eq 19, τ₀ = 1.8), 15–17 cm (Eq 22); errata X4: Eq 19 τ₀ = 2 → 156.0, not 154.3; X5: net 2.61 | n = 3 | `test_mop.c` |
| VAL-02 | BR-005, BR-008, FR-MOP-006, FR-REL-001 | Validation | [23 Eq 7.1.1, 7.2, 7.3] brackets reproduced by `rel_mul(root, ref)` over 1000 seeded draws to 1e-12 (probe `eqs_relational_product.py`) | Anticommuting product | `test_mop.c` |
| VAL-03 | BR-003 | Validation | [02 Eq 3.8] `(n−1)Em`, [02 Eq 3.16–3.17] split, [02 Eq 3.9, 3.12] interaction, [02 rule 4] — through the network engine | — | `test_mop.c` |
| VAL-04 | BR-008, FR-EM-007 | Validation | [02 Eq 3.23–3.26] Fig. 3.4 totals | — | `test_mop.c` |
| VAL-05 | BR-001, BR-008, FR-IDC-006 | Validation | [10 App. Eq 28–34] worked example | — | `test_mop.c` |
| VAL-06 | BR-001, BR-008, FR-IDC-005 | Validation | [06 Eq 3.9] `±√α e^{αt}` | — | `test_duet` (exists) |
| VAL-07 | BR-007 | Validation | The three source-determined verdicts (PLAN R7): First Equation transported, Second Equation and EQS imposed | — | `test_mop_emergence.c` |
| VAL-08 | BR-010 | Validation | Every erratum X1–X12 has a test that the printed form fails its oracle | Rename an erratum's assertion | `test_mop.c`, `test_mop_emergence.c`; `scripts/check_errata.sh` |
| VAL-09 | BR-009 | Inspection | README and docs claim nothing beyond `implemented` requirements; PLAN §6 exclusions are documented | — | Review checklist; last run recorded in PLANLOG.md (`requirements-status`) |

## 8. Entry and exit criteria

- **Enter implementation of a requirement** when its ADR (if any) is merged, its catalogue tests are
  written and fail for the right reason, and the failure is quoted in the PR.
- **Mark `implemented`** when the tagged tests pass on every CI toolchain and the mutation is shown to
  fail them.
- **Release** when every `Must` requirement is `implemented`, every VAL test passes, and
  `make check-trace`, `make coverage-check` and `make ci-local` pass.
