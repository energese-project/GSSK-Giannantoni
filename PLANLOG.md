# Plan log

What changed in [PLAN.md](PLAN.md), when, and why — the reasoning behind each revision, which the plan
itself states only as conclusions. Newest first. The format follows [CHANGELOG.md](docs/CHANGELOG.md):
one entry per revision, findings before decisions.

A revision is logged when the plan's conclusions change, not for wording. Each entry says what was
found, how it was checked, and what the plan does about it, so a later reader can tell a decision
from an assumption.

---

## [Baseline amendments during implementation] — 2026-10-09

Changes to `docs/requirements/` found necessary while implementing it. They are not a plan
revision: the plan's decisions stand. Each one corrects the baseline where it contradicted itself.

### `idc-lde2`

- **T-IDC-04 (i) could not run as written.** The catalogue's case `a₁ = 0, a₀ = −t²` has roots `±t`,
  which coincide at `t = 0`. There `[[1, 1], [r₁(0), r₂(0)]]` is singular, and numerics.md N3 itself
  refuses the initial conditions. The test now uses `a₀ = −(1+t)²`, with roots `±(1+t)`. That keeps
  everything the case was for: variable roots, a zero termwise residual, and a traditional residual
  `c₁e^{φ} − c₂e^{−φ}` (here `φ = t + t²/2`). The original `a₀ = −t²` stays in the test as a refusal
  case (`GIA_E_DOMAIN`).
- **IF-API-002 gains two things it needed.**
  - `gia_lde2_solve` takes `double *t_fail`. A root collision is refused "with the time" (N3,
    T-IDC-04 (iii)), but a static `why` string cannot carry a number, and a file-scope buffer would
    break NFR-REE-001.
  - `gia_lde2_terms` is new. It returns `c`, `r(t)` and `E = e^{∫r}`, so that a test can form the
    defining equation's termwise residual itself. vv-plan.md §2 rule 2 ("residual first") requires
    this, and `f` alone cannot provide it.

### `idc-binary`

- **`gia_binary_eval` returns a `gia_status`.** It returned `void`. But `e^{u²t}` overflows for
  moderate `t`, and NFR-NUM-003 forbids returning infinity, so the function needs a way to say
  `GIA_E_RANGE`.

### `idc-riccati`

- **`gia_riccati_solve` takes `R'` and `t_fail`.** numerics.md N5 says "R' is supplied by the caller
  (a coefficient function), not differenced". The signature in icd.md had no parameter through which
  to supply it. `t_fail` is there for the same reason as in `gia_lde2_solve`: the solve runs through
  the LDE solver, so a root collision must report its time the same way.

### `idc-drift-coupled`

- **IF-API-005 gains `gia_solution_drift`.** FR-IDC-011 requires the kernel to *report* solution
  drift as exactly zero for a constant flow matrix, and to refuse it otherwise. No function in the
  interface did either.

### `emergy-source-terms`

- **IF-API-005 gains the global-balance functions and a limit check.** FR-EM-007 asks the kernel to
  "evaluate the global balance of [02 Eq 3.21] and solve it for one unknown source term", but no
  interface carried it. `gia_emergy_global_balance` and `gia_emergy_balance_solve` now do, with each
  term a `(value, weight)` pair. NFR-LIM-001 needs the refusal to be observable on its own, so
  `gia_emergy_check_limits` reports `GIA_E_LIMIT`.

### `emergy-ordinal-forms`

- **IF-API-005 gains the ordinal forms and the circle product.** FR-EM-005 and FR-EM-006 had no
  interface. A form records its shape as well as its entries: a binary is a column, a duet a row,
  a duet-binary is 2 × 2. It needs the shape because the circle product of [06b Eq 2] is defined
  as a column ∘ a row.

### `mop-first-equation`

- **IF-API-004 gains `gia_quad_gk15`.** T-NUM-03 requires that "an integrable singularity beyond the
  subdivision limit returns `GIA_E_CONVERGENCE`". Through `gia_mop_couple` this cannot happen: a
  piecewise-linear β that avoids 0 gives a bounded, continuous integrand. Exposing the N2 integrator
  makes its error estimate and its refusal testable directly.
- **`mop.h` does not include `relational.h` yet**, and FR-MOP-007 (relational-valued couples) is
  deferred to `mop-relational-algebra`, which creates that header. The seed's `mop` block, the MOP
  CSV and the CLI exit codes (IF-JSON-001, IF-OUT-002, IF-CLI-001, FR-OUT-002, NFR-ROB-001) follow in
  `mop-seed-cli`, which splits from this task to keep each PR reviewable.

### `mop-relational-algebra`

- **`rel_exp` returns a `gia_status`.** It returned an `rel_t`. But `eᵃ` overflows, and
  NFR-NUM-003 forbids returning infinity, so the function needs a status to report the overflow.
- **IF-API-004 gains `gia_mop_couple_rel`** for FR-MOP-007. The requirement is "solve the First
  Equation componentwise on relational elements for k = 1, refuse k > 1", and no interface carried
  it.

### `mop-second-equation`

- **The ordinal power in [23 Eq 6.3] is read as the principal root of unity.** FR-MOP-005 writes
  `A(t) = α₁₂(0) ∘ r + ln(c₁ + c₂t)`, and Eq 6.3 writes `r` as `({ᴺ⁻¹√1})^{↑N N}`. No source defines
  that ordinal power. The engine takes `r = w = e^{2πi/(N−1)}`, the principal (N−1)-th root of
  unity, and the roots of Eq 6.1 as `w^{j−2}`, j = 2…N. `mop.h` says the reading is the kernel's,
  and T-MOP-07 pins it, so a different reading is a visible change rather than a silent one.
- **T-MOP-07's oracle is reconstructed.** No source writes Eq 4.2 explicitly. The time dependence of
  the printed solution is `ln(c₁ + c₂t)` alone, so `u = Ȧ` satisfies `u' + u² = 0`, and that is
  what the test checks, by central differences.
- **`e^B` is computed in closed form.** `B = A·M` with `M = [[1, −1], [−1, 1]]` and `M² = 2M`, so
  `e^B = I + (e^{2A} − 1)/2 · M` exactly. Overflow of `e^{2A}` is `GIA_E_RANGE` (NFR-NUM-003).

### `mop-seed-cli`

- **`k` is required.** IF-JSON-001's grammar marks `reference`, `second_equation` and `eqs` optional,
  and does not mark `k` optional. A block without `k` is a load error, rather than one that silently
  takes a default.
- **A key given twice is a load error.** cJSON keeps both members, and a lookup returns the first
  one. So without this check, `{"k": 1, "k": 2}` would load as `k = 1` and nobody would be told.
- **The loader adds a `detail` buffer to the `why` convention.** IF-API-001's `why` is a static
  string, so it cannot name the offending key. `gia_mop_seed_load` therefore also takes
  `char *detail, size_t cap`, owned by the caller, for the JSON path (`mop.beta[2].samples[1]`).
  `why` stays static.
- **Exit status 2 covers every refusal, not only `GIA_E_UNSUPPORTED` and `GIA_E_DOMAIN`.**
  `GIA_E_RANGE`, `GIA_E_CONVERGENCE` and `GIA_E_LIMIT` are also the kernel declining to compute
  (NFR-NUM-003: never a non-finite value). They are not load or validation errors, and exiting 1
  for them would misreport them as such.
- **IF-OUT-002 stays `planned`.** The CSV has `time` and the couple columns. Its last column, `R_H`,
  is the harmony residual of FR-HAR-001, which lands with `mop-harmony-detector`. That task adds the
  column and flips IF-OUT-002.
- **`eqs` is checked, not evaluated.** [23 Eq 7.1–7.3] need the reference couple's relational
  coordinates `{Σ₀, Φ₀, Θ₀}`. The seed does not supply them, and the First Equation's complex `α₁₂`
  does not determine them. The CLI checks the parameters through `gia_eqs`, so an `ε₂ ≠ ε₃` seed is
  refused (X11, exit 2), and says that no coordinates are computed.
- **`"beta": "network"` loads and is refused at solve time** (`GIA_E_UNSUPPORTED`, naming
  FR-MOP-008) until `mop-network-beta` lands. The grammar accepts it, so treating it as a load error
  would be wrong.
- **T-PERF-01's benchmark lives in `tests/bench_mop.c`**, not `bench/`. `check-trace` looks for
  `Verifies:` tags only under `tests/` and `scripts/`. It measures CPU time, so that another job on
  a loaded runner cannot fail the bound. It is not in CI: NFR-PERF-001 is a Could.

### `mop-ordinality`

- **"Feed one interaction module" means a direct, quantity-carrying leg into it.** ADR 0021 says
  "along quantity-carrying legs" but does not say how far. A direct leg is the reading that keeps
  ADR 0021 §3's hand-derived `input.json` row (`store_1`'s read control does not feed, `consumer_1`'s
  drawn one does) and does not double-count a component upstream of another. A "replicating
  process" is any node with `output_mode: replicate` legs to both components.
- **Closure counts components only.** ADR 0021 §3 gives `closed_loop.json` closure 1.000, so the
  heat sink is not counted. The old fraction counted every non-module node. Three assertions in
  `test_giannantoni.c` pinned the old numbers. ADR 0021's Consequences require revising them in
  this PR, and they now assert closure 1, with the reason in a comment:
  `test_ordinality_invariant_under_respelling` (2/3 → 1, below → at maximum) and
  `test_sink_is_never_closed` (2/3 → 1, 4/5 → 1). The latter is retired with ADR 0015's step in
  `mop-generative-empower`.
- **FR-ORD-004 stays `planned`.** The report labels closure a proxy. But `gia_generate`, the ADR 0015
  heuristic, still chooses what to close from the cycle scan, so closure is not yet "deciding
  nothing". That changes when `mop-generative-empower` replaces the step.
- **IF-OUT-003 stays `planned`.** Its `ordinality`, `maximum_ordinality` and `closure (proxy)`
  lines are written now. Its `harmony.<construction>: <verdict>` lines need FR-OUT-003 (W8).

### `mop-generative-empower`

- **"Total empower" is the sum of the components' empower at `t_end`.** ADR 0021 §2 says "the total
  empower at `t_end` with the emergy pass" but does not say over what. The sum leaves habitat out,
  for the reason ADR 0021 §1 gives: a sink's accumulation is outside the system [10]. It is
  documented in `engine.h`.
- **A tie is within 1e-9 relative, not exact equality.** The emergy pass sums in an order that
  depends on the file. Exact comparison would let rounding break a tie differently when the nodes
  are reordered, which would violate T-ORD-06's byte-identical output.
- **Two refusals the ADR does not spell out.** The step returns the seed unchanged, and says so, in
  two cases: when the emergy pass cannot evaluate any candidate (NFR-LIM-001), and when
  `#sources + #sinks` additions have not reached Maximum Ordinality. The second cannot happen by the
  ADR's termination argument; it is a guard, not a behaviour.
- **The added pathway carries `"flow_type": "max_empower_pathway"`** so a reader of the output can
  see which pathways the step added. The loader ignores `flow_type`.
- **vv-plan §5 dispositions carried out.**
  - `test_emergent_quality_closes`, `test_sink_is_never_closed` and
    `test_emergent_quality_is_a_component` are retired.
  - `test_generative` and `test_control_does_not_close_a_pathway` keep their mode assertions. Their
    ADR 0015 ones ("one component emerged", `emerged_from`) are revised to ADR 0021's: no component
    is added, and the added pathway is named.
  - `test_emergence_is_order_invariant` (T-ORD-06) now compares the appended pathways.
  - `test_mop_threads.c` detects growth by the new pathway label.
- **`label.generative_step` becomes `implemented`** (E7), and `docs/giannantoni_assessment.md`'s rows
  E6–E8 are updated to match.

### `kernel-method-label`

- **T-KER-01 is `tests/kernel_method.sh` (`make test-kernel-method`, in CI).** The catalogue calls it
  a "kernel regression". The golden-file regression suite cannot express "two runs are
  byte-identical, and one prints a notice", and ADR 0018 forbids adding golden files. Two examples
  (`decay_model`, `diffusion_model`) are run under both spellings. The test also checks the
  round-trip that ADR 0022 decision 4 promises.
- **An unknown `method` string still silently becomes `auto`.** This is unchanged, and out of this
  task's scope. Note that before this change `"expm"` fell into exactly that path: the red run's
  serialisation check caught it.

### `mop-network-beta`

- **"Carried from component i to component j" means a direct pathway.** R8 does not say how far
  emergy is followed. A pathway through a module, or through intermediate components, would need a
  rule for splitting the module's output emergy among its products. The sources give no such rule,
  and the engine already has one for direct pathways (the pass's partition and replicate rules). So
  a couple is related exactly when a direct pathway `i → j` moves quantity: not a used leg, and not
  a read control. Several such pathways sum.
- **A loop-closing pathway carries what rides on it.** The emergy pass zeroes back edges in its own
  accounting, so that a loop does not create emergy. Applying the same rule here would make one
  pathway of every 2-cycle "unrelated", which contradicts R8's "couples with no pathway". So
  `gia_emergy_carried` (added to IF-API-005) shares each origin's empower across all of its
  quantity-carrying pathways.
- **β by finite difference.** `β = α' = E'/E`. The derivative is a central difference with
  `h = 1e-5·max(1, t)`, or a second-order forward difference at `t < h`. T-MOP-10 checks it by its
  defining equation, `∫β = Δα`, to 1e-6. A pathway carrying no emergy has no logarithm, so it is
  `GIA_E_RANGE`.
- **The network form refuses k ≠ 1** (`GIA_E_UNSUPPORTED`, exit 2). FR-MOP-008 defines β with
  k = 1.
- **The per-seed verdict table** (PLAN §7, T-HARM-2b) needs the harmony detector. It lands with
  `mop-harmony-detector`.

### `mop-harmony-detector`

- **The §6 perturbation scales each couple; it does not add to it.** Read as an added constant,
  `β_h + δξ` changes the *shape* of an affine-power β. The First Equation's α is β's integral, and
  under that reading the ratios of α would no longer equal the ratios of β. The First Equation
  would then read `absent`, which contradicts R7's derived verdict (`transported`). So
  `β_δ = β_h·(1 + δξ)` scales `a` and `b` together. ξ is the §6 LCG's pair in `[−1, 1]²`, and the
  input's own residual is taken from `β(1)`.
- **The EQS construction reads the root as a complex number.** [23 Eq 7.1–7.5] give relational
  coordinates, while [23 Eq 5.6.5] compares complex ratios. With `ref = {1, 0, 0}`, `ψ₁ = 1` and
  `A = 1`, `gia_eqs` returns `ρ = e^{E B}` and `φ = E C`. From these the construction recovers the
  root `B + C j + C k` and writes it as `B + i√2 C = e^{i√2ψ_l}`. The map preserves the modulus.
  With `ψ₂ = 1` and `ε = −2π`, `√2ψ_l = 2π(l−1)/(N−1)`. No source gives this mapping. It is the
  construction under which [23 §8 ii]'s "assumed harmony" can be put to the detector at all, and it
  reads `imposed` as R7 says.
- **T-HAR-01's symbol rule is "no undefined reference".** The detector is itself named
  `gia_harmony_*`, so "nm mop.o lists no `gia_harmony_*`" cannot be meant literally. The rule
  `check_symbols.sh` enforces is that mop, mop_seed, relational, idc and engine have no undefined
  reference to `gia_harmony_assume_*`. `validation.c`'s `gia_validate_harmony` checks the
  constructor's invariants and does call it, so it is not on that list. `bin/test_mop_emergence`
  is linked from `mop.o` and `relational.o` alone.
- **What the run report evaluates.** It evaluates the three constructions at N = the seed's
  component count, and only for N ≥ 3. It also observes the network's Matrioska on the reference
  row, where a row couple with no pathway makes it `absent`.
  `docs/results/harmony_verdicts.md` records each example's lines. `tests/mop_cli.sh` reads the
  table and checks the report against it. Both examples' networks are "not evaluated", because
  their source has no `quality_input`.
- **IF-API-004 now lists what was built.** `gia_mop_network_beta(…, gia_beta *)` became
  `gia_mop_network(…, gia_matrioska *alpha, gia_matrioska *beta)`, in `mop_seed.h`
  (`mop-network-beta`). The three constructions are added. `api_called_selftest.sh` counts 33.
- **The MOP CSV's `R_H` column** is the residual of the reference row at each time. The cell is
  empty where the residual is undefined: N < 3, a row couple unrelated, or `α₁₂ = 0` at `t = 0`.

### `requirements-status`

After W0–W8 merged (main at `5a0748a`), each BR row of stakeholder.md was checked against its own
acceptance line. Each requirement still `planned` was checked against its verification. All 88 are
now `implemented`. What the check found, and what was changed so that it holds:

- **Two suites were in `CI_TESTS` but not in `deploy.yml`:** `test-mop-emergence` and
  `test-kernel-method`. `test-coverage-gate`, from W1, was in the same state. The Makefile says
  "add a suite here when you add its step to deploy.yml", and the step had not been added. All
  three ran, and passed, locally and on `main`, but CI never ran them. They are now steps.
  `scripts/check_ci_matrix.sh` (DEM-POR-01, in CI) now fails on any `CI_TESTS` suite without a
  step. `test-advanced` stays exempt, as the Makefile records. The check also fails if a toolchain
  leaves the matrix or `-Werror` leaves `CFLAGS`.
- **VAL-08 had no check, and two errata had no named test.**
  - X6 (the co-production/interaction exponents) is now asserted in T-ORD-01: a co-production
    counts as ½, not as an interaction.
  - X9 (the literal j = 1 term of [23 Eq 5.6.5]) is now asserted in `test_mop_emergence.c`.
  - The X1 and X3 assertions now carry their erratum's name.
  - `scripts/check_errata.sh` fails on any erratum in PLAN §5 that no assertion names.
- **IF-API-002 and IF-API-004 were never flipped**, though every function they declare is built
  and called (`check-api-called`). Flipped.
- **VAL-09, the inspection.** README.md, AGENTS.md and `docs/*.md` were read for claims beyond
  `implemented` requirements. None was found. The README instead *under*-stated: it said the
  engine "does not yet implement IDC or the MOP", called ordinality a proxy, and said the generative
  step adds components. It now lists what is implemented, and says the network trajectories remain
  the classical matrix exponential and the harmony constructor is `assumed`. PLAN §6's exclusions
  stay documented in PLAN.md and `docs/giannantoni_assessment.md`.
- **Tags.** BR-004, BR-006, BR-009 and BR-011 had no `Verifies:` tag. They are now tagged on the
  tests their verification line names: T-MOP-01, T-MOP-10, T-ORD-02, T-ORD-04, T-OUT-02, T-OUT-03
  and T-REE-01.

### `engine-bridge`

- **FR-BRG-001 and T-BRG-01 are new.** The baseline had no requirement for running one model
  through both engines. TODO.md §10.5 asked for one, and its blocker (the projection and coverage
  report) had already landed. The bridge is a separate binary, `bin/gia_bridge`, and the only one
  that includes both `gssk.h` and `engine.h`. Neither engine's library units include the other's
  header, so ADR 0011 and NFR-SEP-001 still hold.
- **What it measures.** Both trajectories are classical: the Giannantoni engine solves a network by
  `exp(A t)` (E3). The difference is therefore the kernel integrator's error against the exact
  exponential, labelled `classical`. The TODO item hoped to measure "the drift critique against a
  real integrator". That needs an incipient network solution, which neither the code nor, so far,
  the sources provide. The item says so.
- **Partial projections are refused**, exit 2, and the coverage report is printed. A projection
  that drops a pathway changes the equations of the nodes it touched, so even those nodes' "carried"
  values would differ by the missing pathway. TODO.md's own blocker note warned of exactly this.
- **The oracle is the closed form of each recursion.** On `decay_model`, Euler is `Q0(1−h)ⁿ` and RK4
  `Q0 Rⁿ` exactly, so the reported difference is fixed by arithmetic. The tolerance comes from
  cancellation: the RK4 difference (3.3e-9 out of 0.37) is good to about 1e-5 relative.

### `level1-survey`

- **FR-BRG-002 and T-BRG-02 are new.** TODO.md asked how often the exponential form holds for real
  Odum graphs ("the actual Level 1 claim"), and the baseline had no requirement for answering it.
  `scripts/level1_survey.py` writes `docs/results/level1_survey.md`, and `make check-level1` fails
  the quality gate when the committed table no longer matches a fresh run. It runs only there
  (Linux gcc), not in CI_TESTS, because `supply_chain_30` alone takes about 30 s a bridge run.
- **The survey found a silent projection defect, now fixed.** The projection read the engine's
  forcing key (`kind`) from the kernel's forcing object (`waveform`), found none, and dropped the
  forcing while reporting `forced_source_model` as 100% carried. The engines then disagreed by
  99.8%. The projection now translates sine, ramp (from t_on ≤ t_start), exponential and an
  unbinding clamp exactly, onto the engine's clock (t_start → 0). It refuses square, step,
  sawtooth, impulse, jitter, a late ramp and a binding clamp by name. Before the fix, the red run
  of `[29b]` failed all six of its assertions.
- **The bridge compared a forced source's base value** with the kernel's forced value. It now
  applies `gia_forcing_value` to held nodes before comparing.
- **A scaled metric was added.** The relative difference read 1.00 wherever both engines hold a
  node near zero (Euler on `diffusion_model`), which is noise, not error. The bridge now also
  reports `max_scaled_difference`, the largest absolute difference divided by the node's largest
  magnitude over the run, and the survey tabulates that.
- **The kernel's `expm` freezes forcing over a step.** It is exact only for an autonomous model. On
  the sine-forced fixture at dt 0.5 it is off by about 0.4%, where RK4 is off by about 1e-9.
  `docs/results/level1_survey.md` says so beside the numbers.
- `invalid_model.json` is a negative fixture, so the survey excludes it.

### `h8c-data-driven-forcing`

- **FR-KER-002 and T-KER-02 are new.** TODO.md §7.1 listed the item, and ADR 0006 deferred it.
  Schema v5 §7 had already fixed the shape, so the kernel implements that shape and adds nothing
  to it. Two defaults the spec left open are now set: `interpolation` defaults to `linear` and
  `extrapolation` to `hold`. ADR 0006 gains an addendum rather than a new ADR, because a table is
  one more waveform in the same two places.
- **Decisions the spec did not make:**
  - `times` are absolute, and `t_on` on a table is rejected.
  - Under `cycle`, the last knot is the start of the next pass.
  - A step knot is closed on the right, as `step` is at `t_on`.
  - Table keys on any other waveform are rejected.
  - A one-point table cannot cycle.
- **Red run.** The tests ran with only the enum added to the header. The kernel refused the first
  table model with "unknown key 'times'", and the suite stopped at its FATAL load.
- **Mutations**, each run under an ASan + UBSan build with leak detection on. All nine were caught.

  | Mutation | Caught by |
  |---|---|
  | Knot open on the right | 2 FAILs |
  | No cycle fold | 7 FAILs |
  | No strict-order check | 5 FAILs |
  | `GSSK_AddNode` storage reject leaks the table | LSan leak |
  | `GSSK_Free` misses node tables | LSan leak |
  | `GSSK_AddEdge` reject leaks the table | LSan leak |
  | Serialiser drops `interpolation` | 12 FAILs |
  | Table keys accepted on a sine | 1 FAIL |
  | Linear evaluated as step | 14 FAILs |

- **Baseline amendment.** There is one new pinned trajectory,
  `tests/expected/tabulated_source_model.csv`, for the new example. `make test-update`
  regenerated every baseline, and every existing one came out byte-identical, so the change moved
  no prior trajectory. The values the CSV pins are checked independently in `test_forcing.c`
  against hand-computed integrals: the trapezoid sum under RK4, and the rectangle sum under Euler.
- `docs/results/level1_survey.md` was regenerated for the new example. A table is a forcing the
  projection refuses, so the model is not carried whole.
- **crux.** This container has no `crux` binary and no crux MCP server, so the task could not be
  marked done here. The maintainer should close `h8c-data-driven-forcing` once this merges.

### `sim-report`

- **FR-OUT-004 and T-OUT-04 are new**, at the maintainer's request (2026-10-10): CI should show
  the simulation tests' results in a readable form, not only PASSED. `scripts/sim_report.py` runs
  every example, decides each verdict with `bin/csv_compare`, so the report and the gate cannot
  disagree, and writes markdown for the job summary.
- **What it can and cannot claim.** A golden file pins a trajectory; it does not say that the
  trajectory is right. The report says so, and points to the two checks that do: the hand-computed
  suites, and the Level 1 survey against `exp(A t)`. The quality gate now puts the survey on the
  summary page too.
- **Self-test oracle:** a copy of `tests/expected/` with one value moved by a known amount (1e-3
  must fail and print 1.0e-03; 5e-7 must pass and print 5.0e-07), plus a deleted golden. There are
  two mutations: a comparator that always passes (4 FAILs) and a sparkline that is always flat
  (1 FAIL). Both were caught.
- **Not done here.** Several example models state invariants in their descriptions, such as
  `archetype_price_per_instance`'s price ratios and conservation of money. The report quotes them,
  but checking them mechanically would need a per-model assertion format, and that is a separate
  task.

---

## [Revision 3] — 2026-10-09

A requirements baseline, written from a senior C developer's view of what was missing. The plan
said what to build and in what order. It did not say what the delivered kernel must do, under what
conditions, to what accuracy, or how each claim is proven. In V-model terms, the left arm had been
skipped.

### Added

- `docs/requirements/`:
  - stakeholder and business requirements with acceptance criteria;
  - the SRS (functional and non-functional);
  - an interface control document: C API contracts, the seed's `mop` block, outputs, CLI;
  - a numerical design per algorithm;
  - a V&V plan, with oracle, tolerance and determinism policy and an 81-entry test catalogue.

  In all, 88 requirements, each traced to a source equation and to its verification.
- `make check-trace` (`scripts/check_trace.sh`, POSIX sh + awk) runs in CI. It fails on an untraced
  requirement, an orphan `Verifies:` tag, or an `implemented` claim with no tagged test. Nine planted
  faults each fail it; the unmodified tree passes under macOS awk and Linux mawk. Ten existing tests
  are tagged, which backs the ten requirements already `implemented`.

### Found — writing requirements surfaced what the plan did not

- **The incipient derivative is two operations.** Pointwise, `(f'/f)ⁿ f`, it is the sources'
  definition, and it is not additive. The superposition the sources print as the solution of the
  second-order equation [06 Eq 3.6] solves it only if the derivative acts termwise on exponential
  terms (residual 0 termwise, −1.67 pointwise). The sources use termwise for every superposition they
  print. Decided as R16. It is also the reason E4b holds.
- **Two different drifts had been conflated.** Solution drift of a known model [06 §4] is zero for
  constant coefficients, and the sources define none for coupled nonlinear systems.
  Output-projection drift [09 Eq 13] is defined on any trajectory. Revision 2's T-IDC-11 expected
  the first to be nonzero for an interaction network. The SRS now has one requirement for each, and
  separate CSV columns.
- **[06 Eq 3.7] is wrong under every reading** (erratum X12). Its "second solution at a double root"
  leaves residual −0.7 classically, −0.7 termwise and 0.93 pointwise. Under the termwise derivative a
  double root has only the one-parameter family `c·e^{∫α̃}`, so ICs inconsistent with it are now
  refused.
- **[23 Eq A2.6] gives the real unit i a 4π period** (erratum X11), and Eq 7.4's `√2ψ` holds only for
  ε₂ = ε₃. The EQS uses i's factor only as a scale; unequal angles are refused.
- **The First Equation's closed form needs a numerically stable formulation.** Evaluated as printed
  (even corrected) it cancels near t = 0: relative error 5.9e-5 at t = 1e-11. The unified
  `log1p`/`expm1` form holds 4e-16 and is continuous through b = 0, which also removes the b = 0
  special case (probe `first_equation_numerics.py`).
- **Two defects in the existing engine, by inspection:**
  - `src/engine.c:2830` holds `static const gia_model *sort_model` for a `qsort` comparator. That is
    file-scope mutable state, which `AGENTS.md` forbids, and it makes the generative step
    non-reentrant.
  - `combine_inflows` and `gia_emergy_at` silently ignore inflows beyond 64 and co-production
    masks beyond 64 components, giving a wrong answer instead of an error.

  Both are now requirements (NFR-REE-001, NFR-LIM-001), with tests and a task (`guard-reentrancy`).

### Decided

- The baseline is the specification; the plan is the schedule. Where they disagree, the baseline wins.
- PLAN §8 moved into the V&V catalogue, so there is one test list, not two that drift apart.
- New API surface uses a `gia_status` code with a `why` out-parameter naming the source of any
  refusal. The library never prints or exits.
- The Giannantoni engine is native only in this baseline; WASM stays the GSSK kernel's surface.

---

## [Revision 2] — 2026-10-09

Re-read against the full source set the maintainer added to `docs/`: the 2002 book that began the
programme (scanned; read by OCR and page images), Giannantoni 2006 "Emergy Analysis as the First
Ordinal Theory", Giannantoni & Zoli 2009, both printings of the 2010 MOP paper, the 2010 protein-folding
paper, and the 2022 JAMP paper. Prompted by two questions: whether `AGENTS.md` or `CLAUDE.md` lets an
agent avoid implementing the mathematics, and whether revision 1's open questions survive the earlier
and later sources.

### Found — loopholes in the repository's own process

Nothing in `AGENTS.md` or `CLAUDE.md` forbids implementing IDC or the MOP. Several mechanisms let a
build go green without doing it:

- `make test-update` regenerates golden CSVs from the current binary, and `AGENTS.md` documents that
  as the way to add a test — the implementation becomes its own oracle.
- `make test` prints `SKIPPED` and passes for a model with no golden file (`Makefile:231`).
- The coverage gate is 35% (the README badge says ≥85%), covers only `gssk.c` and `advanced.c` — not
  the Giannantoni engine — and prints `OK` when the percentage cannot be parsed, which is always on
  macOS because BSD `grep` has no `-P`.
- `AGENTS.md` sanctions `(void)` stubs for unimplemented features.
- `TODO.md` and `AGENTS.md` describe the matrix exponential as IDC. This is how the kernel's
  "incipient" solver came to be classical numerics under an IDC label.
- ADRs 0014 and 0015 define ordinality as cycle coverage, which is not Giannantoni's definition;
  redefining a term made the tests pass.
- Harmony tests check a construction against itself.
- Revision 1 of this plan had escape hatches of its own: ADRs allowed to conclude "cannot be settled",
  harmony tests not built, a verdict pinned "whatever it is".

### Found — the sources, checked

- The incipient derivative `(d̃/dt)ⁿ f = (f'/f)ⁿ f` is stable from the 2002 book to 2023.
- The relational product table [23 Eq 5.1.3–5.1.5] is commutative and non-associative as printed —
  and it is the table *as printed*, applied to a De Moivre root, that reproduces the EQS formulas
  [23 Eq 7.1–7.3], to 8.9 × 10⁻¹⁶ over 1000 random draws. The associative reading fails Eq 7.3.
  So the table is not a typo, and revision 1's worry about it is settled.
- Under that same table the EQS roots are *not* roots of unity (N = 4, l = 1 cubes to
  (0.5407, 0, −0.6657)). They are roots of unity only by De Moivre angle arithmetic. The papers use
  two different operations and the code must name both.
- The order of m and n in Ordinality `(m n)` is reversed in the 2023 wording relative to 2022 and to
  the exponents themselves (½ co-production, 2 interaction, 2/2 feedback). The mathematics fixes it:
  m counts interactions, n counts co-productions.
- At maximum ordinality "all the various couples … become of Ordinality {2/2}" [22 §12.1] — every pair
  of components in mutual feedback — which is exactly strong connectivity of the component graph.
  Today's cycle-coverage measure says two disjoint feedback loops are at maximum; they are not.
- The 2002 book's Maximum Em-Power Principle [02 Eq 5.3] gives the generative step an objective —
  maximise empower — in place of a graph heuristic.
- The Second Fundamental Equation is never written explicitly in any source. Its printed solution
  [23 Eq 6.3] depends on time only through `ln(c₁ + c₂t)`, whose derivative is the general solution of
  the Riccati `u' + u² = 0` — the "Riccati of ordinal nature" the paper names. That is now the oracle.
- **Harmony does not emerge from any equation in the sources.** The First Equation transports it from
  its inputs (it is decoupled per couple and homogeneous of degree 1). The Second Equation's solution
  writes the roots of unity in explicitly. The EQS says it assumes them. The "Diffusive Generativity"
  that the papers credit with producing harmony is never formalised, and the 2023 appendix calls
  harmony "not a necessary consequence".
- Giannantoni & Zoli 2009 is reproducible: its incipient Taylor series with n = 2 gives 16.4, 3.01,
  178.0, 172.06 and the 15–17 cm range exactly as printed.

### Found — errata

Ten, each decided by a probe in `docs/sources/probes/` (register in PLAN §5). The material ones:

- [23 Eq 5.5.7–5.5.8] do not satisfy their own defining equation for k ≠ 1 (residual 1.25 and −0.625
  at k = 2, b = 0.25, t = 1; the derived form gives 0). A k = 1 test would never have caught it.
- [06 Eq 3.17] inverts the Riccati linearisation (residual 929 against 1.3 × 10⁻⁸ for the standard
  substitution).
- [06 Eq 3.22] cannot be derived from Eq 3.19, and adds a pure number to a rate.
- [09 Eq 19] at τ₀ = 2 gives 156.0, not the printed 154.3; the minimum scenario's "net 1.91 °C"
  subtracts a different baseline from the maximum scenario's.

### Decided

- Every revision 1 open question and every §8 maintainer decision is closed (PLAN §4, R1–R15).
- What the sources do not define is listed as out of scope with the reason and the engine's behaviour
  (PLAN §6): harmony emergence, the solar-system and Mercury results (parameters not given), the
  three-body and N-body solutions (in untranslated 2007–08 books), the Riccati duet and Abel n-et,
  network phases.
- Guardrails G1–G7 close the loopholes above (PLAN §1, workstream W1). Departures from a printed
  equation require a probe, a register entry, and a test that the printed form fails.
- Sources: text and page renders of the two CC BY 4.0 papers are committed under
  `docs/sources/cc-by/`. Those of works not openly licensed stay in the gitignored
  `docs/sources/local/`, for the same reason the PDFs are gitignored in a public repository.

### Corrected

- Revision 1 said the incipient derivative in `gia_idc_amplitude` was "a definition, not an
  operator", implying a shortcut. It is Giannantoni's own definition in every source from 2002 on.
  The fault in the drift columns is their invented inputs, not the formula.

---

## [Revision 1] — 2026-10-09

First plan, written from the 2006 JCAM paper and the 2023 JAMP paper only.

### Found

The repository's claim to "implement" IDC and the MOP was an overclaim:

- The kernel's `"incipient"` method and the MOP engine's trajectories are a Padé (3,3) matrix
  exponential — classical numerics. For constant coefficients IDC and TDC coincide anyway [06 §4].
- The `_idc`, `_tdc` and `_drift` CSV columns apply a correct identity to per-node φ(t) set by
  invented rules (e.g. storage rate = level/capacity), disconnected from the simulated trajectory.
- The "Harmony Relationships" matrix is built from roots of unity and then checked against its own
  construction.
- "Ordinality" is the fraction of components on a cycle, not Giannantoni's `{k, (m n)}`.
- The generative step is a graph heuristic, and the First and Second Fundamental Equations are absent.
- Odum's emergy algebra (co-production, split, no double counting) is implemented properly — but it is
  Odum's, computed as static accounting, not Giannantoni's incipient formulation of it.

Also found: if drift were "fixed" by taking `φ = ln Q`, it would report nonzero drift for any
constant-coefficient network with more than one mode, contradicting [06 §4 (i)]; the incipient
derivative acts on exponential components, not on the log of their sum.

And a proposition that shaped every harmony test since: the First Fundamental Equation is decoupled
per couple and homogeneous of degree 1, so `α_1j = ω α_12` holds exactly when `β_1j = ω β_12`.
Solving it can only carry harmony through from the inputs; it cannot create it.

### Decided

Remediate the labels first, then implement the First Fundamental Equation with the defining equation's
residual as the oracle, and test harmony with positive and negative controls so that a "no" is
possible. Open questions were left for later ADRs; revision 2 closes them.
