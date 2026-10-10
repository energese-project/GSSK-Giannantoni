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
