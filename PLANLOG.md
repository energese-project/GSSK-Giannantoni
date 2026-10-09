# Plan log

What changed in [PLAN.md](PLAN.md), when, and why — the reasoning behind each revision, which the plan
itself states only as conclusions. Newest first. The format follows [CHANGELOG.md](docs/CHANGELOG.md):
one entry per revision, findings before decisions.

A revision is logged when the plan's conclusions change, not for wording. Each entry says what was
found, how it was checked, and what the plan does about it, so a later reader can tell a decision
from an assumption.

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
