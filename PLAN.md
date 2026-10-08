# PLAN — Giannantoni's IDC and MOP, implemented and tested from the sources

**Status:** proposed · **Revision:** 2 (2026-10-09) · **Scope:** `bin/giannantoni_sim`, the
translation units it gains, the claims made about it, and the repository guardrails that decide
whether a "green" build means anything.

Revision 2 re-reads the full source set: the 2002 book that began the programme, and the 2006, 2009,
2010, 2022 and 2023 papers. **Every question revision 1 left open is now decided** (§4), each either
by an equation in a source or by a probe whose output is quoted. What the sources do not define is
listed in §6 as out of scope, with the reason and what the engine does instead. Nothing is left for a
later ADR to "settle"; the ADRs in §9 record decisions already made here, as the repo's ADR-first
practice requires.

This is a **test-first** plan. Every task starts by writing the tests in §8 that cite its source
equation, running them, and recording the red output, before any implementation.

Citation keys (`[02]`, `[06]`, `[06b]`, `[09]`, `[10]`, `[10p]`, `[22]`, `[23]`) are defined in
[docs/sources/README.md](docs/sources/README.md), along with the extracted text, page renders and
probe scripts.

---

## 1. Loopholes: what lets an agent avoid the real work

Nothing in `AGENTS.md` or `CLAUDE.md` forbids implementing Giannantoni's mathematics. Several
mechanisms do, however, let a build go green without implementing it. Each was checked at `4e96569`:

| # | Loophole | Where | Effect |
|---|---|---|---|
| B1 | A regression model with no golden file prints `SKIPPED` and exits 0 | `Makefile:231` | A new model "passes" by having no expected output |
| B2 | Golden files are regenerated from the current binary, and the documented way to add a test is to run that | `Makefile:257`, `AGENTS.md` §Testing | The implementation becomes its own oracle: whatever it outputs is "expected" |
| B3 | The coverage gate is 35% (README badge says ≥85%), measures only `gssk.c`/`advanced.c` (`SOURCES`, `Makefile:40`), not `engine.c`, and prints `OK` when it cannot parse the percentage. `grep -oP` does not exist in macOS grep, so locally it always parses nothing | `Makefile:599, 623-628` | The Giannantoni engine has no coverage gate at all |
| B4 | "If a feature isn't implemented, return a proper error code or use `(void)` for intentional stubs" | `AGENTS.md` §Coding Standards 3 | A sanctioned stub: compiles warning-free and passes every test that does not call it |
| B5 | The repo's own docs equate IDC with the matrix exponential: "stiff in IDC-eligible form, the matrix exponential already handles it" (TODO Non-Goals); "Incipient-calculus solver as baseline" (TODO Vision); GSSK "using Euler or RK4" (AGENTS Overview) | `TODO.md:20-30`, `AGENTS.md:6` | An agent reading these believes `expm` *is* IDC. This is how E2 (§2) arose |
| B6 | "Enforce physical conservation (clamping Q < 0 to 0.0)" | `AGENTS.md` §Fail-Safe | Applied to Relational Space it silently destroys signed and complex coordinates, and it contradicts the non-conservative algebra the MOP rests on `[10 Eq 17.1]` |
| B7 | "Declared-lossy projection" and "refuse with a named finding" | ADR 0011, `TODO.md` §10.5 | Legitimate, but unbounded: an unimplemented law becomes a *reported* one |
| B8 | "Done" means merged with green CI | `CLAUDE.md` | Combined with B1–B3, green is cheap |
| B9 | **Revision 1 of this plan** let ADRs conclude "cannot be settled from the sources", left T-HARM-4/5 unbuilt, and pinned a harmony verdict "whatever it is" | `PLAN.md` r1 §6–7 | An escape hatch in the plan itself. Closed in this revision |
| B10 | An ADR can redefine a term so the tests pass. ADR 0014/0015 define "ordinality" as cycle coverage, which is not Giannantoni's definition | `docs/adr/0014`, `0015` | Moving the definition replaces doing the mathematics |
| B11 | Tests that check a construction against itself (harmony rows sum to zero *because* they were built from roots of unity) | `tests/test_giannantoni.c:168-205` | Tautologies count as passing tests |
| B12 | A suite runs in CI only if added by hand to `CI_TESTS` and `deploy.yml` | `Makefile:884`, `AGENTS.md` | A suite missing from that list is never run by CI |

**Guardrails this plan adds** (tasks in W1, §7). These rules apply to every Giannantoni task:

- **G1. No golden files for IDC or MOP.** Every oracle is a source equation's residual, a hand-derived
  closed form, or a number printed in a source. `tests/test_mop.c` reads only the fixtures under
  `tests/mop_fixtures/`, each of which carries its citation. `make test` treats `SKIPPED` as a failure
  unless the model is listed in an explicit allowlist.
- **G2. Coverage of the Giannantoni code.** `engine.c`, `mop.c`, `relational.c` and `harmony.c` join
  the coverage build with a ≥ 90% line gate. An unparsable percentage fails, and the parse uses
  portable `awk` instead of `grep -P`.
- **G3. No stubs in planned API.** A script checks that every function declared in `include/mop.h`
  and `include/relational.h` is called from a test, and that no body returns an "unimplemented" code.
  `AGENTS.md` rule 3 is scoped accordingly.
- **G4. Errata protocol.** Any departure from a printed equation needs all three of: a probe in
  `docs/sources/probes/`, a row in §5, and a test asserting that the *printed* form fails the oracle.
  No tolerance is widened without a §5 row.
- **G5. Definitions cite sources.** An ADR may not define a Giannantoni term (ordinality, harmony,
  incipient) except by citing a source equation, or by citing a §5 erratum.
- **G6. Wording.** Remove B5's equations of IDC with `expm`. Scope B6 to the kernel; MOP code never
  clamps, and a test checks that negative and complex coordinates survive.
- **G7. TDD evidence.** Each code PR body maps tests to equations and quotes the red run from before
  the implementation. Each mutation in §8 is run by hand and its failure quoted.

---

## 2. What is wrong today (unchanged from r1, corrected where noted)

| # | Claim | What the code does | Verdict |
|---|---|---|---|
| E1 | README: "a second engine implementing … IDC and MOP" | See E2–E8 | Overclaim |
| E2 | Kernel `"method": "incipient"` is an IDC solver (`src/gssk.c:6-13`) | Padé (3,3) `expm(A·dt)`; `limit` linearised | Classical numerics under an IDC label |
| E3 | MOP engine trajectories `_Q` | `gia_matrix_exp` (`src/engine.c:382`) | Classical |
| E4 | `_idc/_tdc/_drift` columns | Per-node φ from invented rules (`phi_for_node`, `src/engine.c:1381`), decoupled from `_Q` | Identity right, inputs invented |
| E4b | (hazard) drift from `φ = ln Q_i` | Would report nonzero drift for constant-coefficient networks, contradicting `[06 §4 (i)]` | Must not be built; T-IDC-11 guards it |
| E5 | Harmony Relationships (`src/engine.c:172-245`) | Constructed from roots of unity, then checked against the construction | Tautology (B11) |
| E6 | Ordinality (`src/engine.c:2673`) | Fraction of components on a cycle | Proxy; not `[22 Eq 11.1]` |
| E7 | Generative step (`gia_generate`) | ADR 0015 graph heuristic | Not from `[02 Eq 5.3]` or `[22 Eq 2]` |
| E8 | First and Second Fundamental Equations | Absent | Missing |

**Correction kept from r1.** `gia_idc_amplitude`'s `(φ')ⁿ` is Giannantoni's own definition, not a
shortcut. It is stated identically in `[02 Eq 14.8.2, 14.9.4–14.9.5]`, `[09 Eq 12, 25]`,
`[10 Table 1]`, `[22 Eq 5.2]` and `[23 Eq 3.2.5]`.

---

## 3. How the theory developed across the sources

Some of these are stable across twenty years and some shift. The plan relies on the stable ones and
resolves the shifting ones in §4:

| Element | Introduced | Later form | Status |
|---|---|---|---|
| Incipient derivative `(d̃/dt)ⁿ f = (f'/f)ⁿ f` | `[02 Eq 14.9.1–14.9.5]` | Same in `[06]`, `[09]`, `[10]`, `[22]`, `[23 Eq 5.5.2]` | **Stable** |
| Drift = Faà di Bruno − `(φ')ⁿ`, zero iff φ affine | `[02 p.175]` | `[06 §4]`, `[09 Eq 11–13]`, `[10 Eq 14–16]` | Stable |
| Binary function (co-production, order ½) | `[02 App. 7, Eq 14.7.1–14.7.7]` | `[06 Eq 3.10–3.15]` | Stable; ambiguity resolved in R9 |
| Emergy rules as ordinalities: co-production ½, interaction 2, feedback 2/2 | `[02 App. 10]`, `[06b Eq 5–15]` | `[10 §Incipient Derivative]`, `[22 Eq 6–8]` | Stable |
| Circle product | `[02 App. 11, Eq 14.11.5]` | `[06b Eq 2–3]` | Two notations, one rule (R12) |
| Maximum Em-Power Principle, formal | `[02 Eq 5.1–5.9]` | `[22 Eq 2, 4]` | Stable |
| MOP `(d̃/dt)^{(m/n)} {r} [→ {0}` | `[10]` | Two Fundamental Equations `[22 Eq 11–12]`, `[23 Eq 4.1–4.2]` | Second Equation's explicit form only in `[23 Eq 6.1–6.3]` |
| Order of m and n in `(m n)` | `[10]`: exponents ½, 2, 2/2 | `[22]`: "m Interactions, n Co-productions"; `[23]`: reversed | **Shifted**; resolved by the mathematics (R5) |
| Harmony Relationships, roots of unity, EQS | `[22 §12, App. 1–2]` | `[23 §5.6, 6–8, App. A1–A2]` | `[22]` defers definitions to its ref. [26], which is `[23]` App. A, so `[23]` is the definition of record |

---

## 4. Decisions of record

Each entry gives the question, the evidence, and the decision. None is left open.

**R1 — First Fundamental Equation per couple.** `[23 Eq 5.4.2]`: `(d̃/dt)^k α_ij = β_ij`, with the
incipient derivative `(α̇/α)^k α` `[23 Eq 5.5.2]`, `[02 Eq 14.9.5]`.

- *Evidence:* the printed solution is wrong for k ≠ 1 (§5 X1; probe `eq_5_5_8_residual.py`).
- *Decision:* use the solution derived from `[23 Eq 5.5.6]`:
  `α(t) = { (1/k) ∫₀ᵗ β^{1/k} }^k`. For `β = (a+bt)^p` with b ≠ 0 this is
  `α = [((a+bt)^{(p+k)/k} − a^{(p+k)/k}) / (b(p+k))]^k`, and with b = 0 it is `α = (β^{1/k} t/k)^k`.
- The lower limit 0 of `[23 Eq 5.5.5]` fixes `α(0) = 0`; this is not configurable.
- Fractional k uses the principal branch of `β^{1/k}` in ℂ.
- For a general β given as samples, the integral is evaluated by adaptive Gauss–Kronrod quadrature.

**R2 — Relational product and exponential.**

- *Evidence:* `[23 Eq 5.1.3–5.1.5]` (page image p. 3183) gives
  `i∘i=+1, i∘j=j, i∘k=k, j∘i=j, j∘j=−1, j∘k=k, k∘i=k, k∘j=k, k∘k=−1`. Read literally and bilinearly,
  this is commutative and **non-associative**: `(j∘j)∘k = −k ≠ j∘(j∘k) = k`.
- *Test of the reading:* take the EQS root in De Moivre form, `r_l = cos(√2ψ_l) + (j+k)·sin(√2ψ_l)/√2`
  (that is, `B_l, C_l, D_l` of `[23 Eq 7.4]`). Then `r_l ∘ {Σ₀, Φ₀, Θ₀}` under the literal table
  reproduces the three brackets of `[23 Eq 7.1.1, 7.2, 7.3]` to **8.9 × 10⁻¹⁶** over 1000 random
  draws. The associative alternative (`j∘k + k∘j = 0`) does **not**: it misses Eq 7.3's
  `+C_l(Φ₀+Θ₀)` term (probe `eqs_relational_product.py`).
- *Decision:* adopt the literal table as printed. Products are evaluated strictly left to right in
  the order written. Exponentials use the generalised De Moivre form `[23 Eq 5.1.2]`:
  `Exp{a·i ⊕ b·j ⊕ c·k} := eᵃ [cos ρ + (b·j + c·k) sin ρ / ρ]`, with `ρ = √(b²+c²)`; this is never a
  power series.

**R3 — Ordinal roots of unity.**

- *Evidence:* `[23 Eq A2.5–A2.6]` define `(ᴺ⁻¹√1)_l = Exp{(α i ⊕ β j ⊕ γ k)/(N−1)}` with periods 4π
  (i) and 2π (j, k). Unity comes from **angle periodicity**.
- Under repeated *table* products, the EQS root is not a root of unity. For N = 4, l = 1, `r∘r∘r`
  gives `(0.5407, 0, −0.6657)`, not `(1, 0, 0)`; it reaches unity only where sin = 0 (probe
  `ordinal_root_power.py`).
- *Decision:* the root's power is defined by angle multiplication (De Moivre), not by ∘.
  `relational.c` exposes both, named distinctly (`rel_root_pow`, `rel_mul`). Tests pin both
  behaviours, so neither can be mistaken for the other.

**R4 — Second Fundamental Equation.** No source writes `[23 Eq 4.2]` as an explicit differential
equation. `[23 Eq 6.1–6.3]` gives only its solution, which `[23 §6]` calls the solution of "a typical
Riccati's Equation of Ordinal Nature".

- *Evidence:*
  - In `[23 Eq 6.3]` the only time dependence is `ln(c₁ ⊕ {c₂, t})`.
  - `{c₂, t}` is a duet, and its cardinal reduction is `c₂·t` (R12).
  - `d/dt ln(c₁ + c₂t) = c₂/(c₁+c₂t)` is the general solution of the Riccati equation `u' + u² = 0`.
  - `[23 Eq 6.2]` gives `B(t) = [[⊕A, ΘA], [ΘA, ⊕A]]`, the duet-binary (2/2) specular form of
    `[22 Eq 8]`.
- *Decision:* implement the solution as printed:
  - `A(t) = (α₁₂(0) + λ₁₂(0)) ∘ (ᴺ⁻¹√1)^{↑N/N} + ln(c₁ + c₂t)`;
  - `B(t) = [[A, −A], [−A, A]]`;
  - `{r} = e^{B(t)} ∘ (roots 13 … 1N)` per `[23 Eq 6.1]`;
  - λ null per `[23 §8 ii]`.
- The Riccati `u' + u² = 0` that `Ȧ` satisfies is the oracle. It is labelled in code and docs as
  *reconstructed from the printed solution*.

**R5 — Ordinality `{k, (m n)}`.**

- *Evidence:* `[10]` assigns exponents ½, 2 and 2/2 to co-production, interaction and feedback.
  `[22 Eq 6–8]` gives `(d̃/dt)^{1/2}` (binary), `(d̃/dt)^2` (duet) and `(d̃/dt)^{2/2}` (duet-binary).
  So in `q = m/n` the numerator counts interactions (the power) and the denominator counts
  co-productions (the number of root branches). `[22]`'s wording agrees with this; `[23]`'s swaps it
  (§5 X6).
- Maximum is `{2 2} ↑ {N N}` `[22 Eq 11.1]`, and "all the various couples … become of Ordinality
  {2/2}" `[22 §12.1]`.
- *Decision:*
  - **Classify each couple of components:**
    - 2/2 if each reaches the other along quantity-carrying legs (the ADR 0014 walk);
    - otherwise 2 if both feed one interaction module;
    - otherwise ½ if both are products of one replicating process;
    - otherwise unrelated.
  - **k** is the number of components. Boundary nodes (source, sink, constant) are habitat, not
    components: `[10]`'s "including that of the surrounding habitat", and `[23 Eq 6.3]` treats
    habitat as conditions, not matrix entries.
  - **Maximum ordinality** holds iff every couple is 2/2, which is equivalent to the component graph
    being **strongly connected**.
  - Today's measure is renamed **closure** and kept for reporting only. ADR 0014's leg rules stay;
    ADR 0014/0015's definitions of ordinality are superseded.

**R6 — The generative step.**

- *Evidence:* the Maximum Em-Power Principle `[02 Eq 5.3]`, `[22 Eq 2]`:
  `∫Γφ* dV = d/dt ∫ em* dV → Max`. The MOP is "nothing but the re-proposition of the Maximum Em-Power
  Principle" `[10 Abstract]`, `[22 §9]`.
- *Decision:* while the system is below maximum ordinality (R5), add **one** quantity-carrying
  pathway, from a component in a sink strongly-connected component of the condensation to one in a
  source component.
  - Among all such candidates, choose the one that maximises total empower at `t_end`, as computed by
    the emergy pass.
  - Break ties by lexicographic `(from, to)` id.
  - Repeat until strongly connected. This terminates in at most `#sources + #sinks` steps, because
    each addition removes a source or a sink of the condensation.
  - New pathways are `linear` with the seed's mean weight. ADR 0015's emergent component `E` is
    retired, because `[22 §12.1]` has couples *become* 2/2, which is a change of relation, not a new
    component.

**R7 — Harmony: emerge, transported, or imposed?** All three verdicts are determined by the sources.

- **First Equation: transported.** It is decoupled per couple and homogeneous of degree 1, so
  `α_1j = ω α_12 ⟺ β_1j = ω β_12` (r1 §2.4, Proposition). Harmony in the solution is exactly the
  harmony supplied in β.
- **Second Equation: imposed.** `[23 Eq 6.1]` writes the solution with the roots
  `(ᴺ⁻¹√1)_{13…1N}` explicitly. Its only free inputs are `α₁₂(0)`, `λ₁₂(0)`, `c₁` and `c₂`. Couples
  `1j` have no input through which they could differ, so harmony is in the formula.
- **EQS: imposed, by its own statement.** "originate from the assumption that the Harmony
  Relationships…" `[23 §8 ii]`.
- **The "Diffusive Generativity" from which `[22 §12]` and `[23 §5.6, App. A1]` say harmony emerges is
  never written as an equation.** `[23 App. A1]` calls it "not a necessary consequence".
- *Decision:* the engine reports one of `transported`, `imposed` or `absent` for each construction,
  computed by the detector in §8 and never hard-coded. The documentation states plainly that no
  equation in the available sources produces harmony. §6 lists emergence as not reproducible.

**R8 — Network to Relational Space.**

- *Evidence:* `[22 §12.1]` says an emergy diagram "can be transformed into the Relational Space …
  as shown in [26]", and [26] = `[23 App. A]` does not do so.
- *Decision:*
  - Define `e^{α_ij(t)}` := empower carried from component i to component j by the emergy pass. This
    rests on MOP = MEmPP re-proposed (R6).
  - Derive `β_ij` from it by `[23 Eq 5.5.2]` with k = 1.
  - Couples with no pathway are `unrelated` and have no entry.
- *Consequence, stated and tested:* the values are real, so for N ≥ 4 the ratios cannot equal the
  non-real roots of unity, and the detector necessarily returns `absent` (T-HARM-2b). The sources give
  no rule for assigning a phase to a network couple (§6).

**R9 — Binary function `[06 Eq 3.10–3.15]`, `[02 Eq 14.7.1–14.7.7]`.**

- *Evidence:* `f' + A f^{(½)} + B f = 0`, terms `e^{s²t}`, branch characteristics `s² ± A s + B = 0`.
  The − branch's roots are the negatives of the + branch's, so they give the same exponents `s²` and
  the same half-derivative values.
- *Decision:* solve `u² + A u + B = 0` once, giving `u₁, u₂`. Both branch functions are
  `f_σ = c_σ1 e^{u₁²t} + c_σ2 e^{u₂²t}`, with half-derivative `Σ c_σi u_i e^{u_i² t}`. The constants
  come from the four initial conditions `[02 Eq 14.7.6–14.7.7]` through
  `[[1, 1], [u₁, u₂]] c_σ = (f_σ0, f_σ0^{(½)})`. The ± labels are bookkeeping (§5 X7).

**R10 — The incipient Taylor series `[09 Eq 10]`, `[10 Eq 13]`.**

- *Definition:* `f*(t₀+Δ) = f(t₀) Σ_{k=0}^{n} (aΔ)^k/k!`, with `a = f'(t₀)/f(t₀)`.
- *Evidence:* with **n = 2** this reproduces `[09]`'s published values: 16.4 (Eq 16), 3.0125 (Eq 17,
  printed 3.01), 178.0 (Eq 21), 172.067 (Eq 19 at τ₀ = 1.8, printed 172.06), and 15.6 / 17.21 (Eq 22,
  printed "15.0–17.0"). Probe `zoli_2009_reproduction.py`.
- *Decision:* n is a parameter, and the tests use n = 2. Two printed values do not reproduce (§5 X4,
  X5).

**R11 — Riccati `[06 Eq 3.16–3.24]`.**

- Printed Eq 3.17 is inverted (§5 X2). Printed Eq 3.22 cannot be derived from Eq 3.19 by any reading
  (§5 X3), and its source `[06 ref 8]` (Giannantoni 2004, *Differential Bases of Emergy Algebra*) is
  not in `docs/`.
- *Decision:* implement the route the paper also gives. Apply the corrected substitution
  `f = y'/(R y)` to reach `[06 Eq 3.18]`, then solve that with the incipient LDE solver (R13). The
  direct duet form of Eq 3.22 is out of scope (§6).

**R12 — Circle product.**

- *Evidence:*
  - `[02 Eq 14.11.4–14.11.5]`: `l∘l ≠ l²`, and `l∘l = [l, l]`, "a du-et of real numbers".
  - `[06b Eq 2]`: `(a₁; a₂) ∘ [b₁, b₂] = [(a₁b₁; a₂b₁), (a₁b₂; a₂b₂)]`.
  - `[06b Eq 3]` gives the same with the entries left as symbols.
- *Decision:* ∘ is structural. It keeps every pair of factors in the outer arrangement, and the
  **cardinal reduction** maps each pair to its product. `[06b Eq 2]` is the reduced form and
  `[02 Eq 14.11.5]` the unreduced one; `{c₂, t}` in `[23 Eq 6.3]` reduces to `c₂·t`.

**R13 — Incipient LDE of order 2 `[06 Eq 3.3–3.7]`, `[09 Eq 3, 7]`, `[10 Eq 8.1, 10.1]`.**

- Characteristic `α̃² + a₁(t) α̃ + a₀(t) = 0`, solved pointwise.
- `f = Σ c_i e^{∫₀ᵗ α̃_i}`, with `c_i` from `f(0)` and `f̃'(0) = Σ c_i α̃_i(0)`.
- At a double root, `[06 Eq 3.7]` gives the second solution.
- The worked example `[10 App. Eq 28–34]` (a zero root) is a fixture.

**R14 — Emergy algebra fixtures.**

- Co-production source term `Φ(u) = (n−1)·Em(u)` `[02 Eq 3.8]`.
- Interaction in steady state `Em(y) = Em(u₁) + Em(u₂)`, i.e. `Φ(u₁,u₂) = 0` `[02 Eq 3.9, 3.12, 3.15]`.
- Split `Em(y₁) = x·Em(u)`, `Em(y₂) = (1−x)·Em(u)` `[02 Eq 3.16–3.17]`.
- Rule 4 (no double counting of reunited co-products or of feedback) `[02 p. 23]`.
- `[02 Fig. 3.4]` (Brown 1993) labels emergy values only; its internal values are not derivable from
  the four rules without Brown's energy flows, which are not in `docs/`. Only its global balance
  `[02 Eq 3.23–3.24]` (`½·Em(Z) + ½·4·Em(Y) = Em(S) + Em(F) = 30,000`) is used, as an arithmetic
  fixture of the global-balance function. The figure is kept as a transcription in
  `docs/sources/local/`.

**R15 — Kernel `"method": "incipient"`.**

- *Decision:* add `"expm"` as the documented name. Keep `"incipient"` accepted as a deprecated alias
  that emits a one-line notice saying it is the matrix exponential.
- The kernel offers no incipient method: per ADR 0011 the kernel is the classical (TDC) engine.
- This is not a schema break.

---

## 5. Errata register

Every row has a probe in `docs/sources/probes/` and a test asserting the printed form fails (G4).

| # | Printed | Problem | Evidence | Used instead |
|---|---|---|---|---|
| X1 | `[23 Eq 5.5.7–5.5.8]` | Prefactor `1/k` outside the k-th power; 5.5.8 drops the lower limit and misplaces `1/b` | Residual of `[23 Eq 5.5.2]` at a=1, b=0.25, p=1, k=2, t=1: printed 5.5.7 → 1.25, printed 5.5.8 → −0.625, derived → 0 (`eq_5_5_8_residual.py`) | R1 |
| X2 | `[06 Eq 3.17]` `y = f'/(fR)` | Inverted substitution | With Q=R=1, P=2: the standard `f = y'/(Ry)` gives Riccati residual 1.3 × 10⁻⁸; the printed form gives Eq 3.18 residual 929 (`riccati_substitution.py`) | R11 |
| X3 | `[06 Eq 3.22]` | Not derivable from Eq 3.19; `5 + Q(t)` adds a pure number to a rate | Dimensional inspection of the page image (`local/pages/giannantoni-2006-jcam-p331.png`) | Out of scope (§6) |
| X4 | `[09 Eq 19]` τ₀ = 2 → 154.3 cm | Same formula gives 156.0 (τ₀ = 1.8 reproduces) | `zoli_2009_reproduction.py` | Test pins 156.0 |
| X5 | `[09 §2]` minimum scenario "net increase 1.91 °C … 73 %" | Subtracts 1.1, where the maximum scenario subtracts 0.4; consistent net is 2.61 | Same probe | Test pins 3.0125 and 2.6125 |
| X6 | `[23 Eq 4.1.1]` text "m Co-productions and n Interactions" | Reverses `[22]` and the exponent semantics | `[22 Eq 6–8]`, `[10]` | R5 |
| X7 | `[02 Eq 14.7.5]`, `[06 Eq 3.13]` distinct exponents per branch | The two branch characteristics have identical exponent sets | R9 derivation | R9 |
| X8 | `[02 Eq 14.10.2]` "triplet" of solutions | The incipient characteristic of Eq 14.10.1 is `4u² + Au + B = 0`, which has two roots | Substitution `F = e^{ut}` | Two roots; T-IDC-8 |
| X9 | `[23 Eq 5.6.5]` at j = 1 reads `α₁₂ = ω₁ α₁₂` | Literal index collides with the reference couple | `[23 App. A1]` "updates … the same reference couple" | Detector uses ratios (§8) |
| X10 | `[23 A2.5]` "roots of unity" vs `[23 Eq 5.1.3–5.1.5]` | Not roots of unity under the table product | `ordinal_root_power.py` | R3 |

---

## 6. Not reproducible from the available sources (scope, stated)

These are excluded **with a reason**. None is an open decision. The engine refuses them by name, and
the docs say so.

| Item | Why | What the engine does |
|---|---|---|
| Emergence of Harmony from "Diffusive Generativity" | Never written as an equation (R7) | Reports `transported`/`imposed`/`absent` |
| Solar-system distances `[22 Table 3]`, Mercury 42.45″/cy `[06]`, `[09]`, `[10]`, three-body and N-body closed forms `[06b]`, `[10p]` | Parameters (ε, ψ, A, initial data) are not given; the methods are in `[10p]`'s refs (Giannantoni 2007a/2008b, Italian books, "Solving Kernel" theorem) | Not attempted; listed in docs |
| Riccati duet `[06 Eq 3.22]` and Abel n-et `[06 Eq 3.25–3.27]` | X3; `[06 ref 8]` not available | Riccati via R11; Abel refused |
| Phase of a network couple | No source assigns one (R8) | Real empower mapping; `absent` for N ≥ 4 |
| Relational-space α for k > 1 | Division and powers in the non-associative algebra are undefined by the sources | First Equation over ℂ for any k; relational-valued α only for k = 1 (componentwise integration); k > 1 relational refused with a named error and tested |
| `[02 Eq 3.21]` re-normalisation factors (α*, β*, γ*) in general | Chosen per system in the text ("equal to 4 because of the global structure"), with no stated rule | Global balance only as the `[02 Eq 3.23–3.25]` arithmetic fixture |

---

## 7. Workstreams and tasks

Order inside every task: **red** (write the §8 tests, run, quote the failure) → **green** →
**mutations** (§8 column, each must fail) → `make test && make test-advanced && make test-schema &&
make test-mop && make ci-local`. Slugs are proposed crux task names.

### W0 — Remediate the claims (docs only)

| Slug | Deliverable | Acceptance |
|---|---|---|
| `mop-claims-remediation` | README intro, `include/engine.h` header, `src/gssk.c` header comment, `docs/giannantoni_assessment.md` status table (§2, E4b, the correction), TODO 10.3 annotated "constructed", CHANGELOG | No text says the repo implements IDC or MOP until W4–W8 land. Each feature is labelled `implemented`, `classical`, `constructed`, `proxy` or `absent` |

### W1 — Guardrails (§1)

| Slug | Deliverable | Acceptance |
|---|---|---|
| `guard-no-skip` | `make test`: `SKIPPED` fails unless the model is in `tests/skip_allowlist.txt` | Removing one golden file turns `make test` red |
| `guard-coverage-giannantoni` | Coverage build covers `engine.c` + new units; ≥ 90% gate; unparsable fails; portable parse | Gate fails at 89% (mutation: delete a test) and on garbage input |
| `guard-api-called` | `scripts/check_api_called.sh`: every `mop.h`/`relational.h` symbol is referenced in `tests/test_mop.c`; no `UNIMPLEMENTED` returns | Adding an untested declaration fails CI |
| `guard-agents-wording` | `AGENTS.md`: Overview, rule 3 scoped (G3), Fail-Safe scoped (G6), a "Giannantoni work" section pointing at §1 G1–G7; TODO Non-Goals/Vision fixed (B5) | Review |
| `adr-giannantoni-test-protocol` (ADR 0018) | Records G1–G7 and the errata protocol | Merged before any W2+ code PR |

### W2 — IDC (single-variable)

| Slug | Deliverable | Tests |
|---|---|---|
| `idc-general-f` | `gia_idc_of(f, f', n)` for general f per `[02 Eq 14.9.5]` | T-IDC-1, 2, 3 |
| `idc-lde2` | `gia_idc_solve_lde2(a1, a0, f0, f1)` per R13, incl. double root | T-IDC-4, 5, 10 |
| `idc-binary` | `gia_idc_solve_binary(A, B, F0[2], F0half[2])` per R9 | T-IDC-6, 10 |
| `idc-riccati` | `gia_idc_solve_riccati(Q, R, P, f0)` per R11 | T-IDC-7 |
| `idc-taylor` | `gia_idc_taylor(f0, f1, dt, n)` per R10 | T-IDC-9 |
| `idc-nonlinear-14-10` | Characteristic of `[02 Eq 14.10.1]` | T-IDC-8 |
| `idc-drift-coupled` | Modal drift on network trajectories; `_idc/_tdc/_drift` columns re-specified; `phi_for_node` heuristics removed from reported results | T-IDC-11 |

### W3 — Emergy algebra in IDC form

| Slug | Deliverable | Tests |
|---|---|---|
| `emergy-source-terms` | `gia_emergy_source_term()` per R14 on the existing pass | T-EM-1…4 |
| `emergy-ordinal-forms` | Binary/duet/duet-binary constructors `[22 Eq 6–8]`; circle product with cardinal reduction (R12) | T-EM-5, 6 |

### W4 — MOP First Fundamental Equation

| Slug | Deliverable | Tests |
|---|---|---|
| `adr-mop-equations` (ADR 0019) | R1, R4, R7 (verdict semantics) and X1 recorded with probes | Merged first |
| `mop-first-equation` | `src/mop.c`, `include/mop.h` (no `gssk.h`, no `harmony.o`): `gia_mop_solve_couple`, `gia_mop_solve` (Matrioska, internal representation `[23 Eq 5.6.1]`) | T-MOP-1…4 |

### W5 — Relational algebra and EQS

| Slug | Deliverable | Tests |
|---|---|---|
| `adr-relational-algebra` (ADR 0020) | R2, R3, X10 with probes | Merged first |
| `mop-relational-algebra` | `src/relational.c`: `rel_mul` (table, left to right), `rel_exp` (De Moivre), `rel_root`, `rel_root_pow` | T-MOP-5, 6 |
| `mop-eqs` | `[23 Eq 7.1–7.5]` built on `rel_mul` ∘ `rel_root`; tagged harmony-assuming | T-MOP-7 |

### W6 — Second Fundamental Equation

| Slug | Deliverable | Tests |
|---|---|---|
| `mop-second-equation` | R4: `A(t)`, `B(t)`, `{r}` per `[23 Eq 6.1–6.3]` | T-MOP-8 |

### W7 — Ordinality and the generative step

| Slug | Deliverable | Tests |
|---|---|---|
| `adr-ordinality-mop` (ADR 0021) | R5, R6; supersedes ADR 0014's and 0015's ordinality definitions and ADR 0015's emergent component; lists every `examples/giannantoni/` verdict before and after | Merged first |
| `mop-ordinality` | `gia_ordinality_record → {k, n_22, n_2, n_half, n_none}`, `gia_closure` (renamed proxy), `gia_at_maximum_ordinality` = strong connectivity | T-MOP-9 |
| `mop-generative-empower` | `gia_generate` per R6 | T-MOP-10 |

### W8 — Harmony verdicts

| Slug | Deliverable | Tests |
|---|---|---|
| `mop-harmony-detector` | `gia_harmony_verdict()` in `mop.c`. Today's constructor moves to `src/harmony.c` as `gia_harmony_assume_*`; old tests renamed to say they test the constructor | T-HARM-0…5 |
| `mop-network-beta` | R8 mapping; verdict per seed in a committed results table | T-HARM-2b on every seed |
| `kernel-method-label` (ADR 0022 + code) | R15 | Kernel test: `"incipient"` warns and equals `"expm"` byte for byte |

---

## 8. Test suite

Binary `bin/test_mop` (`tests/test_mop.c`) runs via `make test-mop`, with its `.PHONY` declared beside
the rule and added to `CI_TESTS` and `deploy.yml` in the PR that creates it. A second binary,
`bin/test_mop_emergence`, links `mop.o` and `relational.o` **without** `harmony.o`.

Every test function opens with `/* Source: [key Eq n] */`.

**Tolerances** (no widening without a §5 row):

| Kind | Tolerance |
|---|---|
| Closed forms | 1e-12 relative |
| Quadrature | 1e-9 |
| Central-difference residuals (`h = 1e-5·max(1, t)`) | 1e-6 relative |
| Published numbers | the printed precision |

Every "must fail" assertion requires at least 100× tolerance.

### 8.1 IDC

| ID | Source | Oracle | Pass | Mutation that must fail |
|---|---|---|---|---|
| T-IDC-1 | `[02 Eq 14.8.2]`, `[23 Eq 3.2.5]` | `(φ')ⁿ e^φ` by hand for polynomial φ, n = 0…8 | 1e-12 | Use the Bell polynomial instead |
| T-IDC-2 | `[02 Eq 14.9.5]` | `f = 1+t²`: `(2t/(1+t²))ⁿ (1+t²)` | 1e-12 | Drop the `/f` |
| T-IDC-3 | `[02 p.175]`, `[09 Eq 11–13]`, `[10 Eq 16]` | Drift `B_n − (φ')ⁿ`: `ψ₂ = a`, `ψ₃ = 3a²t` for `φ = a t²/2`; zero iff φ affine | Exact | Reverse the sign |
| T-IDC-4 | R13, `[06 Eq 3.3–3.6]` | (i) `a₁=0, a₀=−t²` ⇒ `α̃=±t`; incipient residual 0; ICs met; traditional residual = `c₁e^{t²/2} − c₂e^{−t²/2}`. (ii) Double root `a₁=−2t, a₀=t²` ⇒ `[06 Eq 3.7]` | 1e-12 / 1e-6 | Use the traditional characteristic `[10 Eq 10]` |
| T-IDC-5 | `[10 App. Eq 28–34]` | `u² + ψ_f u = 0`, `ψ_f = 1+t`: `g = C₁ + C₂e^{−∫ψ_f}`; `F = f·g̃'` is **not** constant although the incipient condition holds (the paper's "24.1 ⇏ 24") | F(1) ≠ F(0) by > 100× tol | Force F constant |
| T-IDC-6 | R9, `[02 Eq 14.7.1–14.7.7]` | Residual `f' + A f^{(½)} + B f` per branch; four ICs reproduced | 1e-12 | Swap `u₁ ↔ u₂` in one half-derivative |
| T-IDC-7 | R11, X2 | Constant Q, R, P: Riccati residual of `f = y'/(Ry)`; printed Eq 3.17 gives Eq 3.18 residual > 100× tol | 1e-6 | Printed substitution |
| T-IDC-8 | `[02 Eq 14.10.1]`, X8 | `4u²+Au+B=0`: each root's `F = e^{ut}` has zero incipient (and, φ being affine, zero traditional) residual; exactly two roots | 1e-12 | Coefficient 1 instead of 4 |
| T-IDC-9 | R10, X4, X5 | 16.4; 3.0125; 178.0; 172.0667; 15.6; 17.2067; Eq 19 τ₀=2 → **156.0** | printed precision | n = 3 |
| T-IDC-10 | `[06b]` "linearly dependent on initial conditions" | `solve(ic₁ + ic₂) = solve(ic₁) + solve(ic₂)` for T-IDC-4, 6 | 1e-12 | Clamp a coefficient |
| T-IDC-11 | `[06 §4 (i)]`, E4b | Constant flow matrix: modal drift ≡ 0, n = 1…4, every seed. Interaction module: drift nonzero and stable as dt → 0 | 0 exactly / 1e-6 | `φ = ln Q_i` (must give nonzero on a two-mode network) |

### 8.2 Emergy algebra

| ID | Source | Oracle | Pass | Mutation |
|---|---|---|---|---|
| T-EM-1 | `[02 Eq 3.8]` | Co-production excess `= (n−1)·Em(u)`, n = 2, 3, 4 | 1e-12 | Partition instead of replicate |
| T-EM-2 | `[02 Eq 3.16–3.17]` | Split shares x, 1−x; excess 0 | 1e-12 | Replicate |
| T-EM-3 | `[02 Eq 3.9, 3.12, 3.15]` | Drawn interaction: `Em(y) = Em(u₁) + Em(u₂)`; Φ = 0 | 1e-12 | Product `k·Em₁·Em₂` |
| T-EM-4 | `[02 p.23 rule 4]` | Reunited co-products take the max; feedback not re-injected | 1e-12 | Sum |
| T-EM-5 | `[22 Eq 6–8]`, `[06b Eq 6, 10]` | Binary: 2 branches each `Em(u)`; duet pair; duet-binary `[[a₁,a₂],[a₂,a₁]]` | Exact | Transpose |
| T-EM-6 | R12 | `[06b Eq 2]` reduced form; `l∘l` reduces to l² but is stored as `[l, l]` | Exact | Store l² |
| T-EM-7 | `[02 Eq 3.23–3.26]` | `gia_emergy_global_balance` on the Fig. 3.4 totals (S = 10,000, F = 20,000, Z = 30,000, Y = 7,500). Case A, `1·S + 1·F = ½·Z + ½·4·Y`, balances. Case B, `S + F + Φ_D + Φ_E = Z + 6·Y` with `Φ_E = ½Φ_D`, solves `Φ_D = 30,000` and `Φ_E = 15,000` | Exact | Factor 4 in case B |

### 8.3 MOP

| ID | Source | Oracle | Pass | Mutation |
|---|---|---|---|---|
| T-MOP-1 | R1, `[23 Eq 5.4.2, 5.5.2]` | Residual on a grid `k ∈ {1, 2, 3}` × `p ∈ {0, ½, 1, 2}` × `b ∈ {0, 0.25, 1}` with real and complex a, b; plus `k = ½` with real `a, b > 0` (principal branch, R1); closed form vs quadrature | 1e-6 / 1e-9 | Printed 5.5.8 |
| T-MOP-2 | X1 | Printed 5.5.7 / 5.5.8 fail at (1, 0.25, 1, 2): 1.25, −0.625 at t = 1 | > 100× tol | — |
| T-MOP-3 | R1 | `α(0) = 0`; b = 0 branch | 1e-12 | `α(0) = a^{…}` |
| T-MOP-4 | `[23 Eq 5.6.1]` | Zero diagonal; N(N−1) independent couples | Exact | — |
| T-MOP-5 | R2, `[23 Eq 5.1.3–5.1.5]` | All nine products; `(j∘j)∘k = −k`, `j∘(j∘k) = k` | Exact | Anticommuting j, k |
| T-MOP-6 | R3, `[23 Eq A2.5–A2.6]` | With `√2·ψ_l = 2πl/(N−1)`: `rel_root_pow(r, N−1) = 1` for N = 3…7; `rel_mul(rel_mul(r, r), r)` at N = 4, l = 1 = (0.540721, 0, −0.665721), as `ordinal_root_power.py` prints | 1e-12 / 1e-6 | Power by `rel_mul` |
| T-MOP-7 | `[23 Eq 7.1–7.5]` | Hand-evaluated N = 4 with given Σ₀, Φ₀, Θ₀, ε, ψ; equals `rel_mul(rel_root, X)` brackets | 1e-12 | Anticommuting product (fails Eq 7.3) |
| T-MOP-8 | R4, `[23 Eq 6.1–6.3]` | `Ȧ` satisfies `u' + u² = 0`; `B = [[A,−A],[−A,A]]`; `{c₂,t} → c₂t` | 1e-6 / exact | `ln(c₁ + c₂ + t)` |
| T-MOP-9 | R5 | Six hand graphs: classification counts exact; **two disjoint 2-cycles: closure = 1, not maximum**; boundary nodes excluded; reorder-invariant | Exact | Cycle-coverage definition |
| T-MOP-10 | R6, `[02 Eq 5.3]` | Below max ⇒ reaches max in ≤ #src + #sink steps; each choice equals brute-force argmax of total empower on graphs ≤ 5 components; fixed point at max; reorder-invariant | Exact | First candidate instead of argmax |

### 8.4 Harmony verdicts (`bin/test_mop_emergence`; no `harmony.o`)

The detector measures, in ℂ, `R_H = max_{j=2..N−1} |α_{1,j+1}/α_{12} − e^{2πi(j−1)/(N−1)}|`
(ratios, X9). The verdict is:

- `imposed` if `R_H` is 0 for every admissible perturbation of the inputs;
- `transported` if `R_H(α)` equals `R_H(β)`;
- `absent` otherwise.

| ID | Construction | Expected | Mutation |
|---|---|---|---|
| T-HARM-0 | Link without `harmony.o`; `nm mop.o` has no `gia_harmony_*` or `gia_ordinal_root` | Links | Call the constructor from `mop.c` (link fails) |
| T-HARM-1 | β harmonic, N = 4, 5, 7 | `R_H < 1e-9` | — |
| T-HARM-2 | Random complex β (fixed seed) | `R_H > 1e-2`, `absent` | Widen tolerance to 1 |
| T-HARM-2b | Any real Matrioska, N ≥ 4 | `R_H ≥ |sin(2π/(N−1))|`, `absent` | — |
| T-HARM-3 | First Equation, harmonic β + δ = 1e-3 | `transported` (`R_H(α) = R_H(β)`) | Return the constructed Matrioska (T-HARM-0, 2 fail) |
| T-HARM-4 | Second Equation (R4), random α₁₂(0), c₁, c₂ | `imposed` | — |
| T-HARM-5 | EQS (R2) | `imposed` | — |

---

## 9. Sequencing and PRs

PRs that add a CHANGELOG entry land **one at a time** (`AGENTS.md` §Concurrent branches). ADR PRs
precede their code PRs.

1. **This PR** (`docs/mop-plan`): `PLAN.md`, `docs/sources/` (cc-by text and pages, probes, README),
   `.gitignore` rule for `docs/sources/local/`.
2. W0 `mop-claims-remediation`.
3. W1 ADR 0018, then `guard-no-skip`, `guard-coverage-giannantoni`, `guard-api-called`,
   `guard-agents-wording`.
4. ADR 0019 (MOP equations), ADR 0020 (relational algebra), ADR 0021 (ordinality), ADR 0022 (method
   label): written together (separate files), merged one at a time.
5. Code, in dependency order: W2 → W3 → W4 → W5 → W6 → W7 → W8. `mop-first-equation` owns the
   `test-mop` target and its `CI_TESTS` / `deploy.yml` lines.

```
W0 ─┐
W1 (ADR 0018 → guards) ─┬─► W2 idc-* ─► W3 emergy-* ─┐
                        ├─► ADR 0019 ─► W4 mop-first-equation ─┬─► W6 mop-second-equation ─┐
                        ├─► ADR 0020 ─► W5 relational, eqs ─────┘                            ├─► W8 harmony, network-beta
                        ├─► ADR 0021 ─► W7 ordinality ─► generative-empower ─────────────────┘
                        └─► ADR 0022 ─► kernel-method-label
```

## 10. Definition of done

- §2 rows E1–E8 all read `implemented`, `classical (labelled)` or `out of scope (§6)`, and the README
  matches.
- Every test in §8 passes in `make test-mop` under `make ci-local`, every mutation in §8 fails and its
  failure is quoted in the PR, and the G2 coverage gate passes.
- `docs/` states the three harmony verdicts of R7 and the §6 exclusions.
- Everything is merged to `main`. Crux tasks are closed only after confirming the code is on `main`.
