# Changelog

All notable changes to GSSK are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project adheres to [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

### Changed — breaking for Giannantoni-engine models

- **Module laws leave pathways** ([ADR 0012](adr/0012-where-a-law-lives.md) decision 5). A pathway is now `linear`, `reversible` or `constant`, Odum §III, and nothing else. `interaction`, `limit`, `ratio`, `subtract`, `gain`, `threshold` and `exchange` on a pathway — and the seed labels `generative_production` and `ordinal_feedback`, which meant `interaction` — are load errors that name the module to write instead; the loader never rewrites them, because a gate the author did not draw would vote on ordinality. A module type (`interaction`, `gain`, `switch`, `loop_limited`, `exchange`) without its `module` block is refused too: it used to load as a stock wearing the name.

  **The example seeds are re-modelled as Odum's work gate**, the decision ADR 0014 left to this migration. The sun is the energy input, biomass the autocatalytic control at `use_ratio` 0, and the consumer's feedback a control drawn at `use_ratio` 0.1 into a heat sink. That is a different system from the old seed, whose `interaction_1` was a stock: trajectories change, the gate is no longer a component, and `closed_loop.json` now sits at ordinality 0.8 because a sink is never closed (ADR 0015). Both demo verdicts hold — `input.json` generative, `closed_loop.json` functional.

  **The GSSK projection writes the gate GSSK draws implicitly.** A GSSK edge law becomes an `interaction` (with `divide` or `subtract` for `ratio` and `subtract`), `loop_limited` or `switch` module with named legs, controls at `use_ratio` 0, and the coverage report lists each such gate under *Added in translation*. An n-ary `interaction` edge is now carried as a multi-control gate instead of being refused. A GSSK processing node is dropped with a finding even when it has no params: it was carried as a storage wearing the name. Pathway-law tests were converted to module form; the exchange-pathway tests whose properties the transactor-module tests already check were removed.

- **A module's control input is a flow of energy, and its emergy reaches the product** ([ADR 0017](adr/0017-what-a-control-input-is.md)). A control leg used to be read and never drawn, so it carried no flow and the emergy pass never saw it: a work gate with a transformity-1000 control gave its product transformity 1, the same as a transformity-1 control, where Odum's Fig. 2.6(b) draws the product "moderate". A control leg now states a **required** `use_ratio`, draws `use_ratio · F` from its source, and dissipates it on a new `used` leg leaving the module, one per carrier, to a `sink` of the control's carrier. Its emergy goes into the product and the `used` leg carries none onward: the measured gate now gives transformity 11 (`1 + 0.01 · 1000`), between its two inputs. `use_ratio: 0` keeps the old read-only behaviour, which is right for a price, and the run report names every such leg.

  Load refuses a control with no `use_ratio`, or a negative or non-numeric one; `use_ratio` anywhere but a control leg; a drawn control with no `used` leg for its carrier; and a `used` leg that does not end at a sink, duplicates another's carrier, or has nothing to carry. A control drawn from a source or a constant now makes the system open, and the cycle scan walks a drawn control through its module to the `used` leg only, so no ordinality verdict changes. Every existing control leg in the test suite was migrated to `use_ratio: 0` and passes unchanged; the schema v5 draft is updated to match, which closes its open question on whether a control is consumed.

### Changed — breaking for JavaScript consumers

- **Breaking for JavaScript consumers — the WASM build no longer uses Emscripten.** The kernel is compiled by clang for `wasm32-wasip1` against wasi-libc, from the pinned WASI SDK (`make wasm` fetches it into `tools/` and checks its SHA-256, natively on macOS and Linux). **`dist/gssk.js` is now a hand-written, dependency-free ES module** in place of Emscripten's generated factory, and `package.json` declares `"type": "module"`. File names are unchanged.

  The module keeps the shape callers used: `_GSSK_*` functions, `_malloc`/`_free`, `HEAPU8`/`HEAPU32`/`HEAPF64` and `stringToUTF8`/`UTF8ToString`/`lengthBytesUTF8`. Migration is the import — `import createGSSK from './gssk.js'` — and dropping `ccall`, `cwrap`, `allocate` and `ALLOC_NORMAL`, which nothing in this repository used. `createGSSK({ wasm, print, printErr })` takes the module's location or bytes. A host import the loader does not provide is refused at load, by name.

  **Why:** an archived kernel must stay runnable after today's toolchains are gone. What ships is now a standard `.wasm` and under 200 readable lines of JavaScript, and the same sources build bit-for-bit identically under Guix on aarch64 and x86_64 (`spike/guix`, PR #23).

### Changed

- **TODO.md tidied.** Nine items finished by earlier work are now ticked, each naming what did it:
  the relational algebra, network coupling, the projection and coverage report, Odum's fourth rule,
  carriers, forcing, and the module/harmony-row decision. `mop_network_coupling_requirements.md`
  and the assessment's emergy note now say they are met.

- **All 88 requirements are `implemented`.** The BR-* acceptance lines were checked against `main`,
  and NFR-DET-002, NFR-POR-001, IF-API-002 and IF-API-004 were flipped. The README's status note now
  describes what the second engine implements, where it previously said it implemented none of it.

- **The harmony constructor moves to `src/harmony.c`, renamed `gia_harmony_assume_*`** (FR-HAR-004).
  - Seven functions are renamed: `gia_ordinal_root`, `gia_harmony_init`, `_free`, `_at`,
    `_reconstruct`, `_row_residual` and `_reduction_residual` each become
    `gia_harmony_assume_<name>`.
  - Its output stays labelled `assumed`.
  - No solver or detector calls it (`check-symbols`, T-HAR-01).

- **The kernel's matrix-exponential method is `"expm"`** (FR-KER-001, ADR 0022).
  - `"incipient"` still loads. It runs the same code path, with byte-identical output, and now prints
    one notice per load on stderr saying it is the matrix exponential and naming `expm`.
  - A model serialises back the spelling it was loaded with.
  - The schema enum gains `"expm"`.
  - `make test-kernel-method` (T-KER-01) is in CI.

- **ADR 0015's emergent component is retired.** The generative step no longer adds an
  `emergent_quality_N` node with `emerged_from` and `ordinality_rank`, and it no longer adds
  `ordinal_ascent` or `emergent_feedback_loop` legs. Its tests are retired or revised, as vv-plan.md
  §5 sets out.

- **The generative step follows the Maximum Em-Power Principle** (FR-ORD-005, ADR 0021 §2,
  [02 Eq 5.3]). Below Maximum Ordinality it adds one `linear` pathway at a time, from a sink SCC to a
  source SCC of the component graph. Each time it takes the candidate that maximises total empower
  at `t_end`, and it stops at strong connectivity.
  - `input.json` now gains `consumer_1 → store_1`.
  - `closed_loop.json` is already at maximum and gains nothing.
  - Closure no longer decides anything (FR-ORD-004).
  - `label.generative_step` is now `implemented`.

- **`gia_ordinality` is deprecated, and returns `gia_closure`**: the fraction of components (habitat
  excluded) on a closed pathway. Two consequences:
  - `closed_loop.json` reads closure 1.000 (it was 0.800 with the heat sink counted), and is at
    Maximum Ordinality.
  - Two disjoint 2-cycles have closure 1 but are not at maximum.

  `label.ordinality` is now `implemented`, and `label.closure` is `proxy`.

- **Drift on network trajectories, from the sources: solution drift and output-projection drift** (PLAN.md W2 `idc-drift-coupled`, FR-IDC-011, FR-IDC-014, IF-OUT-001). Breaking for consumers of `giannantoni_sim`'s CSV.
  - **Removed:** the `_idc`, `_tdc`, `_drift` and `psi_network` columns. They applied the drift identity to a per-node φ the engine invented, which was decoupled from `_Q` (PLAN E4).
  - **Added:** a `_drift_proj` column per node, the [09 Eq 13] output-projection drift `(Q'' − Q'²/Q)·Δ²/2` along the solved trajectory. `Q''` is exact (`A²x`) for a constant flow matrix, and Richardson-extrapolated otherwise (numerics N7).
  - **`gia_solution_drift`:** reports solution drift as exactly zero for a constant flow matrix [06 §4 (i)] and refuses it for any other network, which the sources do not define.
  - **Run report:** the per-node "calculi agree / TDC drifts" table and the duet section are gone, since both judged the invented φ. The report now states the solution drift.

- **`AGENTS.md` and `TODO.md` no longer equate the incipient calculus with the matrix exponential** (PLAN.md W1 `guard-agents-wording`, B5 and B6). The changes:
  - The Overview now describes both engines. It no longer says "Euler or RK4" for a repository with four integrators, and it labels the kernel's `"incipient"` method as classical.
  - A new "Giannantoni work: the protocol" section summarises ADR 0018's G1–G7.
  - The clamp `Q < 0 → 0` is scoped to the kernel, since the Giannantoni units never clamp signed or complex coordinates.
  - The `(void)`-stub allowance no longer covers the Giannantoni API.
  - `make test-update` is marked kernel-only.
  - In TODO.md, the Vision, the Non-Goals and the Phase 1 heading stop calling the matrix exponential IDC.

- **A coverage gate on the Giannantoni units, and a kernel gate that can fail** (PLAN.md W1 `guard-coverage-giannantoni`, NFR-COV-001). `make coverage-gia` builds the Giannantoni units with `--coverage`, runs every Giannantoni suite and the CLI system tests, and requires ≥ 90% of lines in `engine.c` (and in `idc.c`, `mop.c`, `relational.c` and `harmony.c` once they exist). `validation.c` and `projection.c` are reported but not gated. The gate is `scripts/coverage_gate.sh`. Its parse is POSIX `awk`, and an unreadable or empty report **fails**: the old kernel gate printed `OK` whenever its `grep -oP` matched nothing, which on macOS was always (PLAN §1 B3). The kernel gate now uses the same portable parse and also fails on an unreadable summary. `make coverage-check` runs both gates. The README badge claimed a 35% gate was ≥ 85%; it now shows both real gates. `make test-coverage-gate` self-tests the gate in CI.

- **`make test` fails on a model with no golden file** (PLAN.md W1 `guard-no-skip`, ADR 0018 rule 1). It used to print `SKIPPED` and exit 0, so a new model "passed" by having no expected output (PLAN §1 B1). A model may now skip only if `tests/skip_allowlist.txt` names it with a reason after `#`. An entry without a reason, or an entry for a model that does have a golden file, also fails, so the list cannot grow silently. `invalid_model` is the one entry: the kernel deliberately refuses it. `make test-guard-no-skip` (`tests/guard_no_skip.sh`, in CI) checks that the guard itself still bites: a golden file removed from a copy of `tests/expected` must turn `make test` red.

- **The Giannantoni engine no longer claims to implement IDC or the MOP, and labels every output** (PLAN.md W0, `mop-claims-remediation`). The README, `include/engine.h`, the `src/gssk.c` header and `docs/giannantoni_assessment.md` said or implied that the repository implements Giannantoni's Incipient Differential Calculus and Maximum Ordinality Principle; the trajectories are the classical matrix exponential, the harmony matrix is assumed, ordinality is a proxy and the generative step is a heuristic (PLAN §2, E1–E8). The assessment now opens with a status-of-record table, TODO §10.3 is marked as testing a constructor, and `giannantoni_sim` prints one `label.<output>: implemented|classical|assumed|proxy|illustrative` line per CSV column and per reported quantity (FR-OUT-001). `make test-mop-cli` (`tests/mop_cli.sh`, in CI) checks the labels against the header the same run writes. No numbers change.

- **Daily work no longer waits on Guix.** `make dev` builds `gssk.wasm` natively with the WASI SDK and runs Vite in the `node:22-slim` container `make test-wasm-container` already uses: seconds, not a Guix fetch. The Guix-backed page is `make dev-guix`. The Guix CI build (`guix.yml`) no longer runs on every PR that edits the `Makefile`, which cost 22–27 minutes per PR; it runs on `guix/` changes, monthly, on demand, and on every version tag, where it still gates the release.

- **Zenodo is disabled until this repository has its own record.** The DOI in the README badge, the citation section and `CITATION.cff` was `energese-project/GSSK`'s, not this repository's; it is commented out (marked `ZENODO:`) to be restored with the new record's DOI, and `.zenodo.json` is parked as `.zenodo.json.disabled`. Whether a release deposits is set by the repository's toggle in Zenodo's GitHub settings, not by these files.

- **The Guix toolchain pack is no longer attached to releases.** At ~450 MB it is produced on demand instead (`guix.yml`, run by hand with `pack-toolchain`). Releases still ship the Guix-built `gssk.wasm`, its hashes and the recipe.

### Added

- **`bin/gia_bridge`: one GSSK-schema model through both engines** (FR-BRG-001).
  - It steps the model with the kernel (`--method euler|rk4|expm|auto|adaptive`), solves the
    model's projection with the Giannantoni engine on the same grid, and reports the largest
    difference per node. `--csv` writes both trajectories and their difference.
  - The difference is labelled `classical`: it is the kernel integrator's error against the exact
    exponential, not incipient drift.
  - A partly projected model is refused (exit 2).
  - `make test-bridge` (T-BRG-01, in CI) checks the Euler, RK4 and exponential cases against
    their closed forms.

- **Three example seeds for the Giannantoni engine**, each checked by hand in `tests/mop_cli.sh` and
  `test_mop.c`:
  - `trophic_chain.json`: the network harmony verdict is computed.
  - `coproduction.json`: the ordinality record shows ½ and 2 couples, and the generative step adds
    three pathways.
  - `harmonic_couples.json`: `--mop-out` writes `R_H = 0`.

  Every seed now runs through the CLI system tests, and `docs/results/harmony_verdicts.md` lists
  them all.

- **`make check-errata`** (VAL-08, in CI): every erratum in PLAN.md §5 has a test assertion naming
  it. X6 and X9 gain theirs.

- **The harmony detector** (FR-HAR-001…003, FR-OUT-003, VAL-07).
  - `gia_harmony_residual` computes R_H of [23 Eq 5.6.5]. `gia_harmony_verdict` runs vv-plan.md
    §6's procedure, and `gia_harmony_observed` gives the verdict for an observed Matrioska.
  - Three constructions give the verdicts R7 derives: First Equation `transported`, Second Equation
    `imposed`, EQS `imposed`.
  - The run report prints one `harmony.<construction>: <verdict>` line per construction.
  - The MOP CSV gains its `R_H` column.
  - `docs/results/harmony_verdicts.md` records each example seed's verdicts, and `mop_cli.sh` holds
    the report to it.
  - `bin/test_mop_emergence` is linked without `harmony.o` (`make test-mop-emergence`, in CI).

- **Boundary conditions from the network** (FR-MOP-008, PLAN R8).
  - `gia_mop_network` sets `e^{α_ij}` to the emergy that the direct pathways `i → j` carry
    (`gia_emergy_carried`, using the emergy pass's partition and replicate rules).
  - It derives `β_ij = α'` by [23 Eq 5.5.2] with k = 1.
  - Couples with no pathway are unrelated.
  - `"beta": "network"` with k = 1 now writes the MOP CSV. With any other k it is refused.

- **Ordinality as ADR 0021 defines it** (FR-ORD-001…003).
  - `gia_ordinality_record` classifies each couple of components as 2/2 (mutual reachability along
    quantity-carrying legs), 2 (both feed one interaction module) or ½ (co-products of one
    replicating process), and counts them as `{k, n₂₂, n₂, n½, n_unrelated}`.
  - `gia_at_maximum_ordinality` now means every couple is 2/2: strong connectivity of the component
    graph.
  - The run report prints `ordinality:`, `maximum_ordinality:` and `closure (proxy):`.

- **The seed's `mop` block, `--mop-out`, and exit codes** (IF-JSON-001, IF-CLI-001, FR-OUT-002,
  NFR-ROB-001).
  - `gia_mop_seed_load` parses the block strictly. Each of these is a load error naming its JSON path:
    an unknown or repeated key, a wrong type, a non-component id, a duplicate couple, and samples
    that do not start at 0 or do not increase.
  - `giannantoni_sim --mop-out PATH` writes `time` and each couple's `_re`/`_im`, in
    `(from, to)` id order.
  - The CLI exits 1 on a load error and 2 on a refusal. A refusal's reason, naming its source, goes
    to stderr.
  - `tests/mop_fuzz/` is a committed corpus of 47 hostile and valid seeds (T-ROB-01).
  - `make bench-mop` times N = 64 at 1,000 output times against NFR-PERF-001's 2 s bound.

- **The Second Fundamental Equation's printed solution** (FR-MOP-005, [23 Eq 6.1–6.3]).
  `gia_mop_second` evaluates `A(t) = α₁₂(0)·w + ln(c₁ + c₂t)`, the specular
  `B = [[A, −A], [−A, A]]`, and the Matrioska row `r₁ⱼ = (e^B)₁₁ w^{j−2}`, with `e^B` in closed form.
  `w = e^{2πi/(N−1)}` is this kernel's reading of the undefined ordinal power in Eq 6.3 (PLANLOG).
  `c₁ + c₂s ≤ 0` on `[0, t]` is refused. T-MOP-07 checks that `u = Ȧ` solves `u' + u² = 0`.

- **The EQS operative form, `gia_eqs`** (PLAN.md W5 `mop-eqs`, FR-MOP-006). It evaluates `ρ₁ⱼ`, `φ₁ⱼ` and `θ₁ⱼ` of [23 Eq 7.1–7.5] for each root `l`, computing the three brackets as the relational product of the De Moivre root `(B, C, C)` with the reference couple's coordinates. Unequal j and k angles are refused, because [23 Eq 7.4] writes one angle for both (erratum X11). The construction is tagged as harmony-assuming [23 §8 ii]. A validation test reproduces the source's printed brackets from the table product over 1000 seeded draws, to 4.4e-16 (VAL-02).

- **The Relational Space algebra, `rel_*`, and relational-valued couples** (PLAN.md W5 `mop-relational-algebra`, FR-REL-001…004, FR-MOP-007; ADR 0020).
  - **Product.** `rel_mul` is the product table of [23 Eq 5.1.3–5.1.5] as printed, which is commutative and non-associative. `rel_mul3` evaluates left to right, so `(j∘j)∘k = −k` while `j∘(j∘k) = k`.
  - **Exponential.** `rel_exp` is the De Moivre form [23 Eq 5.1.2], never a power series.
  - **Roots and powers.** `rel_root` and `rel_root_pow` build the ordinal roots and raise them by multiplying the angle, so that `r^{N−1} = 1`. `rel_mul_pow` is the table power, where the same cube is `(0.540721, 0, −0.665721)` (erratum X10).
  - **Relational couples.** `gia_mop_couple_rel` solves the First Equation componentwise on relational elements for k = 1. It refuses k > 1, because the sources define no division or non-integer power in this algebra.

- **The MOP's First Fundamental Equation, `gia_mop_*`** (PLAN.md W4 `mop-first-equation`, FR-MOP-001…004, NFR-NUM-001…006). New `include/mop.h` and `src/mop.c` solve `(d̃/dt)^k α = β` per couple, with `α(0) = 0`. The solution is the one derived from [23 Eq 5.5.6], `α = {(1/k)∫β^{1/k}}^k`; the printed [23 Eq 5.5.7–5.5.8] leave residuals of 1.25 and −0.625 (erratum X1, ADR 0019).
  - **Affine-power β:** numerics N1's unified form. It holds 1e-12 relative error down to `bt/a = 1e-11`, where the naive closed form loses 5.9e-5, and is continuous through `b = 0`.
  - **Sampled β:** adaptive Gauss–Kronrod 7–15 on a phase continued across the negative real axis.
  - **`gia_mop_solve`:** builds the Matrioska, marking couples with no boundary condition as unrelated rather than zero.
  - **Domain:** the First Equation's domain (FR-MOP-002) is refused by name. Overflow is `GIA_E_RANGE`, never inf, and nothing is clamped.
  - **`gia_quad_gk15`:** the integrator, exposed so that its error estimate and its refusal are testable.

- **The ordinal forms of the emergy processes, and the circle product** (PLAN.md W3 `emergy-ordinal-forms`, FR-EM-005, FR-EM-006).
  - `gia_oform_*` builds co-production as a binary (a column of two `Em(u)` branches), interaction as a duet (a row `[Em(u₁), Em(u₂)]`), and feedback as the specular duet-binary `[[a₁, a₂], [a₂, a₁]]` [22 Eq 6–8].
  - `gia_circle_product` keeps every pair of factors, so `(a₁; a₂) ∘ [b₁, b₂] = [(a₁b₁; a₂b₁), (a₁b₂; a₂b₂)]` before reduction [06b Eq 2]. `l ∘ l` is the du-et `[l, l]`, not `l²` [02 Eq 14.11.5].
  - `gia_circle_reduce` is the cardinal reduction that maps each pair to its product (PLAN R12).

- **Emergy source terms, the global balance, and limits that refuse rather than truncate** (PLAN.md W3 `emergy-source-terms`, FR-EM-004, FR-EM-007, NFR-LIM-001).
  - `gia_emergy_source_term` returns a process's equivalent source term, emergy out minus emergy in [02 Eq 3.6–3.17]: `(n−1)·Em(u)` for a co-production with n products, 0 for a partition, and 0 for an interaction whose inputs are all drawn.
  - `gia_emergy_global_balance` and `gia_emergy_balance_solve` evaluate [02 Eq 3.21] and solve it for the source terms. On the totals of [02 Fig. 3.4] they reproduce `Φ_D = 30,000` and `Φ_E = 15,000` (VAL-04).
  - **Fixed:** the emergy pass silently dropped a node's inflows past the 64th, and ignored co-production ancestry past 64 nodes. It now refuses with `GIA_E_LIMIT` (`gia_emergy_check_limits`), `gia_emergy_at` returns false, and the trajectory CSV leaves the emergy cells empty.

- **The nonlinear equation of [02 Eq 14.10.1], `gia_nl1410_roots`** (PLAN.md W2 `idc-nonlinear-14-10`, FR-IDC-009). Substituting `F = e^{ut}` reduces `F·(d̃²/dt²)F² + A F²·(d̃/dt)F + B F³ = 0` to `4u² + Au + B = 0`. The solutions are its two roots, not the "triplet" [02 Eq 14.10.2] prints (erratum X8). The test forms the equation's incipient residual with `gia_idc_of` for each root, and checks that no third candidate solves it.

- **The incipient Taylor projection, `gia_idc_taylor`, validated against Giannantoni & Zoli 2009** (PLAN.md W2 `idc-taylor`, FR-IDC-010, VAL-01). `f*(t₀+Δ) = f(t₀) Σ_{k≤n} (aΔ)ᵏ/k!` with `a = f'/f` [09 Eq 10], evaluated by Horner. With n = 2 it reproduces [09]'s published 16.4, 3.01, 178.0, 172.06 and the 15–17 cm range. The tests also pin errata X4 (Eq 19 at τ₀ = 2 gives 156.0, not the printed 154.3) and X5 (the minimum scenario's net increase is 2.6125, not 1.91).

- **The Riccati equation by linearisation, and named refusals** (PLAN.md W2 `idc-riccati`, FR-IDC-008 and FR-IDC-013).
  - `gia_riccati_solve` and `gia_riccati_eval` solve `f' + Qf + Rf² = P` with the substitution `f = y'/(R y)`, taking y from the incipient LDE solver ([06 Eq 3.16–3.18], PLAN R11). They report the traditional Riccati residual of the result.
  - The printed substitution of [06 Eq 3.17] is inverted. A test shows that it leaves an Eq 3.18 residual of 18.3 against a tolerance of 1e-6 (erratum X2).
  - `gia_idc_refuse` refuses the Riccati duet [06 Eq 3.22], Abel's n-et and network solution drift, each with the reason the sources give (PLAN §6).

- **The binary function, `gia_binary_*`** (PLAN.md W2 `idc-binary`, FR-IDC-007). It solves `f' + A f^(½) + B f = 0` on both branches [02 Eq 14.7.1–14.7.7] from the single characteristic `u² + Au + B = 0` (PLAN R9; the printed per-branch exponents are erratum X7). The constants come from the four initial conditions. A double root is refused, because no source defines that case. Tests check both branches by hand for A = 1, B = −2. They also check the defining residual with `f'` by central difference, a complex case, linearity in the initial conditions (FR-IDC-012, now `implemented`), and an overflow returning `GIA_E_RANGE`.

- **The second-order incipient LDE with variable coefficients, `gia_lde2_*`** (PLAN.md W2 `idc-lde2`, FR-IDC-006). It solves `f̃'' + a₁(t) f̃' + a₀(t) f = 0` as `Σ cᵢ exp(∫₀ᵗ rᵢ)`, with `rᵢ(t)` the roots of the incipient characteristic [06 Eq 3.3–3.6], labelled by continuity in t.
  - **Method** (numerics N3): the stable root form, plus an adaptive Gauss–Legendre sweep that accepts a step only when its integral agrees with two half-steps to 1e-12 and no root moves more than half the root separation.
  - **Double root:** a root that is double throughout gives the single family `c·e^{∫r}`. Initial conditions inconsistent with that family are refused, and [06 Eq 3.7] is not used (erratum X12).
  - **Coincident or colliding roots:** roots that coincide at t = 0 only are refused, because the initial conditions cannot fix both constants. Roots that collide later are refused with the collision time.

  Tests check the solution by hand for `a₀ = −(1+t)²`, `−(1+sin t)²` and complex `−(1+it)²`. They also check the termwise incipient residual (0), the traditional residual (`c₊e^{φ} − c₋e^{−φ}`, not 0), the worked example of [10 App. Eq 28–34] (VAL-05), and linearity in the initial conditions. The baseline's T-IDC-04 (i) is amended, because its roots coincided at t = 0 where N3 itself refuses (PLANLOG).

- **The incipient derivative of a general function, `gia_idc_of`** (PLAN.md W2 `idc-general-f`, FR-IDC-002). New units:
  - `include/gia_status.h` carries the Giannantoni status codes and the `why` contract (IF-API-001).
  - `include/idc.h` and `src/idc.c` hold the single-variable incipient calculus. `gia_idc_of` computes `(f'/f)^n f` [02 Eq 14.9.5] over ℂ, with integer powers by multiplication. It refuses `f = 0` as outside the domain and an overflow as out of range, and never returns a non-finite value.

  The new V&V suite `make test-mop` (`tests/test_mop.c`, in CI) checks `f = 1 + t²` against `(2t/(1+t²))ⁿ(1+t²)` for n = 0…8 at 1e-12 (T-IDC-03), and the status contract (T-API-01). `idc.c` joins the coverage gate, the symbol checks, the API check and the sanitizer run.

- **ADR 0022: the kernel's matrix-exponential method is named `expm`** ([docs/adr/0022](adr/0022-kernel-method-expm.md)). `"method": "incipient"` runs a Padé matrix exponential, not Giannantoni's incipient calculus (PLAN.md §2 E2). Under ADR 0011 the kernel is the classical engine. `"expm"` becomes the documented name. `"incipient"` stays a deprecated alias with byte-identical output and a one-line notice. This is not a schema break. No behaviour changes until `kernel-method-label` lands.

- **ADR 0021: ordinality as Giannantoni defines it, and a generative step under Maximum Em-Power** ([docs/adr/0021](adr/0021-ordinality-and-the-generative-step.md)). It supersedes the cycle-coverage definition of ordinality in ADRs 0014 and 0015, and ADR 0015's emergent component. ADR 0014's leg rules stand.
  - Each couple of components is classified as 2/2, 2, ½ or unrelated from the exponents in [10] and [22 Eq 6–8]. Maximum Ordinality is "every couple 2/2" [22 §12.1], which is strong connectivity of the component graph. Boundary nodes are habitat, not components.
  - The generative step adds one linear pathway at a time, from a sink to a source of the condensation, choosing the candidate that maximises total empower [02 Eq 5.3]. It repeats until the graph is strongly connected.
  - Cycle coverage survives only as `closure (proxy)`.
  - The ADR lists each example seed's verdict before and after. `closed_loop.json` moves from "below maximum (0.800)" to "at maximum".

  No behaviour changes until `mop-ordinality` lands.

- **ADR 0020: the Relational Space algebra** ([docs/adr/0020](adr/0020-relational-algebra.md)). The engine will use the product table of [23 Eq 5.1.3–5.1.5] as printed, which is commutative and non-associative. Products of three or more factors are evaluated left to right and never reassociated. The evidence is that the literal table reproduces the EQS brackets of [23 Eq 7.1–7.3] to 8.9e-16, while an associative reading does not. Exponentials are De Moivre's closed form, never a power series. A root's power multiplies its angle, and the table power is a separate function: under the table the "roots of unity" are not roots of unity (erratum X10). The EQS refuses unequal angles (X11). The circle product of [02], [06b] keeps factor pairs and reduces cardinally to multiplication. No behaviour changes.

- **ADR 0019: the MOP's two Fundamental Equations, as the engine solves them** ([docs/adr/0019](adr/0019-mop-fundamental-equations.md)). It records three decisions:
  - **First Equation:** solved per couple with the solution derived from [23 Eq 5.5.6]. The printed [23 Eq 5.5.7–5.5.8] leave residuals 1.25 and −0.625 at t = 1 (erratum X1).
  - **Second Equation:** implemented as its printed solution [23 Eq 6.1–6.3], tested against the Riccati equation `u' + u² = 0` that the solution satisfies. That oracle is labelled reconstructed.
  - **Harmony:** reported as `transported`, `imposed` or `absent` from a measured residual, never assumed by the solver. The sources determine the three verdicts: the First Equation transports harmony, and the Second Equation and the EQS impose it.

  Emergence of harmony is recorded as not reproducible from the sources. No behaviour changes.

- **Every Giannantoni API function must be called by a test, and none may be a stub** (PLAN.md W1 `guard-api-called`, NFR-API-001, ADR 0018 rule 3). `make check-api-called` (`scripts/check_api_called.sh`) parses every function declared in `gia_status.h`, `idc.h`, `relational.h` and `mop.h`, plus the functions IF-API-005 adds to `engine.h`. A declaration that no test names fails, and so does a not-implemented marker in `idc.c`, `relational.c`, `mop.c` or `harmony.c`. `AGENTS.md`'s `(void)`-stub allowance never covered these headers, and now cannot. `make test-api-called` runs the check against the C blocks of `docs/requirements/icd.md` itself: it must find exactly the 27 functions declared there, and must name each one when it is left uncalled. Both run in CI.

- **Reentrancy, memory and separation checks for the Giannantoni units** (PLAN.md W1 `guard-reentrancy`). Four checks, all in CI:
  - `make test-mop-threads` runs two models on two threads and requires each to match its solo run (T-REE-01).
  - `make check-symbols` (`scripts/check_symbols.sh`) fails on a writable data symbol, decided by section rather than by `nm` letter so `const` tables pass. It also fails on a call to `exit`, `abort` or `assert`, and on any `#include "gssk.h"` (T-REE-02, T-ERR-01, INS-SEP-01).
  - `make test-mop-asan` runs the Giannantoni suites under ASan, LSan and UBSan, with every finding fatal (T-MEM-01). It runs in the Linux quality gate, because LeakSanitizer does not run on macOS.
  - `tests/mop_cli.sh` checks that two runs give byte-identical output (T-DET-01).

  NFR-REE-001, NFR-MEM-001, NFR-ERR-001, NFR-SEP-001 and NFR-DET-001 are now `implemented`. `make test-mop-threads-tsan` runs the thread test under ThreadSanitizer where the toolchain has its runtime.

- **ADR 0018: the test protocol for Giannantoni work** ([docs/adr/0018](adr/0018-giannantoni-test-protocol.md)). Records PLAN.md's guardrails G1–G7 as a binding decision, together with three rules that close the remaining loopholes from PLAN §1: no golden files for the Giannantoni units; a coverage gate; no stubs in the planned API; the errata protocol (a probe, a PLAN §5 row, and a test that the printed form fails); definitions that cite sources; no clamping of signed or complex coordinates; TDD evidence in every code PR; refusals bounded by PLAN §6; a suite counts only if CI runs it; and a tautology is named as one. No behaviour changes.

- **A requirements baseline for the Giannantoni kernel, with traceability enforced in CI** ([docs/requirements/](requirements/README.md)). It holds stakeholder requirements, an SRS, an interface control document, a numerical design and a V&V plan: 88 requirements, each traced to a source equation, and an 81-entry test catalogue giving each test's oracle and the mutation that must fail it. `make check-trace` runs in CI and fails on an untraced requirement, an orphan `Verifies:` tag, or a requirement marked `implemented` that no test backs; ten existing tests are tagged. Writing the requirements found two new errata in the sources ([06 Eq 3.7], [23 Eq A2.6]), a numerically unstable closed form near t = 0, and two engine defects now tracked as requirements: file-scope mutable state in the generative step (`src/engine.c:2830`) and silent truncation past 64 inflows. `PLAN.md` moves to revision 3, with the reasoning in `PLANLOG.md`. No engine behaviour changes.

- **`make docs-build-container`: the CI documentation build, run locally.** It mirrors the "Documentation builds" job in `deploy.yml` step for step — Node 20 from a container (the Apple `container` CLI), `npm install`, `npm run docs:build`, and the two checks that the theme's design tokens were inlined — so a docs change can be verified on a host without Node before it is pushed. Listed in `CLAUDE.md` beside the other containerised targets.

- **A plan for implementing Giannantoni's IDC and Maximum Ordinality Principle, test-first, from the sources** ([PLAN.md](https://github.com/energese-project/GSSK-Giannantoni/blob/main/PLAN.md), reasoning in [PLANLOG.md](https://github.com/energese-project/GSSK-Giannantoni/blob/main/PLANLOG.md)). It records what the Giannantoni engine actually does today, closes every design question from the 2002–2023 sources, registers ten errata in those sources with the probe that decides each (`docs/sources/probes/`), and audits the repository's own process for ways a build can go green without the mathematics being implemented. Extracted text and equation pages of the two CC BY 4.0 papers are under `docs/sources/cc-by/`. No engine behaviour changes.

- **Divide and subtract are actions of the interaction module** ([ADR 0016](adr/0016-interaction-actions.md)), so Odum's Fig. 2.6(d) and (e) can be written as modules and not only as the `ratio` and `subtract` pathway laws ADR 0012 is retiring. A work gate's `module` block takes `action`: `multiply` (the default, so every existing gate is unchanged), `divide` (`F = k · Q_energy / max(Q_control, ε)`) or `subtract` (`F = max(0, k · (Q_energy − Q_control))`). Both new actions take exactly one control; an unknown action, a non-string one, or `action` on any module but an interaction is a load error. The subtracting action's crossing is located as an event, as a switch's is, so each side is solved exactly. Tests check each module form against its closed form and against the pathway law it replaces, and that its control is drawn under ADR 0017 — including that a clamped subtract draws nothing.

- **`make dev`: run any example model in the browser, on the release toolchain.** Vite serves `dev/`, a page that runs `examples/*.json` through `dist/gssk.js` and charts and tabulates the result, from a long-lived container hosting Guix. Node comes from the pinned Guix (`guix/dev.scm`) and `gssk.wasm` from the release recipe (`make wasm-guix`), so the page runs the bytes a release ships; the host needs only the `container` CLI. This is the loader's first browser test: headless Chromium ran all 23 valid examples, and `invalid_model.json` showed the kernel's own error message. `package-lock.json` is now committed, pinning npm packages by hash.

- **Tagged releases ship a `gssk.wasm` built by Guix, reproducible bit for bit.** [`guix/`](https://github.com/energese-project/GSSK-Giannantoni/blob/main/guix/README.md) pins Guix (`channels.scm`) and builds wasi-libc, wasm32 compiler-rt and the kernel from pinned source (`gssk.scm`), through the Makefile's own `make wasm` and `make test-wasm`. On each version tag CI builds it on x86_64 and aarch64, requires identical bytes and a passing `guix build --check`, and releases that `gssk.wasm` with `gssk-guix.sha256` and the recipe. On demand (`pack-toolchain`), it also packs the toolchain — clang 21, wasi-libc, make and Node, relocatable so `make wasm-toolchain TC=…` rebuilds the identical file on any x86_64 Linux without Guix; at ~450 MB it is not attached to releases for now. Verified before adoption: unpacked in a Debian container with no Guix, the toolchain rebuilt the same bytes and passed all 37 WASM tests. Rolling `latest` keeps the faster WASI SDK build.

  The Makefile's WASM toolchain is now overridable (`WASM_CC`, `WASM_SYSROOT`, `WASM_TOOLCHAIN_FLAGS`, `WASM_TOOLCHAIN_LIBS`), so one rule defines the build for the SDK, Guix and the archived toolchain alike.

- **The WASM artefact is tested like the native one.** `make test-wasm` (or `make test-wasm-container` without host Node) runs three suites through `dist/gssk.js` under Node's built-in test runner: the loader's contract, **every regression model against `tests/expected`** — before, only the forcing evaluator ran under WASM, so a kernel change could break the shipped artefact for every other model unnoticed — and the forcing-parity check, ported from `forcing_parity.cjs` to an ES module. The repository's JavaScript is `.js` (and `.ts`) ES modules only — no `.mjs` or `.cjs`; the VitePress config moves from `config.mts` to `config.ts` accordingly. The CI job's grep for export names is replaced by a test that each declared export is a callable function.

- **The Giannantoni engine's transactor: Odum's §XV exchange as a module with four named legs.** [ADR 0013](adr/0013-module-roles-not-position.md) decision 5.

  `goods_in`, `goods_out`, `counter_in` and `counter_out` are all named and all required; a transactor touched by an unnamed pathway is an error. Position could not serve here and neither could carrier identity — given a goods flow and two money stocks, which one is the payer is ownership, not carrier. Eq (103) is recovered from the trajectory instead of asserted: the two flows counter-flow and their ratio is the price. Because a transactor is linear in the goods it moves, the flow matrix stays constant and `ψ` is exactly zero, unlike a work gate. Barter — all four legs on one carrier — stays legal; a single leg pair straddling two carriers does not.

  **A module holds no carrier**, being a hyperedge drawn as a symbol rather than a stock, so it contributes no carrier class and is exempt from the pairwise cross-carrier check on pathways. Which of its legs must agree is the module's own rule, which also lifts a restriction that was never intended: a work gate metered by a price or a population is an ordinary model, and the pairwise check refused it.

### Removed

- **Makefile targets for trees this fork does not carry.** `test-python`, `demo-python`, `plot-demo`, `demo-native`, `container-image-demo` (with `Containerfile.demo`) and the LaTeX targets `doco`, `whitepaper`, `article`, `conformance`, `doco-clean` all ran against `python/` or `doco/`, neither of which came across from GSSK, so every one of them failed. `make demo` stays, now native and plot-free, printing what the README already described.

### Fixed

- **The harmony constructor's matrix is sized by components, not nodes** (ADRs 0014, 0021). It had one
  row per node, including modules and habitat, under a "components N" label: `input.json` reported
  N = 5 for its two components.

- **CI now runs three suites it had been skipping**: `test-mop-emergence`, `test-kernel-method` and
  `test-coverage-gate` were in `CI_TESTS` but had no step in `deploy.yml`.
  `make check-ci-matrix` (DEM-POR-01, in CI) now fails on any such gap, on a missing toolchain, and
  on `-Werror` leaving `CFLAGS`.

- **The generative step is reentrant** (PLAN.md W1 `guard-reentrancy`, NFR-REE-001). `gia_generate` sorted its open components through a file-scope `static const gia_model *sort_model` read by the `qsort` comparator, so two models generating on two threads could sort by each other's ids. ThreadSanitizer reported the race at `src/engine.c:2918`. Each sort element now carries its id. Two other writable statics, found by the new structural check in a `-O0` build, are now `const`: `RULE` in `src/validation.c` and `prim` in `src/projection.c`. Output is unchanged.

- **`make test-linux` and `make test-linux-clang` run what CI runs.** Each spelled out its own suite list, and both had drifted: gcc was missing `test-limit-logic`, `test-reversible`, `test-node-type-enum`, `test-edge-flows`, `test-giannantoni`, `check-version` and `test-schema`; clang ran only `test` and `test-advanced`. They now share one list, `CI_TESTS`, so `make ci-local` is evidence of what it claims.

- **The MOP ordinal step could move a system away from maximum ordinality, and never stop.** [ADR 0015](adr/0015-what-the-emergent-quality-closes.md). It wired `hub → E → open`, which closes `open` only if `open` already reaches the hub; for a dead end it closed nothing and added `E` as a new open component, so `a ⇄ b, a → c` went 0.667 → 0.500 → 0.400 → … while printing "closed the loop". Four of five measured models never terminated. It also chose what to close by array order.

  The emergent quality now closes every open component in one step, adding only the missing direction(s) through `E`, so one step reaches the fixed point and a second changes nothing. A sink is never closed or drawn from — that would draw on energy already used — so a model whose only open component is a heat sink stops below maximum and says why. Hub ties are broken by id and legs emitted in id order, so a reordered model appends byte-identical output. Closure is printed only after the evolved graph is reloaded and rescanned; otherwise the seed is returned unchanged.

  **Output format:** the emergent component is a `storage` (it was a `gain`, which since ADRs 0012–0013 is a module and could close nothing under ADR 0014), `emerged_from` is an array, and its id and label are `emergent_quality_N` and "Emergent Quality".

- **Ordinality was decided by legs that carry nothing, and by modules that hold nothing.** [ADR 0014](adr/0014-ordinality-over-quantity-legs.md). Since module-hosted laws made a work gate's control an edge, the cycle scan walked it: a pure accumulator was reported at maximum ordinality, so the MOP step returned the seed unchanged when it had an open relationship to close. A module was also counted as a component, so one system scored 0.667 written as a pathway law and 0.75 written as a module. The scan now skips control legs and passes through modules without counting them — a transactor keeping goods and counter-flow apart — and `gia_generate` shares it instead of keeping its own copy, never choosing a module as the component to close or the hub to draw from.

  **`gia_system_is_closed` asked the same wrong question at the boundary.** A gate metered by a constant made a system conserving to 1e-14 report open, and `--print` explained the zero by saying the constant "delivers quantity". A control moves nothing across the boundary; a source on an energy leg still does.

  New: `gia_component_count`. The `--print` header reports modules separately from components.

- **`gia_edge_flow` read a module pathway's flow from the pathway's own law.** A module's pathways carry no law, so that was the default — linear, weight 1 — and the emergy pass therefore carried transformity along a flow the solver never used. It now asks the module. Introduced with module-hosted laws and caught while adding the transactor.

> The rest of this engine's history is in [odum_1972_conformance.md](odum_1972_conformance.md) and in `TODO.md` at the repository root, which were written alongside it; this is the first entry here because the changelog tracks the GSSK kernel's releases.

---

## [5.3.0] — 2026-09-05

### Added

- **The rest of Odum's work gate: an n-ary `interaction` and a new `subtract` logic.** GIP-0001 G1; the decision is [ADR 0008](adr/0008-nary-interaction-and-subtracting-action.md).

  *Modeling for All Scales* Fig. 2.6 shows one interaction glyph computing several ways — (a) a product of two inputs, (c) a product of **three** input forces, (d) a divisor action, (e) a **subtracting** action. GSSK shipped (a) as `interaction` and (d) as `ratio` (ADR 0002, extended by ADR 0005). This adds (c) and (e).

  **(c) `params.control_nodes`** — an array of one to eight node ids, giving `F = k × Q_origin × ∏ Q_control`. The singular `params.control_node` is unchanged and remains the one-control spelling; the two are **mutually exclusive**, and an edge carrying both is rejected rather than having one silently win. Only `interaction` accepts more than one control: `ratio` and `subtract` are binary and `limit` takes one half-saturation constant, so a second control on any of them is a `GSSK_ERR_SCHEMA_VIOLATION`. A quotient or a difference of three things is not defined by the figure, and guessing an associativity would be inventing semantics rather than implementing them.

  The workaround this replaces — an intermediate storage node holding the partial product — was never the same model: the intermediate integrates, so it lags the inputs it multiplies, and it shows up in the state vector, in the emergy accounting and in every CSV column.

  **(e) `subtract`** — `F = max(0, k × (Q_origin − Q_control))`, requiring exactly one control node, read and never consumed. Not the GIP's proposed `params.op` string: the divisor action already shipped as its own logic type, so `op` would have made `logic: "ratio"` and `logic: "interaction", op: "div"` two spellings of one thing, in a published schema. `logic` is also the branch key at eight switch sites, where `-Werror` turns a missed case into a build failure — an `op` string is checked by nothing.

  The clamp is the semantics rather than a guard. `subtract` draws a *barbed* pathway, and a negative flow would drain the target along a line whose barb says it cannot; the gate saturates at zero when the control overtakes the origin. Signed flow across a gradient is `reversible` (ADR 0007), which reads the **target** and draws without a barb. The two compute a difference and are not the same statement.

  **Solver eligibility.** The n-ary product is linearised at the current operating point exactly as the two-input one always was — no new error class, whether there are one or seven controls. `subtract` contributes four flow-matrix entries like `reversible`, and is **exact** rather than linearised where the difference is positive, contributing nothing below the clamp. It is therefore *piecewise* exact: an edge whose difference crosses zero inside a step is integrated as though it stayed on the side it started, the same seam `threshold` has, bounded by `dt`.

  Existing single-control models are **bit-identical** — `control_idx` still holds the one control and the product loop runs zero times — and existing `GSSK_LogicType` values are unmoved, with `GSSK_LOGIC_SUBTRACT` appended at 7. Both are pinned by test. `build_topology_json` emits the spelling an edge actually has, so a one-control model does not come back rewritten into the array form. New: `examples/three_input_gate_model.json`, `tests/test_interaction_nary.c`, `make test-interaction-nary`.

### Fixed

- **`price_node` now resolves for expanded archetype members.** GIP-0002, raised from `gssk-budget`. An `exchange` inside an archetype could not have a per-instance price: every instance traded at the template's constant, and the mechanism that exists to make price a state variable did nothing. A model built on such an archetype validated, loaded, ran to completion and reported success while every transaction diamond in it moved the wrong money — or, with the template constant left at its `0.0` default, no money at all.

  Two things were broken, one behind the other. `GSSK_Init` resolves `price_node` in a second pass over the top-level `nodes` array, and for a composite the entry there is the *instance* — `{"id": "groceries", "type": "purchase_consumed"}` — whose params carry no `price_node`; the expanded members that do carry it are never visited. Behind that, the archetype template struct had no field for it at all, so a template's `price_node` was already discarded at parse time, before expansion could have rewritten anything.

  A template's `price_node` names a **sibling member by its template-local id**, and is now rewritten during expansion exactly as an edge's `origin` and `target` already were. An archetype with a `constant` member `price` and an `exchange` member declaring `price_node: "price"`, instantiated as `groceries`, expands to `groceries__price` and `groceries__deal` with the latter priced from the former. The per-instance value is then supplied through `snapshot.state` on `groceries__price` — a channel that already existed, already applied, and, being part of the model document, inside the content hash. That last property is why this is a parser fix and not a new `GSSK_SetNodePrice`: a price set after `GSSK_Init` would sit outside the hashed document, and two models differing only in price would hash identically. For a ledger the hash is the audit anchor, so that is a correctness objection rather than a matter of taste.

  **A template `price_node` naming no member of its archetype is now rejected** with `GSSK_ERR_SCHEMA_VIOLATION` naming the archetype, the member, the unresolved id and the instance — it does not fall back to the constant. This is the one place the implementation departs from the proposal, on ADR 0004's rule that a name matching nothing is an error rather than a silent re-modelling; a fallback is precisely what let the original defect run green. Nothing can have depended on it, since the field was discarded at parse time. **Top-level `price_node` is unchanged**, second-pass forward references and constant fallback included.

  Note that a price delivered through `snapshot.state` is live state, not a topology initial condition, so it round-trips through `GSSK_SerializeSnapshot` and not through `GSSK_SerializeModel`, which emits each node's `initial_value` by design. Both forms emit the member's namespaced `price_node` reference.

  New assertions in `tests/test_price_node.c`: two instances of one archetype at different prices spend different amounts and each satisfies `spent / inventory == its own price`, the shared buyer's debit equals the sum of the two credits, the template's poison constant is not consulted, an unresolvable template reference fails `GSSK_Init`, and the snapshot round-trip reproduces both instances' spend.

### Changed

- **The repository moved to the `energese-project` organisation, and the citation metadata now reflects it.** GitHub redirects a renamed org for repository URLs but **not** for GitHub Pages, so `sholtomaud.github.io/GSSK` began returning 404 while `github.com/sholtomaud/GSSK` kept working. Three published links were dead rather than merely stale: the README's **Live Docs** and **Interactive Demo**, the **WASM Demo** entry in the docs site's own nav (`docs/.vitepress/config.mts`), and the demo link in the household example. All URLs across the README, the docs site, `package.json`, `web/index.html`, this changelog's compare links and `docs/gip/` now name `energese-project`. `doco/conformance.tex` is deliberately untouched — it carries a scholarly claim about which release its numbers came from, which is the author's to update.

- **`CITATION.cff` is now a complete Zenodo deposit record.** It carried a title, an author and a version, and nothing else; the minted Zenodo record consequently had no keywords at all and was close to undiscoverable by subject. It now carries the concept DOI, an abstract, keywords, the docs URL and a release date. The README gained a **Citation** section distinguishing the **concept DOI** — `10.5281/zenodo.22339312`, which always resolves to the latest release and is what the badge points at — from the per-release **version DOI**, which is what a reproducibility claim should cite. The badge previously pointed at a version DOI, so it would have kept advertising v5.2.0 after every subsequent release.

- **`CITATION.cff` joined the version-sync gate.** `scripts/check_version_sync.py` held `package.json` and `include/gssk.h` to the same version but knew nothing about `CITATION.cff`, and `scripts/release.sh` never touched it — which is how it came to claim 5.1.0 while the archived release was v5.2.0. Both now include it, so a citation that misstates the version fails the build rather than reaching Zenodo.

### Fixed

- **The v5.2.0 tag was cut without running `scripts/release.sh`**, leaving `include/gssk.h` and `package.json` declaring 5.1.0 at a commit tagged v5.2.0, and the changelog's feature work still sitting under `[Unreleased]`. This release restores the invariant: the header, `package.json`, `CITATION.cff` and the changelog all say 5.3.0. The changelog's link-reference definitions had also stopped being maintained at 3.6.0, so `[Unreleased]` compared against v3.6.0 across three majors; 4.1.0, 5.0.0, 5.1.0 and 5.3.0 are now defined.

---

## [5.1.0] — 2026-08-30

### Added

- **`gssk.schema.json` now ships in the npm package.** `package.json` `"files"` listed `dist/`, `include/` and `README.md`, so the only machine-readable statement of the model vocabulary was not published. Downstream consumers hand-maintained their own copy of the node-type enum instead, with nothing to detect drift against the real one. The schema is now in `"files"`, `make check-version` asserts it stays there, and CI asserts `npm pack` actually puts it in the tarball.

- **`GSSK_NodeType` is public, with `GSSK_GetNodeType` beside the existing string getter.** Most of GIP-0001 G7 was already closed on `main`: `GSSK_GetNodeTypeString` returns the type, the primitive set is a closed enum in `gssk.schema.json`, and `GSSK_Init` rejects an unrecognised type naming the node — G7's acceptance criterion verbatim. What remained is that a C consumer had only string comparison for a decision the kernel makes with an integer. That is slower, and worse, typo-able in a way the compiler cannot see: `strcmp(t, "loop_limted")` is a valid program that quietly never matches.

  The enum already existed privately in `src/gssk.c`, and was already called `GSSK_NodeType`. It has been **moved** to `include/gssk.h` rather than copied, so there is no pair to drift; its internal constants are now `GSSK_`-prefixed, which a public header requires and which `-Werror` verified at all 68 use sites.

  The GIP's proposed enum was incomplete against `main`: the primitive set is nine, not four — `storage`, `source`, `sink`, `constant`, `interaction`, `gain`, `loop_limited`, `exchange`, `switch`. The Phase 7 processing nodes are in it. Composites and archetypes are deliberately **not**: they expand during `GSSK_Init`, so by the time a consumer can ask, every node is a primitive. A node authored as part of a `producer` reports `storage`, `interaction` or `sink` — use `GSSK_GetNodeComposite` / `GSSK_GetNodeRole` to recover where it came from.

  Ordinal values are explicit and pinned by test: they cross the WASM boundary as bare integers, so renumbering them is a silent breaking change for every JS consumer. Append new primitives before `GSSK_NODE_INVALID`.

- **`reversible` edge logic — Odum's barb-less pathway.** `F = k × (Q_origin − Q_target)`, signed. Raised as GIP-0001 G3; the decision is [ADR 0007](adr/0007-reversible-pathway.md).

  Odum draws two pathway kinds and the distinction lives in the notation: a barb where the flow depends only on the force behind it, no barb where it depends on the difference between the force at one end and the back force from the other, *"and this pathway may flow in either direction"* (*Modeling for All Scales*, p.23). Every GSSK logic computed forward from origin quantities and **none of them read the target**, so the entire second class — diffusion, exchange across a gradient, any equilibrating process — was inexpressible.

  `reversible` is the only logic that reads both ends, and therefore the only one that can transport backwards along its declared direction. For a barb-less line `origin` and `target` name the two ends rather than a from and a to: **swapping them produces an identical trajectory**, asserted step-by-step rather than at equilibrium.

  This is not the `exchange` node, which is the transaction diamond — two barbed pathways carrying a money/goods counter-flow. A barb-less pathway is one pathway whose sign is a gradient.

  Two properties worth knowing:

  - **It is exactly integrable.** Being linear in the state, the incipient/IDC solver integrates it exactly rather than linearising about the operating point, unlike `limit` and `ratio`. It stays IDC-eligible. It is also the first logic whose flow-matrix contribution touches **four** entries rather than two, which is why `build_jacobian`'s second-variable slot is no longer named for the control node — `reversible` depends on its *target*, which is not a control.
  - **A backward flow carries no transformity.** The quality pass's `flow ≤ 0` clamp stops being incidental and becomes a decision: a flow running back up a gradient is not producing the node it arrives at, so `GSSK_GetEdgeQualityFlow` reads `0.0` while the pathway reverses. `GSSK_GetFlows` reports the signed rate.

  `GSSK_LOGIC_REVERSIBLE` is **appended** at 6; every existing `GSSK_LogicType` value is unchanged and pinned by test, because the value crosses the WASM boundary as a bare integer where nothing recompiles.

- **`examples/diffusion_model.json`** — three tanks in a line, two reversible edges, equilibrating from a step gradient. Both edges are declared *against* the physical gradient on purpose, so the example exercises backward transport rather than merely describing it: a build that clamped the flow at zero produces a visibly different CSV. No existing golden moved.

- **Per-edge flow rates are readable (`GSSK_GetFlows` / `GSSK_GetFlowCount`).** A consumer could read every node quantity and not one rate: `GSSK_GetState` reports the storages, and nothing reported the flows between them. Flow was step-local — computed inside `compute_derivatives` to build `deriv[]`, then discarded. That makes a diagram whose entire subject is flow impossible to annotate, and leaves the heat-sink budget and any pathway-level emergy display with nothing to read. Raised as GIP-0001 G4.

  The new pair mirrors `GSSK_GetState` / `GSSK_GetStateSize` exactly, so it introduces no new idiom: index `i` is the edge at position `i` in the input JSON, matching `GSSK_GetEdgeID` and `GSSK_GetEdgeK`. Both are exported to WASM and declared in `src/gssk.d.ts`.

  Semantics worth knowing before you read them:

  - Refreshed by **every** `GSSK_Step` and `GSSK_StepAdaptive`, whether or not quality accounting is enabled. The pre-existing per-edge flow array only existed inside the quality pass, was freed at the end of it, and what it exposed through `GSSK_GetEdgeQualityFlow` is `Tr × flow`, not flow.
  - Evaluated at the **post-step** state and time, so it is the rate the step just integrated rather than a prediction of the next one.
  - `0.0` before the first step and after `GSSK_Reset` — no flow has been computed yet, which is more useful than reporting a rate no solver has taken.
  - An **inactive** edge reads `0.0`, not the rate it would carry if it were live.
  - **Signed.** `GSSK_GetEdgeQualityFlow` clamps a negative flow to zero following Odum's convention; this does not, because the sign of a flow is information a diagram consumer wants.

### Changed

- **The per-edge flow expression is written once.** The `switch` over edge logic existed twice — in `compute_derivatives` and, in a slightly different dialect, in `compute_quality_pass`. G4 needed it a third time, so it is now one static `edge_flow_rate(e, state, k)`. `k` is a parameter because the callers genuinely disagree about which `k` they mean: the flow cache passes `forced_edge_k`, the rate the derivative integrated; the quality pass passes `e->k`, which is what it has always used. That disagreement is very probably a defect — a forced edge makes emergy accounting disagree with the trajectory it is accounting for — but it moves published emergy numbers, so it is tracked separately as `quality-pass-ignores-edge-forcing` rather than fixed in passing. Behaviour is unchanged: every golden CSV and every focused suite passes untouched.

### Fixed

- **The npm package version was five majors behind the kernel.** `package.json` still said `1.0.0` while `include/gssk.h` said `GSK_VERSION_STRING "5.0.0"`; `scripts/release.sh` had only ever bumped the header. That is not cosmetic — the package ships `include/`, so a consumer who pinned `gssk@1.0.0` from npm was reading a header out of `node_modules` and getting the current one, or pinning a version that never described the kernel it was served. GIP-0001 was written this way: it quotes `INTERACTION /**< Multiplier flow (k * Q1 * Q2) */` and `GetStateSize /* Number of storage nodes */`, neither of which has been the comment for several majors.

  `package.json` is now `5.0.0`, `scripts/release.sh` bumps it alongside the header (and re-reads the file to confirm it is still valid JSON before committing), and `scripts/check_version_sync.py` fails the build when the two disagree — or when `GSK_VERSION_STRING` disagrees with the `GSK_VERSION_MAJOR`/`MINOR`/`PATCH` macros beside it. It is stdlib-only, runs as `make check-version`, and is a prerequisite of `make test`, so a skew cannot survive a local test run.

- **Out of bounds is now reportable.** `GSSK_GetNodeTypeString` returns `"storage"` for a NULL instance or an out-of-range index, which is indistinguishable from a genuine storage node and cannot be checked for. `GSSK_GetNodeType` returns `GSSK_NODE_INVALID` instead, following `GSSK_GetNodeID`'s contract. The string function's behaviour is unchanged — it is load-bearing elsewhere — but it now carries an `@warning` saying so, and a test asserts the warning is still true.

- **`GSSK_GetNodeTypeString`'s header comment listed eight of the nine primitives**, omitting `constant`.

### Documentation

- **The schema now says where limit logic's saturation constant C comes from.** `F = k × Q_origin / (1 + Q_origin/C)` has been implemented since the primitive was added and the formula is stated in `include/gssk.h`, but nothing told a consumer how to *supply* C. `gssk.schema.json` described `control_node` as a node that "modulates the flow" without saying it **is** the denominator constant, and described `threshold` for threshold and ratio logic only, never mentioning that it doubles as C when no control node exists. Raised as GIP-0001 G6.

  Both `EdgeParams` descriptions and the `EdgeLogic` enum now state the rule: `control_node`'s current Q supplies C if that node is named, otherwise `params.threshold` when it is above zero, and **`control_node` wins when both are given** — `threshold` is the source used in its absence, not a runtime fallback.

  Two behaviours that were not written down anywhere are now stated, in the schema and in the header:

  - A C taken from `control_node` is a **state variable, not a constant.** If that node is itself a store, the edge has a moving half-saturation point.
  - When C falls to `1e-9` or below the flow becomes **exactly 0.0 rather than an error.** A control node that decays toward zero therefore closes the pathway mid-run, in a model that loaded without complaint and whose file says nothing about it. That is deliberate — a saturation constant of zero means the pathway saturates at zero throughput — but it was invisible.

  The GIP filed the second point as "flow is silently 0.0 rather than an error". That is only half true, and the half matters: a limit edge with **neither** source of C is rejected at load with `Logic Error: Edge N (limit) requires control_node or threshold > 0`. Only a C that was valid at load and decayed afterwards reaches the solver.

  `threshold` logic's comparand is documented too: always `Q_origin`, never the control node, and the comparison is strict.

- **`GSSK_LIMIT_C_EPSILON`** is now a named public constant carrying that explanation, following the `GSSK_RATIO_EPSILON` precedent. The eight bare `1e-9` literals in the limit paths use it, so the documentation cannot drift from the threshold it documents. (The loop-limited node's `node_C` guard is a different constant with a different fallback and is untouched.)

- **`tests/test_limit_logic.c`** pins all four facts, so the documentation stays true. The precedence assertion is mutation-tested: swapping `control_node` and `threshold` priority in the kernel makes it fail. The decaying-control test asserts `Q_origin` is *exactly* frozen afterwards, and separately that the origin still holds most of its contents — otherwise a fully-drained store would sit still and pass.

- **A per-instance-price regression fixture for the full CLI path (GIP-0002).** The archetype `price_node` fix landed with a unit test; `examples/archetype_price_per_instance.json` is the artefact that stops it coming back through the route that test does not cover — `bin/gssk`, a golden trajectory, and the serialised-model fixture set.

  One `purchase` archetype holds a `constant` price member and an `exchange` naming it through `params.price_node`, instantiated as `a` and `b` and priced independently through `snapshot.state` on the expanded ids `a__price` (5.0) and `b__price` (12.0). The forward flow `k x Q_seller x Q_buyer` does not depend on price, which makes the fixture read as a controlled experiment: **the two inventories are bit-identical at every step while the two money stocks diverge**, ending at `a__spent` 293.79 against `b__spent` 705.10. `spent / inventory` is each instance's own price for the whole run, and `buyer + a__spent + b__spent` holds at 1000.0 because money is conserved across both diamonds.

  **Both regressions it guards against are poisoned so they cannot look plausible.** If the template `price_node` stops resolving, the exchange falls back to a scalar `price` of `999.0`: the buyer is drained to exactly 0.0 and both instances spend an identical 500.0. If instead a `snapshot.state` id stops resolving — which is *silent*, `GSSK_Init` skips an id it cannot find — both prices keep the template's declared `0.0`, no money moves at all, and the inventories run to 400.0 unpriced. Neither can be mistaken for the golden, and neither is a plausible-looking number.

### Fixed

- **`gssk.schema.json` no longer says a template's `price_node` is ignored.** It documented `price_node` as "accepted and ignored on a template", which stopped being true when expansion began rewriting it, and would have described the new fixture as relying on behaviour the schema disclaimed. The description now states what actually happens: a template's `price_node` names a sibling member by its template-local id, expansion rewrites it to `{instance}__{member}` exactly as it does an edge origin or target, and naming no member of the archetype is a schema error rather than a fallback to the scalar price.
- **Phase D.1 — money is a closed conserved loop, so `conserved: true` finally asserts something.** `examples/odum_gnp_loop.json` circulates money `households → exchange → firms → households` with no money source and no money sink, which is Odum's Fig. 3: energy flows *through* and leaves as heat, money goes *round*, counter-current the whole way. `make test-gnp-loop` pins that asymmetry — money's per-step `GSSK_GetCarrierConservationError` stays under the model's `solver_tolerance` (worst 2.2e-16) while energy's does not over the very same run (worst 1.5e-03, as the reserve depletes from 1001 to 42).

  **The declaration was previously vacuous, and the suite says so with a negative control.** Conservation is summed over *storage* nodes only, and C.4's `examples/odum_countercurrent.json` modelled money as a `source`/`sink` pair — no money storage at all. Its reported error is a flawless `0.0` computed over nothing, and `test_open_loop_conservation_is_vacuous` asserts both halves of that: zero money storages **and** a zero error. A check that cannot fail is worse than no check, so the new suite refuses to be one — it also asserts that the loop actually turns (`firms` peaks at 30% of the money stock, so the exchange really cleared) and that the topology contains no money source or sink by construction.

  Added to the PR quality gate and to `make test-linux`, not just the local chain.

## [5.0.0] - 2026-08-26

### Breaking

GSSK now **rejects models it previously accepted**. Both changes are listed in full below; they are collected here because they are the reason this is a major bump and not `4.2.0`.

- **An unrecognised model key is rejected** instead of silently ignored — `GSSK_Init` returns `GSSK_ERR_SCHEMA_VIOLATION` where it used to return `GSSK_SUCCESS`. Now enforced at every level of the model, not the five originally named.
- **An unrecognised node `type` is rejected** instead of silently becoming a `storage` node.

Neither accepts anything `gssk.schema.json` ever declared valid, so a model that validates against the published schema loads unchanged. What breaks is a model carrying a typo'd or extra key that the parser used to swallow — including, specifically, a node whose `type` was misspelled and which has therefore been silently simulated as a `storage` node, possibly for a long time. That is the case worth checking on upgrade: the error is the first time GSSK has ever told you about it.

**To upgrade:** run your models through `GSSK_Init` and read the messages, which name the offending key directly:

```
Schema Error: Node 'A' has unknown key 'bogus_key'.
Schema Error: Node 'A' has unknown type 'stroage'.
```

`make test-schema` validates a model corpus against the schema without running it.

### Added

- **`GSSK_EnsembleResult` no longer crosses the WASM boundary as a raw struct.** `GSSK_EnsembleForecast` returned a pointer and nothing could read it safely from JS, so every consumer decoded the fields by hand. `web/index.html` did it through `HEAPU32`, under a comment admitting it had gone and read `src/advanced.c` to find the `s * node_count + n` stride; a downstream user independently pinned the same offsets in a golden test, having concluded the layout was an undocumented ABI.

  Five flat accessors close it, following the `GSSK_GetCarrierID` / `GSSK_GetCarrierUnit` precedent that already exists for `GSSK_Carrier`: `GSSK_GetEnsembleNodeCount`, `GSSK_GetEnsembleStepCount`, `GSSK_GetEnsembleMin`, `GSSK_GetEnsembleMax` and `GSSK_GetEnsembleMean`. They take `(step, node)` and apply the stride internally, so it lives in the kernel once instead of in every caller. All five are exported to WASM and declared in `src/gssk.d.ts`; `web/index.html` now goes through them and its pointer arithmetic is gone.

  **Hand-decoding was never portable.** The fields sit at `0/4/8/12/16` under wasm32 but `0/8/16/24/32` in a native 64-bit build, where `sizeof(GSSK_EnsembleResult)` is 40 — any pinned-offset reader is wrong on one of the two targets, and `-sMEMORY64` moves them again.

  **The semantics were documented in `include/gssk.h` the whole time**, and are now harder to miss. The three arrays are pointwise statistics *across runs* — min, max and mean — not three sampled trajectories, so `min <= mean <= max` holds everywhere by construction. The downstream report that two of the arrays "exchange places from sample to sample" was ties, not disorder: a constant node ties at every step, and step 0 of every node ties because perturbation only touches edge `k`. On `examples/economic_model.json` that is 3,006 of 6,006 points. `test_ensemble_getters` in `tests/test_advanced.c` asserts the ordering pointwise, asserts the getters agree with the fields they expose, pins both tie classes, and checks that a non-constant node genuinely spreads so the ordering assertions cannot pass on a degenerate ensemble.

- **Phase E.1 — the countercurrent is a first-class example, with an annotated twin.** C.4 shipped `examples/odum_countercurrent.json` carrying its own commentary inside `_`-prefixed fields. That is the wrong file for it: every other reference model here is a plain model beside a separate annotated twin (`household_model.json` / `household_model_annotated.json`), and someone reading for the topology should not have to read four essays to find thirteen nodes. The prose now lives in `examples/odum_countercurrent_annotated.json`, and the plain file is the model alone.

  **The golden CSV did not change by a single byte.** That is the evidence that this was a re-shelving and not an edit — `tests/expected/odum_countercurrent.csv` is bit-for-bit what C.4 generated, and `make test-net-energy` still asserts the same claim against the stripped file.

  `make test` now also checks **annotated twins against the models they document**. Each `X_annotated.json` gets its own golden CSV, byte-identical to `X.csv`, and the harness compares the two trajectories directly. Golden CSVs alone cannot catch this class of drift: edit the annotated model, run `make test-update`, and its golden is regenerated to match the model that has just diverged, so everything passes while the documented model and the running model are different systems. Verified by negative control — perturbing the exchange `k` in the annotated file to `0.0031` and regenerating its golden passes every per-model comparison and fails the twin check.

  The fuzz corpus gains `tests/fuzz_corpus/seed_closed_money_loop.json`: money circulates `buyer → exchange → till → buyer` with no money source and no money sink, so the parser is seeded with the closed-loop topology `d1` will build on, and money is conserved over the run. `seed_exchange_node.json` stays, as the open-loop counterpart.

- **Phase C.4 — inflation emerges from net-energy decline.** `examples/odum_countercurrent.json` builds Odum's Figure-2 structure: the fuel reserve is a depleting **storage**, extraction is a work gate needing both the reserve and the machinery, and the energy spent *getting* energy is fed back as a cost inversely proportional to what is left. As the reserve depletes the feedback fraction rises, the net energy reaching the economy collapses faster than the gross, and price — money over real throughput — rises.

  Measured over the run: fuel `1000 → 20`, structure booms `1 → 391` by `t≈29` then busts to `22`, net-energy-per-gross falls `0.85 → 0.27`, and the price index rises **`25 → 18,708`, a 739× inflation**, monotonically across the whole depletion phase.

  **The money supply is exactly constant for the entire run** — `buyer` is a `source`, so `compute_derivatives` pins it — and `make test-net-energy` asserts that at every step. That is the control that makes the claim mean something: with `M` fixed, the rise in `P = M/W` is attributable to `W` alone and cannot be a monetary effect. This is Odum (1973) points 1–3, and it is the piece that makes inflation emerge from net-energy decline rather than from supply and demand.

  **No kernel change was required, which was not expected.** The task was filed against `src/gssk.c`. What unblocked it was C.3: `params.numerator_node` ([ADR 0005](adr/0005-price-relaxation-and-named-numerator.md)) lets a `ratio` edge name its numerator instead of taking it from the origin, so the gross-energy tap can read `k·Q_yield/Q_fuel` from a **pinned** origin without consuming either node. Without it that tap would have had to originate at `yield` and would have double-debited the very flow it exists to measure — the hazard ADR 0003 named and ADR 0005 closed.

  `make test-net-energy` asserts the **claim**, not a trajectory — the golden CSV in `make test` already pins the digits and can say nothing about whether the mechanism is the one described. It checks boom/bust, one-way depletion, the feedback fraction rising against its closed form, monotone net-per-gross and price across the depletion phase, the pinned money supply, and that price *tracks* `M/W` (worst departure 4.3%, the relaxation lag) rather than merely correlating with it. Both halves are verified by negative control: breaking the `α` equality produces a 51% departure, and removing the depletion feedback holds net/gross at exactly `1.0000`.

- **Forcing functions — Odum's eleven, as one waveform vocabulary attached in two places.** A source node held its declared value for the whole run, so GSSK expressed exactly two of the eleven forcing functions in *Systems Ecology* Fig. 7-2 (constant force and constant flow) and had no representation for the other nine.

  `forcing` on a **node** drives its held value (Odum X/N, a *force*); `forcing` on an **edge** drives its rate `k` (Odum J, a *flow*). Eight waveforms: `step`, `impulse`, `ramp`, `sawtooth`, `square`, `sine`, `exponential`, `jitter`. Odum's eleven are a node-value-versus-edge-rate distinction crossed with a carrier distinction, and carriers were already modelled — so eleven node types would have been the wrong shape. [ADR 0006](adr/0006-forcing-one-vocabulary-two-attachments.md) records that decision and what was rejected. Worked model in `examples/forced_source_model.json`; full vocabulary and formulas in [docs/concepts.md](concepts.md#forcing-functions).

  **Waveforms are evaluated at solver STAGE times, not once per step.** This is the requirement that is easy to get wrong and invisible when you do: sampling once per step leaves the forcing first-order while the state is fourth- or fifth-order, and the run still completes and still looks smooth. Integrating a sine-forced source against its closed form, the error ratio per halving of `dt` is **16.02× / 16.01× / 16.00×** with stage times and **1.97× / 1.98× / 1.99×** without. Every other test in the suite passes either way. `h8a-thread-time-through-ode-core` landed first to make the stage times available.

  **`jitter` is latched once per accepted step**, drawn from the instance SplitMix64 stream (`GSSK_SetSeed`), never libc `rand()`. A per-stage draw would make the trajectory depend on solver internals — the same model answering differently under `rk4` and `dopri5` for reasons that are not physics. The test asserts the *draw sequence* is bit-identical across `rk4`, `incipient` and `adaptive`.

  **A storage node cannot be forced** — its value is the integral of its flows, so forcing it asserts two things about one quantity. `GSSK_Init` and `GSSK_AddNode` both reject it naming the node, rather than ignoring the block, which is the failure mode `h8b` was landed to remove. A periodic waveform without a positive `period` is likewise rejected rather than silently treated as a constant.

  **`phase` is a time offset, not an angle** — same units as `t`, and *subtracted*, so a positive phase delays the waveform. **`impulse` is area-normalised** over one nominal `dt`, so its integral is `area` at any step size. Both conventions are stated in the header, the schema, `docs/concepts.md` and the example, because an ambiguous convention is how two implementations diverge while both look right.

  New accessors, all exported to WASM and typed in `src/gssk.d.ts`: `GSSK_GetNodeForcingKind`, `GSSK_GetEdgeForcingKind`, `GSSK_EvaluateNodeForcing`, `GSSK_EvaluateEdgeForcing`. They are **the same evaluator the derivative path uses** — a test drives each of the eight waveforms and asserts what the evaluator reports equals what the kernel integrated, so the two cannot drift. Flat scalars rather than a struct pointer, for the reason the flat carrier getters exist. `js/gssk.js` gains the matching wrappers.

  Forced node values are written into the live state after each step, so `GSSK_GetState` and the CSV show the waveform. Previously the derivative was correct while a sine-forced source appeared as a flat line next to the storage it was visibly driving.

### Changed

- **`GSSK_Reset` does not rewind the random stream, and now says so.** Behaviour is unchanged; the contract is newly documented because forcing makes it reachable. Rewinding was tried and is wrong: `GSSK_EnsembleForecast` and `GSSK_CalibrateMonteCarlo` perturb with the instance RNG and then reset once per run, so rewinding collapses an ensemble to one trajectory (`test_advanced`'s calibration caught it). `GSSK_Reset` means "back to `t_start`", not "back to the start of the stream". To repeat a `jitter` run exactly, call `GSSK_SetSeed(inst, GSSK_GetSeed(inst))` first.

### Known limitation

- **Trajectories are not bit-identical across platforms**, because the kernel uses the platform's `libm`. Measured three ways on the same eight sample points: the WASM build agrees with Linux GCC **exactly**, and macOS Apple clang's `sin()` is the outlier, differing by 1 ULP at two of the eight. `sin`/`exp`/`pow` are not required by IEEE-754 to be correctly rounded, so this is inherent rather than a defect in any toolchain. It pre-dates forcing — `exp()` was already on the Riccati duet path and `pow()` on the adaptive step controller — but forcing makes it common, since a sine-forced model hits a transcendental on every stage of every step. `make test-forcing-wasm` therefore asserts agreement to 4 ULP rather than bit-equality, with the measurement recorded in the check itself. Whether GSSK should ship its own correctly-rounded transcendentals is a Phase G reconstruction question, tracked as `deterministic-transcendentals-cross-platform`.

### Changed

- **Simulation time is now threaded through the ODE core.** `compute_derivatives` takes an explicit `double t`, and every call site passes the correct *stage* time: classical RK4 at `t`, `t+h/2`, `t+h/2`, `t+h`; DOPRI5 at `t + c_i·h` for `c = (0, 1/5, 3/10, 4/5, 8/9, 1, 1)`. The whole derivative surface is covered, not just RK4 — both `rk4_step` variants, DOPRI5's seven stages, the IDC path (`build_flow_matrix`, `build_forcing_vector` and the processing-node helpers), `build_jacobian`, `compute_param_deriv`, `compute_quality_pass`, `compute_quality_sensitivity`, threshold sub-stepping, and the adjoint's backward integration.

  **Nothing consumes `t` yet, and every trajectory is bit-identical.** `tests/expected/` is untouched and `make test` passes against it unchanged — that is the acceptance criterion, not a side note. Landing the threading on its own means any later trajectory change is attributable to the feature that uses `t` rather than to a mistake in threading it through 16 call sites and seven solver stages.

  The DOPRI5 c-nodes were *already written down* — `/* Stage 2 at c2 = 1/5 */` and so on — and thrown away, because there was no `t` for them to offset. They are now used rather than described.

  This exists for forcing functions, and it is the trap that feature would otherwise fall into: a waveform sampled once per **step** instead of once per **stage** is first-order while the state is fourth- or fifth-order. The run completes, the trajectory looks smooth, and the order loss is invisible without a convergence study.

  `make test-stage-times` records the time the solver actually hands each derivative evaluation and asserts the sequence, rather than re-deriving the arithmetic in the test — which would just be the same mistake written twice. It covers the two sharp cases: adaptive sub-stepping, where stage 1 of sub-step 2 must be at `t + h₁` and not `t`, and the adjoint, which runs time backwards and must land on `t_start`. The recorder is compiled in only under `-DGSSK_STAGE_TIME_PROBE`, which only that target defines; it is not in `gssk.h`, not in the shipped library, and not in the WASM export list. **No public API change and no schema change.**

### Added

- **Flat carrier accessors that carry no struct layout across the WASM boundary.** `GSSK_GetCarrierID`, `GSSK_GetCarrierUnit`, `GSSK_GetCarrierConserved` and `GSSK_FindCarrierIdx`, all four exported to WASM and typed in `src/gssk.d.ts`.

  The data was always reachable from C — `GSSK_GetCarrier` returns the whole `GSSK_Carrier`. The gap was JS: across WASM that same call is a bare heap pointer, so reading `unit` or `conserved` meant assuming field offsets, the width of `bool` and the absence of trailing padding, none of which is an ABI contract and all of which break by returning plausible garbage rather than by failing. `unit` is also the y-axis label a plotting consumer would otherwise hardcode, and `unit` plus `conserved` is what a consumer needs to decide two series may not share a scale (ADR-6, ADR-8).

  `GSSK_GetCarrier` is unchanged and stays for C consumers. The out-of-range conventions differ deliberately and are documented on both: the string getters return `""` and never `NULL`, following `GSSK_GetNodeCarrier`, while `GSSK_GetCarrier` still returns `NULL`; `GSSK_GetCarrierConserved` returns `0`, which is indistinguishable from a declared non-conserved carrier, so bound-check against `GSSK_GetCarrierCount` first; `GSSK_FindCarrierIdx` returns `-1`, matching `GSSK_FindNodeIdx` / `GSSK_FindEdgeIdx`. `make test-carrier-api` covers each getter, every out-of-range path, the lookup round-trip, and — the check that matters — that the flat path and `GSSK_GetCarrier` agree at every index, so the two cannot drift.

- **Phase C.3 — price is now a state variable that relaxes toward Odum's ratio.** `dP/dt = α(M/W − P)`, with `α` a settable per-model adjustment rate rather than a hard-coded one. `examples/price_dynamics_model.json` wires circulating money `M` and the C.2 delivered-work signal `W` into the C.1 `ratio` primitive and feeds the result back into the transaction diamond through the C.0 `price_node` hook. The fixed point is `M/W` itself, not something proportional to it: `make test-price-dynamics` asserts convergence to `M/W` to 1e-9 and the approach to `(M/W)(1 − e^{−αt})`, so both the level and the time constant are hand-checked rather than golden. `examples/emergent_price_model.json` is unchanged and stays in the suite, so the difference between the Tier 1 proportional anchor and the Tier 2 true ratio remains visible.

  **What this is, precisely:** the *imposed* transactor with a dynamically computed price. Odum draws two price notations — one where the transaction generates price (*Systems Ecology* Fig. 23-3c) and one where price is determined outside and imposed on the pathway (Fig. 23-3d) — and GSSK's exchange node is the second: `F_money = P × F_goods`. His endogenous price is a ratio of two independently driven *flows*, `p₁ = J₁/J₄` (Fig. 23-2c), which cannot be reached by measuring this diamond's own flows, since `J_money := P·J_goods` makes the quotient an identity. `M/W` is dimensionally a price (AUD/kg) and proportional to `J₁/J₄`, with the constant set by β. So price here is endogenous in where the number comes from and exogenous in how it acts on the trade — a real and common configuration, and the one the task specified, but not the price-generating transactor. [ADR 0005](adr/0005-price-relaxation-and-named-numerator.md) records both the relaxation-over-algebraic decision and this fidelity boundary in full.

- **`ratio` edges accept `params.numerator_node`.** The numerator is now an optional *named* operand, read by id and not consumed — the same contract `control_node` has had as the denominator. Omitted, the numerator is `Q_origin` and behaviour is bit-for-bit unchanged.

  This exists because an edge debits its origin, so before it the only way to put a stock in the numerator was to drain it: a price mechanism reading `M` would have eaten the money supply it was observing. ADR 0002 called for a ratio whose numerator and denominator are "both named and distinguishable" and then named only the denominator; this supplies the other name. The edge remains a flow from origin to target — pin that origin with a `source` or `constant` node when the flow must cost nothing, as [ADR 0003](adr/0003-delivered-work-signal.md) established for the `W` tap.

  An unknown `numerator_node` id is a linkage error, and `numerator_node` on any logic other than `ratio` is a logic error rather than a silently ignored key.

### Fixed

- **Deactivation did not survive serialise → reload.** `GSSK_DeactivateEdge` clears `active` *and* sets `k` to `0`, and `GSSK_Init` accepted the emitted `"active": false` without acting on it — so a round-trip produced an edge that was **active with `k = 0`**, not an inactive edge. The trajectory matched either way, which is why it survived: `k = 0` kills the flow regardless, so the only thing that differed was the flag.

  The flag is not decorative. It is read at roughly twenty sites in `src/gssk.c`, and the ones that matter **count elements rather than sum flows** — which is exactly what `k = 0` cannot stand in for. `network_is_isolated_duet` requires exactly one active edge before the Riccati closed form is used, motif detection skips inactive nodes and edges, and the closed-system conservation check sums only active nodes. A reloaded model could therefore be treated differently from the model it was serialised from while producing identical numbers.

  **The node half was worse and was lost outright.** `build_topology_json` emitted no `active` for nodes at all, and a node has no `k` to carry the deactivation the way an edge accidentally did, so `GSSK_DeactivateNode` did not survive in any form. `Node.active` is now emitted (only when false), declared in `gssk.schema.json` and listed in `NODE_KEYS`.

  **The flag is never inferred from `k`.** An edge authored with `k: 0` is present and carrying nothing; a deactivated edge has been taken out of the network. Reading the flag back from the conductance would reproduce the trajectory and lose the topology — the same defect one level down. `tests/schema_fixtures/deactivated_elements.json` pins both directions, and puts an `active` key into `tests/results/serialized/` for the first time: that corpus never contained one, which is why the original defect went unseen.

  `make test-deactivation` (`tests/test_deactivation_round_trip.c`) asserts on motif count and on what happens when `k` is restored — an inactive edge stays dead, an authored-zero edge starts flowing — because a trajectory assertion cannot distinguish the two by construction. It fails eight assertions against the pre-fix kernel. `GSSK_AddNode` and `GSSK_AddEdge` honour `active` too, being separate parsers. [ADR 0004](adr/0004-schema-advisory.md) is amended with the closure, including a correction: the flag is *not* read by `GSSK_ReclassifyNetwork`, which inspects no edges at all.

- **A `ratio` edge serialised as a `linear` edge.** `logic_type_str` was never given a `ratio` case and fell through to its `"linear"` default, so every round-trip through `GSSK_SerializeModel` or `GSSK_SerializeSnapshot` silently replaced the division with a proportional flow — permanently, and with no error. This reached snapshots and the Phase G archival dumps, where the serialised form *is* the artefact. ADR 0002 tabulated the eight sites a new logic type must be added to; the serialiser was not among them, which is how it was missed.

- **A `ratio` edge's denominator floor did not survive serialisation.** `params.threshold` is the floor override for `ratio` logic, but it was emitted only for `threshold` logic. A deliberate floor of `0.01` round-tripped back to `GSSK_RATIO_EPSILON` (1e-9) — a 1e7× change in the price a model reports once `W` empties.

### BREAKING

- **Unknown-key rejection now covers the whole model, not the five levels `h8b` named.** `GSSK_Init` returns `GSSK_ERR_SCHEMA_VIOLATION` for an unrecognised key in node `params`, `metadata`, a `carriers[]` entry, `snapshot` (and its nested `state[]`, `edge_k[]`, `solver`, `rng_state` and `mutation_log[]` objects), an `archetypes` template (and its node, edge and `params` objects), and the root-level `mutation_log[]`. Messages name the key **and** its container, as before — `Schema Error: Node 'mkt' params has unknown key 'pric'.`, `Schema Error: metadata has unknown key 'authorr'.`, `Schema Error: Carrier 'energy' has unknown key 'unti'.`, `Schema Error: Archetype 'widget' edge 'bleed' has unknown key 'origen'.`

  **Node `params` is the sharpest of these** and is the direct analogue of edge `params`, which `h8b` did close. `{"type":"exchange","params":{"pric":10}}` used to load, ignore `pric`, and run the transaction at the default price — the same "plausible completed run, quietly different model" failure `h8b` exists to prevent, on the parameter a Phase C price model turns on. `GSSK_AddNode` applies it too, being a separate parser; `node_keys_ok` is shared between the two.

  **`metadata` matters for a different reason**: it carries `model_hash`, which the kernel round-trips and never computes. A typo'd provenance key was silently dropped from an artefact whose whole purpose is provenance.

  **`^_` keys are accepted at every new level**, as everywhere else. Every one of these levels was **already** declared with `additionalProperties: false` and `^_` `patternProperties` in `gssk.schema.json`, so this is the same argument as `h8b` — bringing the kernel into agreement with a contract the project already publishes — not a new rule.

  **Two things the sets record rather than paper over.** `snapshot.dt` is **emitted by `GSSK_SerializeSnapshot` and never read back**: reload takes `dt` from `config`, as it does for a model with no snapshot at all. Deriving `SNAPSHOT_KEYS` from the parser alone would have omitted it and rejected every snapshot the kernel has ever written — the drift direction that breaks *working* models, which is worse than the bug being fixed — so it is in the set with a comment saying why, and has its own assertion. Separately, `archetypes` node templates are declared with the full `NodeParams` `$def` while `parse_user_archetypes` reads only `k`, `C`, `threshold` and `price`; the wider published set is honoured, because narrowing the kernel to the parser would reject models the schema calls valid.

  **Migration**: identical to `h8b`'s. A model carrying a stray key at any of these levels now fails to load; correct the key or `_`-prefix it. Nothing in `examples/`, `tests/schema_fixtures/` or `tests/results/serialized/` changed — all 21 loadable models still load and all 42 serialised artefacts still validate. `make test-unknown-keys` covers each new level with one rejection naming key and container, one `_`-prefixed acceptance, and the corpus regression; it fails 30 assertions against the previous kernel. [ADR 0004](adr/0004-schema-advisory.md) is amended with the closure.

- **An unrecognised model *key* is now rejected instead of silently ignored.** `GSSK_Init` returns `GSSK_ERR_SCHEMA_VIOLATION` for a key it does not recognise at the root object, in a node object, in an edge object, in edge `params`, or in `config`. The message names both the key and its container so an authoring UI can highlight the element that is wrong — `Schema Error: Edge 'e1' params has unknown key 'bogus_param'.`

  The hazard is sharpest for a feature the kernel does not yet have. This model used to load, return `GSSK_SUCCESS` and run to completion with `forcing`, `nonsense_top_level_key` and `bogus_param` all ignored in silence:

  ```json
  {"metadata":{"schema_version":4},
   "forcing":{"sun":{"waveform":"sine","amplitude":5,"period":24}},
   "nonsense_top_level_key":123,
   "nodes":[...], "edges":[{"...":"...","params":{"k":0.5,"bogus_param":9}}]}
  ```

  A model authored against a kernel that *has* forcing, loaded by one that does not, therefore produced a plausible completed run with no diagnostic. Its JSON — and any external content hash of it, which is what `metadata.model_hash` carries, since the kernel round-trips that field and never computes it — says "forced". Its trajectory says "constant". Nothing reconciles the two. This is the same hazard class as the node-`type` fallback below: a wrong key is not a smaller mistake than a wrong type, it is the same mistake one level up.

  **Any key beginning with `_` is accepted everywhere, at every level.** That is load-bearing rather than a courtesy: `examples/household_model_annotated.json` and `examples/price_dynamics_model.json` carry `_note` and `_mechanism` throughout, and `gssk.schema.json` documents the convention.

  `GSSK_AddNode` and `GSSK_AddEdge` apply the same check — they are separate parsers. A rejected add is a true no-op: the check runs before the `realloc`, so counts are unchanged and the instance is still steppable.

  **This makes the kernel agree with a contract the project already publishes.** `gssk.schema.json` has always set `additionalProperties: false` at the root and on `Node`, `Edge`, `EdgeParams` and `Config`, and permitted `^_` keys via `patternProperties`. `make test-schema` caught violations in `examples/`; nothing caught them at runtime for a *consumer's* model, which is where it matters. See [ADR 0004](adr/0004-schema-advisory.md).

  **Migration**: a model carrying a stray key now fails to load; correct or `_`-prefix the key. Nothing in `examples/`, `tests/schema_fixtures/` or the fuzz corpus changed — all 18 loadable models still load, and the two fuzz seeds that fail still fail for their original, unrelated reasons (one is not JSON, one has no `nodes` array).

### Fixed

- **The serialiser emitted `"active"`, a key the published schema forbids.** `build_topology_json` writes `"active": false` for an edge deactivated via `GSSK_DeactivateEdge`, but `active` appeared nowhere in the parser and `Edge` sets `additionalProperties: false` — so the kernel was emitting output its own schema rejects. It went unnoticed because no model in `examples/` or `tests/schema_fixtures/` has a deactivated edge, so the `tests/results/serialized/` corpus never contained one. Found by the stricter parser above, which is exactly what it is for. `active` is now declared in the schema and accepted on load, and `tests/test_unknown_keys.c` covers the deactivated-edge round-trip directly.

  Note what is *not* fixed here: `GSSK_Init` accepts `active` but does not act on it. Deactivation survives the round-trip through `params.k`, which `GSSK_DeactivateEdge` sets to `0.0`, so the trajectory is reproduced — but `edges[i].active` is not, and that flag is read by topology classification, so a reloaded edge is *active with k = 0* rather than *inactive*. That is a behavioural change and is tracked separately.

### BREAKING

- **An unrecognised node `type` is now rejected instead of silently becoming a `storage` node.** `parse_node_type` returned `NODE_STORAGE` for any string it did not recognise, so `"storge"`, `"Source"` or `"producer_"` produced a *different model* that ran to completion and reported success. `GSSK_Init` now returns `GSSK_ERR_SCHEMA_VIOLATION` for any type that is neither one of the nine primitives (`storage`, `source`, `sink`, `constant`, `interaction`, `gain`, `loop_limited`, `exchange`, `switch`), a built-in composite (`producer`, `consumer`, `misc_box`, `system_frame`), nor an archetype declared in the model's own `archetypes` block. The message names the node id and the offending string — `Schema Error: Node 'grasss' has unknown type 'storge'.` — so an authoring UI can highlight the element that is wrong.

  `GSSK_AddNode` rejects the same strings, and additionally rejects composite and archetype names: it performs no expansion, so `{"type":"producer"}` added at runtime had been becoming a single storage node rather than the producer subgraph. It now fails with a message saying composites can only be added at `GSSK_Init`. A rejected add is a true no-op — nothing is allocated or grown before the check, so the instance a drag-and-drop editor is mutating is left exactly as it was and remains steppable. Expanding composites at runtime is a separate change and is not attempted here.

  **Migration**: a model relying on the old fallback now fails to load. The fix is to correct the type string; every previously-accepted string that was actually a primitive, a built-in composite, or a declared archetype is unaffected. Nothing in `examples/`, `tests/schema_fixtures/` or the fuzz corpus changed.

  This closes the hazard [ADR 0004](adr/0004-schema-advisory.md) left open, using the fix that ADR named. Schema validation could never have caught it: `Node.type` cannot be a closed enum, because archetype names are user-defined, so a validator cannot tell a typo from a legitimate archetype reference. The parser can — it reads the `archetypes` block before the node list — which is why the check belongs there and why the schema can stay advisory.

### Fixed

- **`docs/api-reference.md` documented `GSSK_Carrier.id` as `id[64]`.** It is and always has been `char id[32]`. Harmless in C, but the API reference is exactly where a JS consumer would go to work out the struct offsets to decode by hand — a reader who trusted it would have read `unit` 32 bytes past where it lives. Noticed while writing up the flat getters; corrected, and the reference now tells JS not to decode the struct at all.

- **`GSSK_AddNode` read freed memory when reporting a duplicate id.** The error message formatted `id->valuestring` after `cJSON_Delete` had already released the tree that owned it. Found while adding the type check next to it; the message is now formatted before the delete.

- **GCC 11 can build the kernel again.** `append_mutation` used `strncpy` followed by an explicit terminator, a pattern GCC's `-Wstringop-truncation` rejects under `-Werror` even though it is correct. Every `strncpy` in `src/gssk.c` (107 sites) is now `safe_str_copy`, which always NUL-terminates and derives its bound from `sizeof(dst)` so the bound cannot drift from the field width. This removes the limitation recorded against 4.1.0.

- **Schema conformance is now tested.** `make test-schema` (and `make test`) validates three corpora against `gssk.schema.json` and fails the build on a mismatch: the hand-written models in `examples/`, the corner-case fixtures in `tests/schema_fixtures/`, and — via the new `bin/dump_serialized` — the JSON that `GSSK_SerializeModel` and `GSSK_SerializeSnapshot` actually emit for every one of them. The schema and the parser can no longer drift apart unnoticed in either direction. CI installs `jsonschema` to make the gate real; locally it skips with a message when the dependency is absent. See [ADR 0004](adr/0004-schema-advisory.md) for why the schema stays advisory rather than being enforced inside `GSSK_Init`.

- **Fixed: the schema rejected the kernel's own adaptive-solver config.** `config.rel_tol`, `abs_tol`, `h_min` and `h_max` are read by `GSSK_Init` and written back by `GSSK_SerializeModel`, but `Config` did not list them and set `additionalProperties: false` — so any model using DOPRI5 tolerances, including one the kernel had just serialised, failed validation against the project's own schema. All four are now described. The root-level `mutation_log` block is also documented for what it is: an archival copy, which `GSSK_Init` restores only from `snapshot.mutation_log`.

---

## [4.1.0] — 2026-08-17

The first release published under an immutable version tag. Everything below shipped since 3.6.0; the intervening 4.0.0 was distributed only through the rolling `latest` pre-release and was never tagged or recorded here. Consumers who pinned artifacts from that rolling tag should move to `v4.1.0`, which differs from it — this release adds nine exported functions and changes how the stochastic entry points draw randomness.

### Added

- **Phase 7 — complete ESL node type taxonomy.** All seven Odum symbols are now implemented: `interaction`, `gain`, `loop_limited`, `exchange` and `switch` join `source`, `storage`, `sink` and `constant`. Processing nodes are configured through the node's `params` block (`k`, `C`, `threshold`, `price`) rather than through edge parameters, consistent with the ESL topology rule. Schema v4 with a `gssk migrate --from 3` path.
- **Phase 8 — composite node types and archetypes.** Built-in `producer`, `consumer`, `misc_box` and `system_frame` composites expand at `GSSK_Init` time into namespaced primitives (`{instance}__{member}`). User-defined templates may be registered via a top-level `archetypes` block and used as node types. Ports define external attachment.
- **Phase 9 — runtime pattern discovery.** Recurring 2–3 node subgraph motifs are detected after each step and, once stable, promoted to named archetypes via `GSSK_ProposeArchetype`. `GSSK_GetGenerativityIndex` reports the rate of emergence.
- **Composite membership API** (GH #29 item 1): `GSSK_GetNodeComposite`, `GSSK_GetNodeRole`, `GSSK_GetCompositeMemberCount`, `GSSK_GetCompositeMemberIndex`, `GSSK_GetCompositeArchetype`. Membership is recorded during expansion, so consumers no longer have to infer it by string-matching the `{instance}__{member}` prefix — an inference that is unsound in both directions and silently corrupts aggregation.
- **Seedable randomness** (GH #29 item 5): `GSSK_SetSeed`, `GSSK_GetSeed`, `GSSK_NextRandom`, `GSSK_NextRandomUniform` and `GSSK_DEFAULT_SEED`. Snapshots now carry `{seed, state}` in place of the previous null placeholder.
- **Containerised Linux toolchains**: `make wasm-container`, `make test-linux`, `make test-linux-clang`, `make ci-local` build WASM and run the suite under real GCC from macOS.
- **Whitepaper and article** under `doco/`, built with `make doco`.

### Changed

- **Stochastic entry points are now reproducible.** `GSSK_EnsembleForecast` and `GSSK_CalibrateMonteCarlo` draw from an instance-owned SplitMix64 generator seeded at init instead of libc `rand()`. Same model plus same seed now gives bit-identical results across platforms and under WASM. **This is a behavioural change**: `srand()` no longer influences either function, so callers that relied on it must use `GSSK_SetSeed`.
- **`gssk.schema.json` regenerated for v4** (GH #29 item 2). The published schema rejected models the kernel accepts. Beyond the missing `archetypes` block and node types, it also required `logic` and `params.k` on every edge, had no node `params` block, declared `snapshot` as a closed empty object so no serialised snapshot could validate, and rejected the `_`-prefixed annotation convention used in the bundled examples. `Node.type` is now an open string, since user archetype names are open-ended and unrecognised types are not rejected by the kernel.
- **`docs/concepts.md` corrected** (GH #29 item 3). It described composites as future work although they shipped, and the composite table was wrong on the facts. The same stale future tense applied to Phase 7 and Phase 9.
- **Releases now publish immutable version tags** (GH #29 item 6) with `gssk.schema.json` attached alongside `gssk.js`, `gssk.wasm` and `gssk.d.ts`, and a SHA-256 table in the release notes. The rolling `latest` pre-release continues, now labelled as republished in place.
- **CI builds WASM on pull requests.** Previously only the deploy job built WASM, so export changes were unverified until after merge.

### Fixed

- `GSSK_AddNode` zeroed the node struct, which would have made every runtime-added node report membership in composite 0.
- Monte Carlo calibration selected population indices with `rand() % n`, biasing toward low indices; now rejection sampling.
- `src/gssk.d.ts` reported the schema version range as "2 or 3" (GH #29 item 4).
- The release workflow force-pushed the bare version tag onto an orphan dist commit lacking `Package.swift` and `src/`, which would have broken SPM resolution and made every version tag mutable. The dist tree now uses a `dist-vX.Y.Z` namespace.

### Known limitations

- An unrecognised node `type` is not rejected — the kernel falls back to `storage`, so a typo yields a silently incorrect model. Validate against `gssk.schema.json` before calling `GSSK_Init`.
- `system_frame` is structural only: it reserves a name but expands to no subgraph. The ESL switching-box composite is unimplemented.
- Structural capacities are fixed at compile time: 32 archetypes, 128 composite instances, 16 nodes and 32 edges per archetype.
- Motif detection is skipped above 64 nodes.
- Building with GCC 11 fails on a `stringop-truncation` diagnostic in `append_mutation`. GCC 13, Clang and Apple Clang are unaffected.

---

## [3.6.0] — Phase 6: Productionisation

### Added
- **CI/CD**: GitHub Actions matrix (`ci.yml`) covering gcc/clang × linux/macos, AddressSanitizer + UBSan, lcov coverage gate (≥ 85%), Valgrind, LibFuzzer (30 s), perf regression, Swift 5.10 + 6.0
- **Python binding** (`python/gssk.py`): full ctypes wrapper with `GSSKSimulator`, `from_file`, sensitivity, carriers, serialisation, optional pandas `run_dataframe()`
- **Python tests** (`python/test_gssk.py`): 31 tests via `make test-python`
- **JavaScript/TypeScript wrapper** (`js/gssk.js`): `GSSKSimulator` ES module with async `create()`, carrier access, mutation log
- **Benchmark suite** (`bench/`): `run_bench.sh` + models at 10/100/1 000 nodes; `make bench` and `make bench-check` targets
- **VitePress docs restructure**: `docs/concepts.md`, `docs/api-reference.md`, `docs/cookbook.md`, `docs/CHANGELOG.md`; sidebar reorganised into Guides / Reference / Examples sections
- **Release script** (`scripts/release.sh`): bumps `GSSK_VERSION` in `include/gssk.h`, tags git
- **SECURITY.md**: vulnerability reporting process
- **LICENSE**: MIT

### Changed
- `Makefile`: added `shared`, `test-python`, `asan`, `test-asan`, `coverage-build`, `coverage-report`, `coverage-check`, `test-valgrind`, `fuzz-build`, `fuzz-run`, `bench`, `bench-check` targets
- WASM exports list updated to include carrier functions

---

## [3.5.0] — Phase 5: Multi-Carrier Networks

### Added
- `carriers` array in schema (v3): `{id, unit, conserved}`
- `GSSK_GetCarrierCount`, `GSSK_GetCarrier`, `GSSK_GetNodeCarrier`, `GSSK_GetEdgeCarrier`, `GSSK_GetCarrierConservationError`
- Swift binding: `GSSKCarrier` struct, `carrierCount`, `carrier(at:)`, `carriers`, `nodeCarrier(at:)`, `edgeCarrier(at:)`, `carrierConservationError(for:)`
- TypeScript declarations for carrier functions in `src/gssk.d.ts`
- Household model (`examples/household_model.json`): 4 carriers × 23 nodes
- Household model documentation (`docs/examples/household/README.md`)
- Annotated household model (`examples/household_model_annotated.json`)
- Jupyter notebook (`examples/household_notebook.ipynb`)
- Browser interactive demo (`docs/examples/household/demo.html`)

### Fixed
- LIMIT edge validation: `threshold > 0` is accepted as the saturation constant without requiring `control_node`

---

## [3.4.0] — Phase 4: Advanced Analysis

### Added
- Forward (tangent-linear) sensitivity: `GSSK_EnableForwardSensitivity`, `GSSK_DisableForwardSensitivity`, `GSSK_GetSensitivity`
- Adjoint sensitivity: `GSSK_RunAdjoint`, `GSSK_GetTransformitySensitivity`
- Gradient calibration: `GSSK_CalibrateGradient`
- Monte Carlo calibration: `GSSK_CalibrateMonteCarlo`
- Ensemble forecast: `GSSK_EnsembleForecast`, `GSSK_FreeEnsembleResult`
- Event detection: `GSSK_GetEventCount`, `GSSK_GetEventTime`, `GSSK_GetEventEdgeID`, `GSSK_GetEventDirection`

---

## [3.3.0] — Phase 3: Structural Mutation

### Added
- Mutation log with `GSSK_GetMutationCount`, `GSSK_ExportMutationLog`, `GSSK_ClearMutationLog`, `GSSK_SetMutationCause`
- Dynamic topology: `GSSK_AddNode`, `GSSK_AddEdge`, `GSSK_DeactivateEdge`, `GSSK_DeactivateNode`, `GSSK_ReclassifyNetwork`
- Replay: `GSSK_Replay`

---

## [3.2.0] — Phase 2: Adaptive Stepping & Serialisation

### Added
- Adaptive RK4 solver: `GSSK_StepAdaptive`, `GSSK_GetLastStepSize`, `GSSK_GetNextStepSize`
- Model serialisation: `GSSK_SerializeModel`, `GSSK_SerializeSnapshot`, `GSSK_FreeString`
- Solver diagnostics: `GSSK_GetSolverConfidence`, `GSSK_GetEdgeErrorEstimate`, `GSSK_GetStepErrorEstimate`
- Diagnostic hooks: `GSSK_SetDiagHooks`

---

## [3.1.0] — Phase 1: Core Kernel

### Added
- C99 kernel: `GSSK_Init`, `GSSK_Step`, `GSSK_Reset`, `GSSK_Free`
- RK4 and Euler integrators
- Edge logic: `linear`, `constant`, `michaelis_menten`, `limit`, `interaction`
- Node types: `storage`, `source`, `sink`
- State/node/edge accessors: `GSSK_GetState`, `GSSK_GetStateSize`, `GSSK_GetNodeID`, `GSSK_GetEdgeID`, etc.
- cJSON bundled parser
- Swift Package: `CGSSK` system module + `GSSK` wrapper
- CLI tool: `bin/gssk <model.json> <output.csv>`
- Regression test suite with CSV comparison

---

[Unreleased]: https://github.com/energese-project/GSSK/compare/v5.3.0...HEAD
[5.3.0]: https://github.com/energese-project/GSSK/compare/v5.2.0...v5.3.0
[5.1.0]: https://github.com/energese-project/GSSK/compare/v5.0.0...v5.1.0
[5.0.0]: https://github.com/energese-project/GSSK/compare/v4.1.0...v5.0.0
[4.1.0]: https://github.com/energese-project/GSSK/compare/v3.6.0...v4.1.0
[3.6.0]: https://github.com/energese-project/GSSK/compare/v3.5.0...v3.6.0
[3.5.0]: https://github.com/energese-project/GSSK/compare/v3.4.0...v3.5.0
[3.4.0]: https://github.com/energese-project/GSSK/compare/v3.3.0...v3.4.0
[3.3.0]: https://github.com/energese-project/GSSK/compare/v3.2.0...v3.3.0
[3.2.0]: https://github.com/energese-project/GSSK/compare/v3.1.0...v3.2.0
[3.1.0]: https://github.com/energese-project/GSSK/releases/tag/v3.1.0
