# Stakeholder and business requirements

What the Giannantoni kernel is for, for whom, and how each need is accepted. Every `BR` is verified by
validation tests (`VAL-…`), inspections or demonstrations in [vv-plan.md](vv-plan.md); the functional
requirements that realise it are in [srs.md](srs.md).

## 1. Purpose

A C99 kernel that computes what Giannantoni's Incipient Differential Calculus (IDC) and Maximum
Ordinality Principle (MOP) define — from the equations in the 2002–2023 sources, corrected where a
probe shows the printed form wrong — and says plainly what it computed, with which calculus, and where
the sources stop. It is a research instrument: its value is that a result can be traced to an equation
and checked, not that it produces numbers.

## 2. Stakeholders

| ID | Stakeholder | Needs |
|---|---|---|
| S1 | **Researcher** in emergy, systems ecology or Giannantoni's programme | Solve IDC and MOP problems; compare incipient and traditional results; run emergy accounting on Odum networks |
| S2 | **Sceptical reviewer** checking a claim made with the kernel | Trace any output to a source equation; reproduce published results; see where the sources are wrong or silent |
| S3 | **Maintainer** | A spec to accept work against; a build that cannot be green without the work being done |
| S4 | **Integrator** embedding the library | Stable C API contracts, error codes, reentrancy, no surprise process exits |
| S5 | **Implementing agent** (human or LLM) | Unambiguous requirements, each with an oracle, so "done" is decidable |

## 3. Constraints

- **C1.** C99 (`-std=c99 -Wall -Wextra -Werror`), built by GCC and clang on Linux and Apple clang on
  macOS (`AGENTS.md`).
- **C2.** The Giannantoni units do not include `gssk.h` (ADR 0011).
- **C3.** Native library and CLI only. The Giannantoni engine is not part of the WASM artefact in this
  baseline; WASM remains the GSSK kernel's surface.
- **C4.** Sources are those listed in [docs/sources/README.md](../sources/README.md). A behaviour no
  source defines is not invented; it is refused by name (BR-009).

## 4. Business requirements

### BR-001 — Solve the incipient calculus the sources pose
The kernel shall compute the incipient derivative of integer and fractional order, and solve the
incipient equations the sources work through: the second-order linear equation with variable
coefficients, the binary (half-order) equation, the Riccati equation by its linearisation, the
nonlinear equation of [02 Eq 14.10.1], and the incipient Taylor projection.
- **Source:** [02 App. 7–9], [06 §3], [09], [10]
- **Verification:** VAL-01, VAL-05, VAL-06 · **Priority:** Must · **Status:** planned · **Task:** idc-lde2

Acceptance: every FR-IDC requirement is `implemented`, and VAL-01, VAL-05 and VAL-06 pass.

### BR-002 — Quantify the drift between the two calculi
The kernel shall report the difference between incipient and traditional results — both the solution
drift where the model is known [06 §4] and the output-projection drift where only a trajectory is
known [09 Eq 13] — and shall report zero where the sources say the calculi coincide.
- **Source:** [06 §4], [09 Eq 11–13], [10 Eq 14–16]
- **Verification:** VAL-01, T-IDC-11 · **Priority:** Must · **Status:** planned · **Task:** idc-drift-coupled

Acceptance: FR-IDC-004, FR-IDC-011 and FR-IDC-014 `implemented`.

### BR-003 — Emergy algebra and its ordinal forms
The kernel shall compute emergy and transformity on an Odum network under Odum's rules, and express
co-production, interaction and feedback in Giannantoni's ordinal forms.
- **Source:** [02 Ch. 3], [06b], [22 Eq 6–8]
- **Verification:** VAL-03 · **Priority:** Must · **Status:** planned · **Task:** emergy-source-terms

Acceptance: every FR-EM requirement `implemented`; VAL-03 passes.

### BR-004 — Solve the MOP's First Fundamental Equation
The kernel shall solve the First Fundamental Equation for every couple of a system, from boundary
conditions given in the seed or derived from the network.
- **Source:** [23 Eq 4.1, 5.4.2–5.5.8], PLAN R1, R8
- **Verification:** T-MOP-01, T-MOP-10 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

Acceptance: FR-MOP-001 to FR-MOP-004 and FR-MOP-008 `implemented`.

### BR-005 — Evaluate the Second Fundamental Equation and the EQS
The kernel shall evaluate the printed solution of the Second Fundamental Equation and the operative
EQS form, in the Relational Space algebra as printed.
- **Source:** [23 Eq 6.1–6.3, 7.1–7.5, 5.1.3–5.1.5]
- **Verification:** VAL-02, T-MOP-07 · **Priority:** Must · **Status:** planned · **Task:** mop-eqs

Acceptance: FR-MOP-005, FR-MOP-006 and every FR-REL requirement `implemented`; VAL-02 passes.

### BR-006 — Ordinality, Maximum Ordinality and the generative step
The kernel shall compute a network's Ordinality as the sources define it, decide whether it is at
Maximum Ordinality, and below it perform the generative step under the Maximum Em-Power Principle.
- **Source:** [22 Eq 11.1, §12.1], [02 Eq 5.3], PLAN R5, R6
- **Verification:** T-ORD-02, T-ORD-04 · **Priority:** Must · **Status:** planned · **Task:** mop-ordinality

Acceptance: every FR-ORD requirement `implemented`.

### BR-007 — Report harmony honestly
The kernel shall measure whether a Relational Space satisfies the Harmony Relationships, and classify
each construction as imposing, transporting or not producing harmony — by measurement, never by
assertion.
- **Source:** [23 Eq 5.6.5, App. A1], PLAN R7
- **Verification:** VAL-07 · **Priority:** Must · **Status:** planned · **Task:** mop-harmony-detector

Acceptance: every FR-HAR requirement `implemented`; VAL-07 passes.

### BR-008 — Reproduce published results where the sources allow it
The kernel shall reproduce every published numerical result whose inputs the sources give.
- **Source:** [09 Eq 16–22], [23 Eq 7.1–7.4], [02 Eq 3.23–3.26], [10 App. Eq 28–34], [06 Eq 3.9]
- **Verification:** VAL-01, VAL-02, VAL-04, VAL-05, VAL-06 · **Priority:** Must · **Status:** planned · **Task:** idc-taylor

Acceptance: the five VAL tests pass on every CI toolchain at the printed precision; the two `[09]`
values that do not reproduce are asserted as errata (X4, X5), not skipped.

### BR-009 — Never claim more than the sources support
Every output shall say which calculus or construction produced it; every behaviour the sources do not
define shall be refused by name with the reason; the documentation shall not describe the kernel as
implementing more than its `implemented` requirements.
- **Source:** PLAN §1 (B5), §6
- **Verification:** VAL-09, T-OUT-02, T-OUT-03 · **Priority:** Must · **Status:** planned · **Task:** mop-claims-remediation

Acceptance: FR-OUT requirements `implemented`; VAL-09 inspection passes.

### BR-010 — Every claim checkable from the repository alone
Each requirement shall trace to a source equation and to the tests that verify it; each departure from
a printed equation shall carry a probe and a test that the printed form fails.
- **Source:** PLAN §1 (G4, G5, G7)
- **Verification:** T-TRC-01, VAL-08 · **Priority:** Must · **Status:** planned · **Task:** guard-trace

Acceptance: `make check-trace` passes in CI; VAL-08 passes.

### BR-011 — An embeddable, dependable library
The library shall be reentrant, deterministic, free of process exits, and portable across the CI
toolchains, with contracts an integrator can program against.
- **Source:** `AGENTS.md` §Architecture, §Coding Standards; C1
- **Verification:** T-REE-01, T-REE-02, T-ERR-01, DEM-POR-01 · **Priority:** Must · **Status:** planned · **Task:** guard-reentrancy

Acceptance: every NFR requirement `implemented`.
