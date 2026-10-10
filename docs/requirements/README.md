# Requirements baseline — Giannantoni kernel

This directory specifies what the Giannantoni kernel (`bin/giannantoni_sim` and the library units
behind it) must do, and how each requirement is proven. [PLAN.md](https://github.com/energese-project/GSSK-Giannantoni/blob/main/PLAN.md)
decides *how and in what order* to build it; this baseline decides *what counts as built*. Where the
two disagree, this baseline wins and the plan is revised.

| Document | Contents |
|---|---|
| [stakeholder.md](stakeholder.md) | Stakeholders, business requirements (`BR-nnn`) and their acceptance criteria |
| [srs.md](srs.md) | Functional (`FR-…`) and non-functional (`NFR-…`) requirements |
| [icd.md](icd.md) | Interface control: C API contracts, the seed's `mop` block, outputs, CLI (`IF-…`) |
| [numerics.md](numerics.md) | Numerical design for each algorithm: formulation, accuracy, failure modes (`N1`…) |
| [vv-plan.md](vv-plan.md) | Verification and validation plan, tolerance and determinism policy, and the test catalogue |

## Verification and validation, distinguished

The specification here is a theory with known errata, so the two are kept apart:

- **Verification** — the code implements this baseline. The baseline states the *corrected*
  equations (PLAN §4–5 and the errata register), so a verification test's oracle is a source
  equation's residual, a hand-derived closed form, or an analysis.
- **Validation** — this baseline is what Giannantoni meant. It is shown by reproducing results the
  sources print (`VAL-…` in the catalogue), and by showing each erratum's printed form fails the
  defining equation it claims to solve.

## Conventions

**Requirement IDs** are `BR-nnn`, `FR-<AREA>-nnn`, `NFR-<AREA>-nnn` and `IF-<AREA>-nnn`. IDs are never
reused or renumbered; a withdrawn requirement keeps its heading with `Status: withdrawn`.

| Area | Meaning |
|---|---|
| `IDC` | Single-variable incipient calculus |
| `EM` | Emergy algebra and its ordinal forms |
| `MOP` | The two Fundamental Equations, EQS, network mapping |
| `REL` | Relational Space algebra |
| `ORD` | Ordinality and the generative step |
| `HAR` | Harmony detection and verdicts |
| `OUT` | Labelling and reporting |
| `KER` | The GSSK kernel's method label |
| `NUM`, `DET`, `REE`, `MEM`, `ERR`, `LIM`, `ROB`, `POR`, `PERF`, `SEP`, `TRC`, `COV`, `API` | Non-functional areas, listed in srs.md §3 |
| `API`, `JSON`, `OUT`, `CLI` (in `IF-`) | Interface areas |

**Every requirement** is a `###` heading `ID — title`, a statement using **shall**, then the fields

```
- **Source:** …                        equation, erratum, or PLAN decision it rests on
- **Verification:** <IDs> · **Priority:** Must|Should|Could · **Status:** planned|implemented|withdrawn · **Task:** <crux slug>
```

**Verification IDs** name entries in the catalogue in [vv-plan.md](vv-plan.md): `T-…` (test),
`VAL-…` (validation test), `INS-…` (inspection), `ANA-…` (analysis), `DEM-…` (demonstration).

**Status `implemented`** means a test in the repository carries a `Verifies:` tag naming the
requirement and passes on `main`. It is a claim `make check-trace` checks, not a judgement.

**Tagging tests.** A test that verifies requirements says so in a comment immediately above it:

```c
/* Verifies: FR-IDC-001, FR-IDC-004 (T-IDC-01) */
static void test_drift(void) {
```

## The trace check

`make check-trace` (`scripts/check_trace.sh`) fails when:

1. a requirement ID is duplicated, or malformed;
2. a requirement has no `Verification:` IDs, or names one that is not in the catalogue;
3. a catalogue entry verifies no requirement, or names one that does not exist;
4. a requirement is `implemented` but no `Verifies:` tag under `tests/` or `scripts/` names it;
5. a `Verifies:` tag names a requirement that does not exist.

It runs in CI beside `make test-giannantoni`. It cannot tell whether a test is *good*; that is what the
oracle policy in vv-plan.md §2 and the mutations in the catalogue are for. It makes an untested
requirement, an orphan test, or an unbacked `implemented` claim impossible to merge silently.

## Changing the baseline

A change to a requirement is a PR that edits it here, explains why in
[PLANLOG.md](https://github.com/energese-project/GSSK-Giannantoni/blob/main/PLANLOG.md), and — when it
changes engine semantics — lands its ADR first. A departure from a printed source equation needs a
probe in `docs/sources/probes/` and an entry in PLAN §5 before any requirement may rely on it.
