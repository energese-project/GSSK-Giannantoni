# ADR 0018 — The Giannantoni test protocol: no golden files, sourced oracles, enforced gates

- **Status**: accepted
- **Date**: 2026-10-09
- **Task**: `adr-giannantoni-test-protocol`
- **Depends on**: [ADR 0011](0011-two-engines-declared-lossy-projection.md) — two engines;
  [PLAN.md](https://github.com/energese-project/GSSK-Giannantoni/blob/main/PLAN.md) §1 (loopholes B1–B12, guardrails G1–G7);
  [docs/requirements/](../requirements/README.md) — the baseline and its V&V plan
- **Precedes**: every W2–W8 code PR in PLAN.md §7

## Context

PLAN.md §1 audited the repository, at `4e96569`, for ways a build can go green without the
Giannantoni mathematics being implemented. It found twelve (B1–B12). They are not hypothetical. B5
(documents that equate IDC with the matrix exponential) is how a classical solver came to be labelled
"incipient" (PLAN E2), and B11 (a test that checks a construction against itself) is how the harmony
matrix came to be reported as verified (E5).

The kernel's testing practice is right for the kernel and wrong for this work. For the kernel, a
golden file pins numerical stability: the behaviour is specified elsewhere (Odum), and the regression
suite catches unintended change. For the Giannantoni units the behaviour is what is under
construction. A golden file made from the first implementation turns that implementation into its
own oracle (B2). Whatever it printed becomes "expected", including its mistakes.

The sources add a second difficulty: they contain errors. PLAN §5 registers twelve (X1–X12). Some
have probes showing that the printed equation does not solve the equation it claims to solve. A
protocol that let an implementation silently depart from a printed equation would invite the
converse failure: "fixing" the source until a test passes.

## Decision

Every Giannantoni task (PLAN W2–W8, and any later work on `src/engine.c`, `src/idc.c`, `src/mop.c`,
`src/relational.c`, `src/harmony.c`) follows these rules. Each states which loophole it closes and the
check that enforces it.

1. **G1 — No golden files.** A test's oracle is a source equation's residual, a hand-derived closed
   form, a number printed in a source, or an independent computation in a probe script
   (`docs/sources/probes/`). `make test-update` never writes a Giannantoni oracle. In `make test`, a
   model with no expected file fails unless `tests/skip_allowlist.txt` names it, with a reason.
   *Closes B1, B2. Enforced by `guard-no-skip`.*
2. **G2 — Coverage.** The Giannantoni units are in the coverage build, with a gate of ≥ 90% line
   coverage. A coverage figure that cannot be parsed fails the gate. The parse uses portable `awk`,
   not `grep -P`. *Closes B3. Enforced by `guard-coverage-giannantoni` (NFR-COV-001).*
3. **G3 — No stubs in planned API.** Every function declared in `idc.h`, `mop.h` and `relational.h`
   is called from a test, and no body returns a not-implemented code. `AGENTS.md`'s rule permitting
   `(void)` stubs does not apply to these headers. *Closes B4. Enforced by
   `scripts/check_api_called.sh` (`guard-api-called`, NFR-API-001).*
4. **G4 — The errata protocol.** A departure from a printed equation needs all three of:
   - a probe in `docs/sources/probes/`;
   - a row in PLAN §5;
   - a test asserting that the *printed* form fails the oracle (VAL-08).

   No tolerance in [vv-plan.md](../requirements/vv-plan.md) §3 is widened without such a row and a
   PLANLOG entry. *Closes the converse failure described under Context.*
5. **G5 — Definitions cite sources.** An ADR may define a Giannantoni term (ordinality, harmony,
   incipient, drift) only by citing a source equation or a PLAN §5 erratum. It may not redefine a term
   so that tests pass. *Closes B10.* ADR 0014's and ADR 0015's definitions of ordinality are
   superseded on these grounds by ADR 0021.
6. **G6 — Wording.** No document equates IDC with the matrix exponential. `AGENTS.md`'s clamp
   (`Q < 0 → 0`) applies to the GSSK kernel only. The Giannantoni units never clamp, floor or take
   absolute values of signed or complex coordinates, and a test checks that they survive
   (NFR-NUM-006, T-NUM-04). *Closes B5, B6. Enforced by `guard-agents-wording`.*
7. **G7 — TDD evidence.** Each code PR body:
   - maps each test to the equation and requirement it verifies;
   - quotes the red run from before the implementation;
   - quotes the failure of each catalogue mutation, run by hand.

   Each test carries a `Verifies:` tag. A requirement becomes `implemented` only when a tagged test
   passes. *Closes B8, B11. Enforced by `make check-trace` (NFR-TRC-001) and by review.*

Three further rules close the remaining loopholes:

8. **Refusal is bounded.** A feature may be refused (`GIA_E_UNSUPPORTED`) only if it is listed in
   PLAN §6, with the reason the sources do not define it. A refusal not in §6 is an unimplemented
   requirement, and it stays `planned`. *Closes B7, B9.*
9. **A suite that is not in CI does not exist.** The PR that creates a Giannantoni test binary adds its
   target to `CI_TESTS` and to `deploy.yml`, with `.PHONY` declared beside the rule. *Closes B12.*
10. **A tautology is not a test.** A test whose oracle is computed by the code under test is a test of
    that code's consistency, not of the theory. It is named as such: for example, `test_harmony`
    becomes a test of the constructor (FR-HAR-004). *Closes B11.*

## Consequences

- Tests for the Giannantoni units are slower to write: each needs a derivation or a citation. That
  cost is the point. A test that is cheap because it compares against the previous output verifies
  nothing about the theory.
- Some tests assert failure: the printed form of an erratum must *not* pass. A reader who sees a
  test expecting a large residual should look for the §5 row it cites.
- The kernel's golden-file workflow is unchanged for the kernel. `make test-update` remains the way
  to pin a kernel model, and the allowlist in G1 is the only change to it.
- These rules constrain this repository's own agents as much as anyone. `AGENTS.md` points here
  (`guard-agents-wording`).

## Alternatives considered

- **Golden files plus review.** Rejected. Review of a 10-digit CSV cannot tell a correct number from
  a plausible one, and B2 is exactly that workflow.
- **Leave the protocol in PLAN.md only.** Rejected. A plan is revised and then archived; an ADR is
  the repository's record of a decision that binds later work, and G5 itself constrains ADRs.
