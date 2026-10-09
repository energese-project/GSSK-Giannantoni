# ADR 0022 — The kernel's matrix-exponential method is named `expm`; `incipient` is a deprecated alias

- **Status**: accepted
- **Date**: 2026-10-09
- **Task**: `kernel-method-label`
- **Depends on**: [ADR 0011](0011-two-engines-declared-lossy-projection.md) — the kernel is the classical
  (TDC) engine; [ADR 0018](0018-giannantoni-test-protocol.md) rule 6 (G6) — wording
- **Records**: PLAN.md §4 R15; srs.md FR-KER-001

## Context

The GSSK kernel accepts `"method": "incipient"` (`src/gssk.c`, `GSSK_METHOD_INCIPIENT`). What it runs
is a Padé (3,3) matrix exponential, `exp(A·dt)`, with `limit` edges linearised through an effective
conductance (PLAN §2 E2). That is a sound classical method, but it is not Giannantoni's incipient
calculus. The name is how the repository came to believe it implemented IDC (PLAN §1 B5): a user
choosing `"incipient"` reasonably expects the incipient calculus, and the schema's description says
"force IDC".

ADR 0011 already assigns the kernel the classical role and the Giannantoni engine the incipient one.
The method name contradicts that decision.

## Decision

1. **`"method": "expm"` is the documented name** of the kernel's matrix-exponential method. The
   schema enum gains `"expm"`, and its description says what the method computes.
2. **`"incipient"` stays accepted as a deprecated alias.** It selects exactly the same code path,
   with byte-identical output, and the kernel prints one notice per load on stderr, saying that the
   method is the matrix exponential and naming `expm`. The notice goes to stderr, not to the CSV, so
   no output a consumer parses changes.
3. **The kernel offers no incipient method.** Per ADR 0011, it is the classical engine. The incipient
   calculus is built in the Giannantoni engine (PLAN W2) and nowhere else.
4. **This is not a schema break.** Every model that loaded before loads after, and runs to the same
   numbers. `GSSK_METHOD_INCIPIENT` keeps its enum value. `GSSK_METHOD_EXPM` is added as an alias with
   the same value. Serialisation writes back the spelling the model was loaded with, so a model's
   round trip is unchanged.
5. **Verification** (T-KER-01): the same model run with `"incipient"` and with `"expm"` produces a
   byte-identical CSV, and the `"incipient"` run prints the notice exactly once.

## Consequences

- The one observable change for an existing model is the notice on stderr, and the CHANGELOG lists
  it under "Changed".
- A later release may remove the alias. This ADR does not schedule that; removal would be a schema
  break and needs its own decision.

## Alternatives considered

- **Leave the name and fix the documentation.** Rejected. Every model file that says `"incipient"`
  repeats the claim to whoever reads it.
- **Remove `"incipient"` outright.** Rejected. It breaks every model that uses it, for a renaming.
