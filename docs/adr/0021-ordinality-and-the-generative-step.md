# ADR 0021 — Ordinality as Giannantoni defines it, and a generative step under Maximum Em-Power

- **Status**: accepted
- **Date**: 2026-10-09
- **Task**: `adr-ordinality-mop`
- **Supersedes**: the definitions of *ordinality* in [ADR 0014](0014-ordinality-over-quantity-legs.md)
  and [ADR 0015](0015-what-the-emergent-quality-closes.md), and ADR 0015's emergent component `E`
  (decisions 1, 2 and 5). ADR 0014's **leg rules** (which legs carry quantity, how a walk passes
  through a module) are kept unchanged.
- **Depends on**: [ADR 0018](0018-giannantoni-test-protocol.md) — rule 5 (G5): definitions cite sources
- **Records**: PLAN.md §4 R5, R6 and erratum X6; srs.md FR-ORD-001…005
- **Precedes**: `mop-ordinality`, `mop-generative-empower`

## Context

ADRs 0014 and 0015 define ordinality as the fraction of components that lie on a closed pathway, and
"maximum ordinality" as that fraction reaching 1. That is a reasonable graph property, but it is not
Giannantoni's ordinality. ADR 0018 rule 5 forbids defining a Giannantoni term except by citing a
source (PLAN §1 B10). The proxy also gives wrong answers on structures the source decides
differently:

- **Two disjoint 2-cycles.** Every component is on a closed pathway, so the proxy reports maximum.
  But no couple across the two cycles is related at all, so they are not "all the various couples …
  of Ordinality {2/2}" `[22 §12.1]`.
- **Boundary nodes.** The proxy counts sources and sinks as components, so the example
  `closed_loop.json` sits at 0.800 because its heat sink is "open". `[10]` places the surrounding
  habitat outside the system ("including that of the surrounding habitat"), and `[23 Eq 6.3]` treats
  habitat as conditions, not as matrix entries.

ADR 0015's generative step is a graph heuristic: a hub, a new component `E`, and legs to close every
open component. No source contains it (PLAN E7). The sources' generative principle is the Maximum
Em-Power Principle `[02 Eq 5.3]`, `[22 Eq 2]`, and the MOP is "nothing but the re-proposition of" it
`[10 Abstract]`, `[22 §9]`.

## Decision

### 1. Ordinality `{k, (m n)}`, from the exponents the sources give (R5, X6)

`[10]` assigns exponents ½, 2 and 2/2 to co-production, interaction and feedback. `[22 Eq 6–8]` gives
`(d̃/dt)^{1/2}` (binary), `(d̃/dt)^2` (duet) and `(d̃/dt)^{2/2}` (duet-binary). So in `q = m/n` the
numerator counts interactions (the power), and the denominator counts co-productions (the number of
root branches). `[22]`'s wording agrees with this; `[23 Eq 4.1.1]`'s text, "m Co-productions and n
Interactions", reverses it, and is erratum **X6**.

**Components.** Non-module, non-boundary nodes. `source`, `sink` and `constant` are habitat, and
modules hold nothing (ADR 0014).

**Each unordered couple of components is classified by the first rule that applies:**

| Ordinality | Condition |
|---|---|
| `2/2` (duet-binary, feedback) | each reaches the other along quantity-carrying legs. The walk is ADR 0014's, unchanged: a read control carries nothing, and a module passes energy to its products and a drawn control to its used leg |
| `2` (duet, interaction) | both feed one interaction module along quantity-carrying legs |
| `½` (binary, co-production) | both are products of one replicating process |
| unrelated | none of the above |

**The record** is `{k, n₂₂, n₂, n½, n_unrelated}`, with k the number of components `[23 Eq 3.2]`.

**Maximum Ordinality** is `{2 2} ↑ {N N}` `[22 Eq 11.1]`: every couple is 2/2 `[22 §12.1]`. Since 2/2
is mutual reachability, this holds **exactly when the component graph is strongly connected**, and
the engine decides it that way. A model with fewer than two components has no couple and is not at
maximum.

**Closure.** The old measure, the fraction of components on a closed pathway, is renamed `closure`.
It is reported as `closure (proxy)` and decides nothing (FR-ORD-004). `gia_ordinality` keeps its name
for one release, returns closure, and is marked deprecated.

### 2. The generative step: one pathway at a time, by maximum empower (R6)

Below Maximum Ordinality, the step repeats the following until the component graph is strongly
connected:

1. Form the condensation of the component graph (its strongly connected components, as a DAG).
2. **Candidates** are pathways from a component in a *sink* SCC (no out-edge) to a component in a
   *source* SCC (no in-edge).
3. Each candidate is a `linear` pathway with the seed's mean edge weight. For each, compute the total
   empower at `t_end` with the emergy pass. Choose the candidate that **maximises total empower**.
   This is the discrete form of `∫Γφ* dV = d/dt ∫ em* dV → Max` `[02 Eq 5.3]`.
4. Ties are broken by lexicographic `(from, to)` id, never by array order.

**Termination.** Each addition gives a sink SCC an out-edge, or a source SCC an in-edge, or merges
them, so the number of sources plus sinks of the condensation strictly falls. The step therefore ends
in at most `#sources + #sinks` additions. At Maximum Ordinality it adds nothing.

**No new component.** ADR 0015's `E` is retired. `[22 §12.1]` has couples *become* 2/2, which is a
change of relation between existing components. It does not add a component.

### 3. Every example's verdict, before and after

Hand-derived from the walk above. `mop-ordinality` and `mop-generative-empower` assert these rows as
tests, so this table cannot silently go stale.

| Seed | Before (ADRs 0014/0015) | After (this ADR) |
|---|---|---|
| `examples/giannantoni/input.json` | ordinality 0.000, below maximum. Generative: adds component `E` and its legs | components `store_1`, `consumer_1`; record `{2, 0, 0, 0, 1}` (`store_1 → consumer_1` is one-way; `store_1`'s control is read, so the two do not both feed `interaction_1`); **below maximum**; closure 0.000. Generative: adds `consumer_1 → store_1` (the only candidate), which makes the record `{2, 1, 0, 0, 0}`, at maximum |
| `examples/giannantoni/closed_loop.json` | ordinality 0.800 (heat sink counted and open), below maximum. Functional: "fixed point below maximum" | components `store_1`, `consumer_1`, `recycler_1`; all three mutually reachable (`… → recycler_1 → source_1 → interaction_1 → store_1 → …`); record `{3, 3, 0, 0, 0}`; **at maximum**; closure 1.000. Functional: nothing to add |

## Consequences

- Two kinds of test change: those that assert a cycle-coverage *number* as ordinality, and those that
  assert an emergent `E`. vv-plan.md §5 lists them. They are retired or revised in the PR that
  implements this ADR, not deleted silently.
- `examples/giannantoni/closed_loop.json` changes verdict from "below maximum" to "at maximum". The
  heat sink was never part of the system; it was counted only because the proxy counted every
  non-module node.
- The generative step now needs the emergy pass for each candidate, so it costs one emergy evaluation
  per candidate per addition. For N ≤ 64 (NFR-LIM-001) this is small.

## Alternatives considered

- **Keep cycle coverage and call it ordinality.** Rejected by ADR 0018 rule 5, and wrong on two
  disjoint 2-cycles.
- **Choose the generative pathway by a graph rule (shortest closing pathway, highest-degree hub).**
  Rejected. The sources give the selection principle, Maximum Em-Power, and the engine already
  computes empower.
- **Close all open components in one step through a new component.** Rejected (ADR 0015's design).
  The sources describe couples becoming 2/2, not a new component.
