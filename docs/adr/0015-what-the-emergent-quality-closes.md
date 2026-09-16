# ADR 0015 — What the emergent quality closes

- **Status**: accepted
- **Date**: 2026-09-13
- **Task**: `mop-step-closes-nothing`
- **Supersedes**: nothing
- **Depends on**: [ADR 0012](0012-where-a-law-lives.md) — where a law lives;
  [ADR 0013](0013-module-roles-not-position.md) — roles, never position;
  [ADR 0014](0014-ordinality-over-quantity-legs.md) — ordinality over legs that carry quantity

## Context

Mode 2 is the engine's generative mode. When a graph is below maximum
ordinality, `gia_generate` grows one component — the emergent quality, `E` — to
close a relationship that is open, and the evolved graph differing from the seed
is how a generative run is detected. The whole mode rests on one promise, which
the existing seed test states in a comment: the step "must actually raise
ordinality, not merely add something."

It does not keep that promise. Today the step:

1. takes the **first** component, in array order, that is not on a closed
   pathway, and calls it `open`;
2. takes the component with the most arriving flows, the **first** on a tie, and
   calls it the hub;
3. wires `hub → E → open`, and prints *"closed the loop"*.

That closes `open` only if `open` can already reach the hub. The seed passes
because its open component is a source feeding the loop, which does reach the
hub. The seed test is the only test of the step.

### Reproduction 1: the step moves away from maximum ordinality, and never stops

Each model below was run through the step repeatedly, and ordinality measured
before each step. The same five models were run under the rule this ADR decides.

| model | today | decided rule |
|---|---|---|
| seed (`input.json`) | 0.750 → 1.000, stops | 0.750 → 1.000, stops |
| dead end: `a ⇄ b`, `a → c` | 0.667 → 0.500 → 0.400 → 0.333 → 0.286 → … | 0.667 → 1.000, stops |
| isolated: `a ⇄ b`, `c` | 0.667 → 0.500 → 0.400 → 0.333 → 0.286 → … | 0.667 → 1.000, stops |
| chain: `s → x → y → z` | 0.000 → 0.600 → 0.500 → 0.429 → 0.375 → … | 0.000 → 1.000, stops |
| heat sink: `a ⇄ b`, `a → heat` | 0.667 → 0.500 → 0.400 → 0.333 → 0.286 → … | 0.667, stops — see decision 3 |

In four of five models the step lowers ordinality every time and never
terminates. For a dead end `c`, wiring `a → E → c` leaves `c` with nothing
leaving it, and adds `E` as a second component with nothing returning to it. The
step reports *"closed the loop"* each time.

ADR 0014 surfaced this. Before it, the cycle scan walked control legs, which
marked an accumulator in a module model as already closed, so the step never
reached it.

### Reproduction 2: which component is closed depends on array order

One model — `a ⇄ b`, `a → c`, and an unconnected `d` — with its nodes listed in
two orders:

| | order `c, d, a, b` | order `b, a, d, c` |
|---|---|---|
| today | `a → E → c`, ordinality 0.400 | `b → E → d`, ordinality 0.400 |

A different component is closed, from a different hub, because of the order of
lines in a file. ADR 0013 rejected exactly this for module inputs.

### Reproduction 3: fixing the direction is not enough

The obvious repair adds whichever direction is missing — `open → E → hub` when
`open` cannot reach the hub, `hub → E → open` when the hub cannot reach `open`.
A prototype of that terminates on every model above, and gets two things wrong:

- **It closes the heat sink** by wiring `heat → E → a`, returning degraded
  energy to do work. The engine's own description of a sink is Odum 1972 §V: it
  "absorbs used energy and is never drained".
- **It is still decided by array order.** With the two node orders above it
  wired `c → E → a` in one and `b → E → d`, `d → E → b` in the other.

## Decision

1. **The emergent quality closes every open component, in one step, by adding
   only what is missing.** With hub `h`, for each open component `o`:

   - if `h` cannot reach `o`, add `E → o`;
   - if `o` cannot reach `h`, add `o → E`;
   - add `h → E` if any `E → o` was added, and `E → h` if any `o → E` was added.

   Every such `o` is then on a closed pathway: `h → E → o ⇝ h`, or
   `o → E → h ⇝ o`, or both. `E` is on one too. Adding legs only adds
   reachability, so what the scan found before the step still holds after it. If
   the hub is itself open and nothing else is, the step adds `h → E` and
   `E → h`.

2. **One step never lowers ordinality, and a second step changes nothing.** After
   one step every closable component is on a closed pathway, so the step reaches
   its fixed point at once. This is asserted, not inferred: the evolved graph is
   reloaded and rescanned, for every model in Reproduction 1.

3. **A sink is never closed, and is never the hub.** Closing a sink means drawing
   from it, and drawing from it recycles degraded energy. When the only open
   components are sinks, the step changes nothing and says why. The heat-sink
   model therefore stops at 0.667, below maximum. That is a stated fixed point,
   not a failure to converge, and the report says which case it is.

4. **The step is decided by the model, not its serialisation.** The open
   components are a set, not the first one found. Hub ties are broken by id. The
   legs the step adds are written out in id order. The same model with its nodes
   or edges reordered **appends a byte-identical component and legs**. The rest
   of the output is the seed, copied in the order it was given, so two orderings
   of a seed still serialise differently there — the claim is about what the
   step adds, which is the only part it decides. The prototype already added the
   same set of legs in both orders of Reproduction 2 and differed only in the
   order it wrote them, which is what the last requirement removes.

   *Corrected on implementation.* As first merged this decision said the whole
   output was byte-identical, which cannot hold for a copy of a reordered seed.

5. **The emergent quality is a component, not a module.** Under ADR 0014 a
   module is not counted toward ordinality and its control closes nothing. An `E`
   written as a module could close nothing, and would not count toward the
   ordinality it exists to raise. Its legs are Odum §III linear pathways — the
   `ordinal_ascent` and `emergent_feedback_loop` labels already load as `linear`,
   which ADR 0012 keeps. Its type is `storage`, not `gain`: since ADRs 0012 and
   0013, a gain is a module whose second input is a control. `emerged_from`
   becomes an array of the ids it closes, in id order, because one step may now
   close several.

6. **The report claims closure only after checking it.** The step rescans the
   graph it produced and prints *"closed"* only if the rescan confirms it. If the
   rescan disagrees, the step returns the seed unchanged and reports that
   instead. Decision 1 means this should never fire. The check stays anyway,
   because a claim the engine never verified is the defect this ADR records.

## What this does not decide

- **Whether sources and sinks count toward ordinality at all.** Both are Odum's
  boundary: §II holds a source at its value, §V never drains a sink. Today both
  count, so a source is closed by feedback into a held value — the seed does
  this — and, after decision 3, a model with a heat sink can never reach 1.0.
  Leaving the boundary out would let such a model reach maximum and would stop
  the step feeding back into sources. It would also change every seed's
  ordinality, and it is its own decision.
- **How the hub is chosen,** beyond the tie-break. "Most arriving flows" is kept.
  Preferring a hub already on a cycle, or one emergent quality per disconnected
  cycle, is not decided.
- **The weights of the emergent legs.** They stay the mean pathway weight and
  `1 − ordinality`, which are templated rather than derived. Deriving the
  emergent quality is already listed as not implemented.
- **Whether closing every open relationship at once is Giannantoni's reading.**
  This is the engine's rule, not a claim about his algorithm. Closing one
  component per step would need a rule for which goes first. Every such rule
  found is arbitrary — array order, id order — and array order is what
  Reproduction 2 rejects. Closing them all at once is the only choice found that
  needs no such rule.

## Consequences

- The seed still has one open component and gains the same two legs, with the
  same weights. Its `output.json` still differs from before: `E`'s type becomes
  `storage`, `emerged_from` becomes an array, and — so the output does not call a
  storage a gain — its id becomes `emergent_quality_N` and its label
  "Emergent Quality".
- `examples/giannantoni/closed_loop.json` is a hand-written evolved graph whose
  emergent component is a `gain`. It still loads, and its ordinality is
  unchanged. Rewriting it belongs to the model migration, which now depends on
  this ADR.
- Tests: each model in Reproduction 1 reaches its fixed point in one step, and a
  second step is functional. Ordinality never decreases across a step. Both
  node orders append byte-identical components and legs. No sink appears on any
  emergent leg.
  Test [45], which reads `emerged_from` as a string, changes with decision 5.
