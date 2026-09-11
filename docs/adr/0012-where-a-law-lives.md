# ADR 0012 — Where a law lives: pathway or module

- **Status**: proposed
- **Date**: 2026-09-11
- **Task**: `adr-where-laws-live`
- **Supersedes**: nothing
- **Blocks**: nothing

## Context

The projection built for [ADR 0011](0011-two-engines-declared-lossy-projection.md)
put a number on how much of a GSSK model this engine can carry: **92.3% mean
over the 24 models in `examples/`**. The largest single blocker, 3 of those 24,
it reported as

> a processing node carries its law in node params, and this engine puts laws
> on pathways

That diagnosis is true and too shallow. Reading `compute_interaction_node` in
`src/gssk.c` shows what is actually different:

```c
double F = inst->nodes[ni].node_k;
for (each active incoming edge)   F *= state[e->origin_idx];
if (energy_orig >= 0)             deriv[energy_orig] -= F;   /* only the FIRST */
double share = F / (double)out_count;                        /* over ALL outputs */
```

Three things, none of which is about where parameters are written down:

1. **Arity.** `F` is a product over *every* incoming edge. This engine's
   `interaction` is `F = k·Q_origin·Q_control` — binary, with at most one named
   control.
2. **Asymmetry among inputs.** The first incoming edge is the energy input and
   is consumed; the rest are controls, read but not depleted. A pathway law has
   one origin, so it has nowhere to put "these three are read and that one is
   drained".
3. **Fan-out.** The output is partitioned evenly over *every* outgoing edge. A
   pathway has one target.

So a GSSK processing node is a **hyperedge**: one relation over many
components, whose parameters naturally live on the relation rather than on any
of its endpoints. This engine's pathway is a **binary relation**. That is a
difference in the shape of the graph, not in where a `k` is stored.

It is also the same root as a blocker the projection reported separately —
`three_input_gate_model.json`, blocked on `control_nodes` having two entries
([ADR 0008](0008-nary-interaction-and-subtracting-action.md)). Both are arity.
The projection counted them as two findings because it reads syntax; they are
one problem.

### Odum puts laws in two places, and so should this

The temptation is to pick a side. Odum does not.

- **The pathway carries a law.** Odum 1972 §III distinguishes the barbed
  pathway, where flow depends only on the force behind it, from the barb-less
  one, where it depends on the difference between the forces at each end. That
  distinction lives on the *line*, not in any symbol, which is why it is drawn
  as a property of the line — and it is why `reversible` is an edge logic here
  ([ADR 0007](0007-reversible-pathway.md)).
- **The module carries a law.** The work gate (§X), constant gain amplifier
  (§IX), switch (§XI), cycling receptor (§XIII) and economic transactor (§XV)
  are *modules*. The transformation happens inside the symbol. Odum gives each
  its own mathematical definition, and the volume editor introduces the set as
  "mutually exclusive and exhaustive" — a claim that the typing is load-bearing.

This engine collapsed the module laws onto pathways. That is exact for the
binary two-node case, which is why `interaction`, `limit`, `gain` and the rest
work and are tested against closed forms. It cannot express n-ary inputs,
asymmetric inputs, or fan-out, and no amount of care with a binary relation
will make it.

### The n-ary product is not blocked on new mathematics

Worth settling now, because it removes the reason this looks harder than it is.

The flow matrix already linearises a work gate by folding the control quantity
into the conductance: `g = k·Q_control`, so that `F = g·Q_origin` is linear in
the state being drained. The n-ary case is the same move with more factors:

```
F = k · Q_energy · Π Q_control        ⟹     g = k · Π Q_control
```

All controls but the energy input fold into `g`. Fan-out is a coefficient:
each outgoing pathway takes `g·Q_energy / out_count`. Both fit the existing
augmented matrix without a new solver, and both keep the property that matters
— `psi` is exactly zero when the matrix is constant, non-zero at the
multiplicative junction, which is precisely where a work gate is.

So the decision below is a modelling decision, not a numerical one.

## Decision

**1. A law lives where Odum draws it.**

Pathway laws express Odum §III and nothing else: `linear` (barbed) and
`reversible` (barb-less), plus `constant` for a rate that reads no state.

Module laws are hosted by the component that *is* the module — the work gate,
amplifier, switch, cycling receptor and transactor of §IX–§XV — and are
configured on it, as GSSK already does.

**2. A module is a hyperedge over its neighbourhood, not a node with extra
fields.**

Its law reads every incoming pathway, distinguishes the energy input from the
controls, and partitions its output over every outgoing pathway. That is the
whole of the difference, and writing it down as a hyperedge rather than as
"a node with a `params` block" is what keeps the arity and fan-out visible.

**3. The projection reports one finding, not two.**

`control_nodes` with several entries and a configured processing node are the
same limitation seen from two syntaxes. The coverage report should say so,
because the current output implies two independent problems and overstates how
scattered the gap is.

**4. Nothing about `psi` or the exactness guarantee changes.**

An n-ary gate is a multiplicative junction, so it makes the flow matrix
state-dependent exactly as a binary one does. The composition and its reported
integration error already handle that case.

## What this does not decide

- **Whether the pathway-level `interaction`, `limit`, `ratio` and `subtract`
  logics are kept.** They are the binary special cases of module laws, they are
  shipped, tested against closed forms, and used by every example under
  `examples/giannantoni/`. Removing them is a breaking change that wants its
  own ADR and its own migration; keeping both is the two-spellings drift ADR
  0008 argued against. This ADR deliberately does not resolve that tension, and
  it should be resolved before the module work starts rather than after.
- **The concrete schema** for declaring a module and its role assignments.
- **Whether `exchange` becomes a module here.** It is currently a pathway law
  with named counter-flow legs, which works and is tested (Eq 103 recovered from
  the trajectory), but leg discovery from the diamond's shape is exactly the
  neighbourhood-reading a module does naturally. See ADR 0001.

## A finding for GSSK, not a decision for this repository

GSSK admits an interaction in **both** places: as `GSSK_NODE_INTERACTION`,
configured through the node's `params`, and as `GSSK_LOGIC_INTERACTION` on an
edge. The two are not equivalent — the node form is n-ary with a consumed first
input and partitioned output; the edge form is binary — so they are not two
spellings of one thing. But a reader cannot tell that from the names, and
choosing between them requires knowing which semantics the kernel attaches to
each.

ADR 0008 rejected a `params.op` field on precisely this ground: *"an `op` field
would mean `logic: "ratio"` and `logic: "interaction", op: "div"` are two
spellings of one thing — the drift this repository's standards exist to
prevent, arriving in the published schema where it cannot be corrected
quietly."* The node/edge interaction pair is the same shape of hazard, arrived
at from a different direction. Recorded here because this repository noticed it;
it is GSSK's to act on.

## Consequences

**The coverage number should rise, and by a knowable amount.** Four of the 24
models are blocked by this and by ADR 0008's n-ary case together. Nothing else
in the measured set turns on it.

**`docs/odum_1972_conformance.md` §1 gains a distinction it currently lacks.**
It scores modules as implemented or not without recording *where* the
implementation put the law, which is why this difference survived several
passes over that document unnoticed.

**Risk: a module is a second place for behaviour to live, and two places is how
vocabularies drift.** The mitigation is the one Odum's own set uses — modules
are a closed, enumerated set with individual mathematical definitions, not an
open extension point. If a sixth module is ever wanted, that is an ADR, not a
config field.

**Risk: hosting a law on a component makes its behaviour depend on the
neighbourhood.** A pathway law can be read locally; a module law cannot, because
adding an edge to a gate changes the gate's arity. That is a real cost, and it
is the reason `compute_interaction_node` has to scan the edge list twice. It is
accepted because it is what Odum's diagram means: the arrows into a work gate
are inputs to *it*, not independent facts about the pairs they connect.
