# ADR 0012 — Where a law lives: pathway or module

- **Status**: accepted — amended 2026-09-11 to decide that module laws leave pathways (decision 5)
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

**5. Module laws leave pathways.** *(Amendment, 2026-09-11.)*

As first written, this ADR left open whether the pathway-level module laws
survive once module-hosted laws exist, and said the question should be settled
before any module code rather than after — because building the module forces
an answer whether or not anyone gives one, and whatever the implementation
happened to do would become the decision by default.

Settled: **they are removed.** Pathways carry Odum §III and nothing else —
`linear`, `reversible` and `constant`. Seven laws become module-hosted:
`interaction`, `limit`, `ratio`, `subtract`, `gain`, `threshold` and
`exchange`. There is one way to write a work gate, and it is the way Odum draws
one.

The alternative weighed was to keep both spellings behind one shared evaluator,
the edge form being the module with one energy input and one output. That is
what ADR 0001 decided for the transaction diamond — *"Keep both authoring
forms. Deprecate neither. Back them with one shared primitive"* — and it would
have broken nothing. It is not followed here, for a reason specific to where
each decision was made. ADR 0001 was taken in GSSK, against a published schema
with existing models and users, where a breaking change has a real cost paid by
other people. This repository is young, its two example models are its own, and
it has no external users; the cost of a single spelling is a migration done
once, while the benefit — that the two forms cannot drift apart, because there
are not two forms — is permanent. The balance ADR 0001 struck is right for
GSSK and is the wrong one here.

### Emergence is not desugaring

Removing the edge forms invites a confusion worth heading off, because both of
the things involved can be described as "synthesising a node".

**Emergence** is a node produced by the *dynamics*. In this engine, a graph
below maximum ordinality grows a component to close an open pathway, and that
component is Giannantoni's emergent quality. It means something, and it must
appear in the output: `output.json` differing from the seed is precisely how a
generative run is detected.

**Desugaring** is a node invented by the *loader*, as a device for representing
something the author wrote a different way. It means nothing. It is the thing
ADR 0001 rejected — a synthesised node leaks into the CSV as a state column
nobody wrote, appears in serialised output, and gets in the way of runtime
mutation.

They are not in tension, and there is a reason specific to this framework to
keep them firmly apart. Ordinality is computed as the components on a closed
pathway over the total number of components, and ordinality is what decides
**whether emergence happens at all**. A bookkeeping node would count toward that
ratio. Desugaring edge laws into hidden gate nodes would let an implementation
device vote on whether the system generates — the one place in the engine where
bookkeeping must have no say.

So removal is carried out by *migrating models* to write modules explicitly,
never by quietly converting edge laws into modules at load time.

## What this does not decide

- **The concrete schema** for declaring a module and its role assignments.
- **How `exchange` finds its counter-flow legs once it is a module.**
  Decision 5 makes it one. As a pathway law it needs its legs named, because
  which component pays and which receives cannot be recovered from carrier
  identity alone — the reason the carrier work declined to infer them. As a
  module it can read its whole neighbourhood, which is exactly the ownership
  the diamond's shape encodes (ADR 0001). So module form is what makes leg
  discovery possible at all. Whether to do it, and by what rule, is not
  settled here.

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

**Migration, from decision 5.** Seven pathway laws become modules. In the
suite, `exchange` accounts for four of the test models that write a module law
on an edge, `threshold` for two, and `interaction`, `ratio` and `subtract` for
one each; both models under `examples/giannantoni/` use `interaction` through
their `generative_production` and `ordinal_feedback` labels. `linear`, by far
the most used pathway law, is untouched. The migration is to write these as
modules explicitly — not to have the loader rewrite them, for the reason given
under *Emergence is not desugaring*.

**The projection has to translate, and must say so.** GSSK writes
`interaction`, `limit`, `threshold`, `ratio` and `subtract` on edges. Reading a
GSSK model therefore means turning each such edge into a gate component the
GSSK author did not write. That is a node added in translation, and it is
legitimate for the same reason desugaring is not: it happens in a new document
in a different vocabulary, and the coverage report records it openly rather than
presenting the result as though the node had always been there. It also moves
node counts, and so ordinality, which is one more reason the report must name
it.

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
