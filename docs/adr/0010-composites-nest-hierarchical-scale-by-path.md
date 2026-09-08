# ADR 0010 — Composites nest: hierarchical scale by path, not by flattening

- **Status**: proposed
- **Date**: 2026-09-07
- **Task**: unassigned — no crux task yet
- **Supersedes**: nothing
- **Blocks**: nothing

## Context

Odum's diagram vocabulary has two tiers. There are **fundamentals** — source,
storage, sink, interaction, gain, loop_limited, exchange, switch — and there are
**composites** that are nothing but recurring arrangements of fundamentals:
`producer`, `consumer`, `misc_box`.

The tiers are not fixed by nature. `loop_limited` is the proof: Odum's D-shape
is conceptually a store plus an interaction plus an internal recycle loop, and
he promoted it to a symbol of its own because he drew it often enough that
keeping it compositional stopped paying. Naming is lexicalisation — a
combination earns a symbol by recurring, exactly as a frequent collocation
earns a dictionary entry. Phase 9 already automates this at runtime: motifs
that recur and persist become archetype candidates, and `GSSK_ProposeArchetype`
promotes one to a named type.

So the vocabulary is open-ended by design. The representation is not.

### What Odum 1972 actually defines

The definitive early statement is Odum, H.T. (1972), *An Energy Circuit Language
for Ecological and Social Systems: Its Physical Basis*, in Patten (ed.), *Systems
Analysis and Simulation in Ecology* Vol. II, ch. 4 — held in `docs/`. Odum states
the language as *"a dozen basic modules, each having a mathematical definition"*
(§I, p.141), and the volume editor introduces them as *"a mutually exclusive and
exhaustive set sufficient for the functional representation of macroscale
systems"*.

Three of those modules bear directly on this ADR:

- **§XII Self-Maintaining Module** (Fig. 1g) — *"potential energy input being in
  part restored against potential generating forces, and then routed to increase
  flow from the upstream source by at least one kind of work on the upstream
  flow."* Growth is logistic; the multiplicative feedback loop is *"autocatalytic
  by chemists and logistic by population biologists."* Specifying it needs a
  storage (`C`), a work gate (`k`), and the ratio of potential-generating work to
  dissipational flows.
- **§XIII Cycling Receptor Module** (Fig. 1e) — multiplicative positive feedback
  *"with the additional constraint of being dependant on the recycling of a
  material which is constant over short periods"*, giving *"output as a hyperbolic
  function of the input, the output leveling off as the recycling process becomes
  limiting."* Odum's examples are the photocell and **photosynthesis**.
- **§XIV Production and Regeneration Module (P-R)** (Fig. 1h) — *"the module
  formed by combining a cycling receptor module, a self-maintaining module which
  it feeds, and a work feedback loop which controls the inflow process by
  multiplicative and limiting actions. **An example is the green plant**, which
  has a respiratory system."*

So Odum's own green plant is **P-R**, and he defines it explicitly as a
combination of two other named modules. That is a composite of composites, in
the source text, in 1972.

Mapping the 1972 modules onto what GSSK ships:

| Odum 1972 module | GSSK today |
|---|---|
| §II Energy Source | `source` |
| §V Heat Sink | `sink` |
| §VI Passive Energy Storage | `storage` |
| §VII Potential-Generating Storage (Fig. 1d — storage plus heat dispersal against backforce) | approximately the `consumer` **composite** |
| §IX Constant Gain Amplifier | `gain` |
| §X Work Gate | `interaction` |
| §XI Switch | `switch` |
| §XIII Cycling Receptor | `loop_limited` |
| §XV Economic Transactor | `exchange` |
| §XII Self-Maintaining | the `producer` **composite** |
| §XIV Production and Regeneration (P-R) | **absent** |

Two things fall out. What GSSK calls `producer` is Odum's *Self-Maintaining
Module* — a primitive in 1972, demoted here to a composite. And Odum's actual
producer, the green plant, is **P-R, which GSSK does not have at all.**

### The concrete thing that does not work

So: write P-R. Odum's definition is a cycling receptor, plus a self-maintaining
module which it feeds, plus a work feedback loop — in GSSK terms `loop_limited`
plus whatever carries the self-maintaining body.

It cannot be written. A composite template may only name fundamentals — not by
convention but by schema: `ArchetypeDefn.nodes.items.type` is
`{"$ref": "#/$defs/PrimitiveNodeType"}`, whose enum is closed over the nine
primitives. `producer` and `consumer` are both composites, so neither may appear
inside another template — and P-R is precisely a template that needs to.

This is the bug, and Odum's canonical example is what it blocks. Every composite
must be inlined to fundamentals by hand, so definitions cannot be built out of
each other, and the vocabulary cannot compound the way Odum's own explicitly
does.

### The same limit, at the other end of the scale

The identical constraint blocks hierarchical modelling. A 64-node subsystem
might be a farm; a farm is one component in a nation; a nation is one component
in a trading bloc — and downward, a prime minister's coffee decomposes to
molecules. Modellers expect to package a subsystem into one component, wire
packages together, and open one to see inside, with computation unaffected
throughout. Grasshopper clusters, Blender node groups, Houdini subnets and
Simulink virtual subsystems all afford this.

Four things in the current design prevent it:

1. **Expansion is one-way and boot-time only.** Composites expand at
   `GSSK_Init`; nothing composite survives into the solver. There is no repack,
   and no way to close a package again once open.
2. **Membership is one level deep.** `GSSK_GetNodeComposite` returns a name, not
   a path. A node can say it belongs to `producer_2`; it cannot say it belongs
   to `world/eu/uk/farm_3/producer_2`.
3. **`system_frame` — the one symbol meant for encapsulation — encapsulates
   nothing.** The schema records it as a single `constant` node and the docs
   call it a placeholder.
4. **Motif scanning is capped at 64 nodes** by
   `bool adj[GSSK_MOTIF_SCAN_NODE_LIMIT * GSSK_MOTIF_SCAN_NODE_LIMIT]` — a
   fixed stack array — plus an O(N³) triple loop for 3-node motifs. An
   implementation artifact that reads as a modelling limit.

### The distinction that has to be got right first

Two different operations were being conflated in the discussion that produced
this ADR, and conflating them makes the problem look mathematically hard when
it is not.

**Packaging is syntactic.** Group N nodes, name the group, open it to look
inside. The solver expands it and evaluates the contents. Nothing is replaced,
nothing is approximated, and the computed trajectory is bit-identical whether
the package is open or closed. There is **no mathematical obstruction of any
kind**: any subgraph, linear or not, any depth, with feedback, feedforward and
recursion.

**Surrogate folding is semantic.** Replace a subgraph with an effective
transfer function so its interior is never computed. Exact for a linear
subgraph — the boundary restriction of `exp(A·t)`, which the Padé (3,3) path
already computes — and approximate when interaction edges make the subgraph
bilinear, which is Giannantoni's n-et problem and open for n ≥ 3.

Simulink names this split: **virtual** versus **atomic** subsystems. Grasshopper
offers only the virtual kind, which is why it composes to any depth without
qualification.

Everything a modeller means by "package this and open it later" is the virtual
kind. Surrogate folding is an optional performance feature that can arrive
later or never.

## Decision

**1. Composite membership becomes a path.**

`GSSK_GetNodeComposite` returns `world/eu/uk/farm_3/producer_2` rather than
`producer_2`. Storage stays flat — nodes remain a single array, edges remain
flat references — and hierarchy is carried entirely in the path string. A level
is a path-prefix group-by.

Flat-plus-path is chosen over nested JSON for the reasons the flat design was
chosen originally: every operation (integration, mutation, ordinality,
serialisation) stays uniform over one kind of node; tooling never has to
recurse; cross-level edges stay directly expressible; and a node can carry
membership in overlapping patterns, which a tree cannot represent.

`GSSK_GetNodeRole` continues to return the leaf template id (`body`, `gate`).

Querying a level is the operation this exists to serve, so it gets named API
rather than leaving every caller to parse paths by hand:

- `GSSK_GetNodeCompositeDepth(inst, node_idx)` — number of path segments.
- `GSSK_GetNodeCompositeAt(inst, node_idx, depth)` — the segment at a depth, so
  a caller can group without string-splitting.
- `GSSK_CountNodesUnder(inst, prefix)` and
  `GSSK_GetNodeUnder(inst, prefix, i)` — iterate the members of a frame,
  transitively. A `NULL` or empty prefix means the root level.

These are read-only accessors over the existing flat node array; none of them
allocates, and none changes what the solver integrates.

**2. Archetype templates may reference other archetypes.**

`ArchetypeDefn.nodes.items.type` widens from `PrimitiveNodeType` to the same
union `Node.type` already accepts: a fundamental, a built-in composite, or a
declared archetype. A template node's type may then name a fundamental *or*
another archetype, expanded
recursively at `GSSK_Init`. Expansion detects cycles and enforces a depth
bound, rejecting a self-referential archetype at parse time with a named error
rather than recursing. Generated ids derive from the expansion path, so
uniqueness is structural rather than a prefixing convention.

This is what makes `producer = loop_limited + consumer` *expressible*. Whether
`producer` should be redefined that way is deliberately **not** decided here —
see Non-decisions.

**3. `system_frame` encapsulates a named node set.**

It becomes real virtual packaging: a frame names a set of members, contributes
no equations, and expands to exactly its contents. Closing and opening a frame
is a view operation with no numerical consequence. This is exact by
construction, so it carries no tolerance and no error estimate.

**4. Packaging is virtual. Surrogate folding is a separate, later feature.**

When it arrives it will **fold with a computed error bound**, not refuse. An
earlier draft of this ADR proposed refusing to fold any subgraph containing an
interaction edge. That was wrong twice over. It would have been useless in
practice — a `producer` contains an autocatalytic interaction gate, so the most
common Odum composite could never be folded — and it contradicts the precedent
this codebase already set in Phase 1, where interaction edges are handled by
Padé linearisation with `GSSK_GetEdgeErrorEstimate` and
`GSSK_GetStepErrorEstimate` reported every step.

The rule this repo already follows is: **approximate where you must, compute
the error, expose it, and let the caller set a tolerance.** What is forbidden is
*silent* approximation, not approximation.

**5. Motif scanning becomes level-local, with a heap adjacency.**

Scan within a frame rather than across the whole graph. The performance
argument is secondary; the modelling argument is that a motif spanning a
molecule and a nation is not a motif but a coincidence. Recurrence only carries
meaning inside a transformity band. With per-level scanning the node cap stops
binding in practice, and the fixed stack array is replaced by a heap allocation
so the remaining limit is O(N³) cost rather than 4 KB of stack.

**6. Ordinality is reported per level.**

A component may sit on a closed pathway at its own level and not at its
parent's. That is not an inconsistency to be resolved — it is what drawing a
system boundary means, and it is what `system_frame` is for. A single global
`ordinality = 1.0` with no scale attached is the incomplete part.

**7. `producer` keeps its name.**

It implements Odum 1972 §XII, the Self-Maintaining Module, and §XIV names P-R
as the green plant — so a literal reading would rename this to
`self_maintaining`. It is not renamed. `producer` is shipped, schema-visible,
covered by golden CSVs, and is the word modellers actually use for a
self-maintaining autotrophic unit; the 1972 module names are a physical
taxonomy, not a modelling vocabulary. The mapping is recorded in
`docs/odum_1972_conformance.md` so the correspondence is documented rather than
lost, and a future P-R must therefore take a different type name.

**8. Cross-scale edges are typed as control, not flow.**

Interaction happens within a transformity band. Coupling between bands is the
amplifier relationship — a small high-transformity signal modulating a large
low-transformity flow — which is Odum's `gain` triangle. An edge crossing a
frame boundary between levels is expected to be control-typed; an ordinary
flow edge crossing scales is a modelling error worth reporting.

## Why transformity is the scale coordinate

Levels are not arbitrary nesting depth. Odum's claim is that energy hierarchy
*is* spatial and temporal scale: many joules of sunlight, fewer of biomass,
fewer still of human work, and tiny flows of information and control. Moving up
in transformity is moving up in scale.

The kernel already computes and propagates transformity every step, so the
scale axis does not need inventing — a level is a transformity band, and the
existing `quality_input` machinery already populates it.

## Non-decisions

- **Whether to add P-R** as a built-in composite once nesting makes it
  expressible. It is the more faithful "green plant" and the obvious first user
  of archetype nesting, but adding a built-in is a vocabulary decision, and this
  ADR is about mechanism. Note that if P-R is added, it needs a name that is not
  `producer` — see the decision below.
- **Whether exogenous `source` and `sink` count toward ordinality.** They are
  boundary conditions rather than system components, which argues for excluding
  them; that interacts with per-level ordinality and is deferred.
- **The surrogate fold operator itself**, including which error norm bounds it.

## Consequences

**Schema.** `composite` becomes a path-valued string. Existing single-segment
values remain valid as depth-1 paths, so serialised models and the
`GSSK_GetNodeComposite` contract stay readable; consumers that string-compare
the whole value must switch to a prefix or leaf comparison.

**Recursive expansion is a new failure surface.** Cycles and unbounded depth
must be rejected at parse time with a named error, in the style of the existing
`GSSK_ERR_SCHEMA_VIOLATION` messages that name the offending node id.

**`GSSK_AddNode` keeps accepting only fundamentals.** Runtime emergence produces
fundamentals; naming stays a recognition layer above. This is unchanged and is
the right split — it is also the split the Giannantoni engine's generative step
already follows.

**The generative engine gains a scale.** `src/engine.c` currently reports one
global ordinality and spawns one untyped `regulator`. Per-level ordinality and
level-local motif recognition are what let an emergent component be *recognised*
as a producer or a recycler rather than named by a hardcoded string.

**Golden CSVs are unaffected** by (1)–(3) and (5)–(7), because none of them
changes what the solver integrates. Only a future surrogate fold would, and
that is out of scope here.

**Risk: paths encourage deep hierarchies that nothing validates.** A model may
declare `a/b/c/d/e/f` with no member at intermediate levels. A depth bound and
a well-formedness check on the path are required, not optional.
