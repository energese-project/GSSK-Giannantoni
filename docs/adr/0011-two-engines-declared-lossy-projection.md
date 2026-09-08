# ADR 0011 — The two engines do not converge; a declared-lossy projection instead

- **Status**: proposed
- **Date**: 2026-09-08
- **Task**: unassigned — no crux task yet
- **Supersedes**: nothing
- **Blocks**: nothing
- **Depends on**: [ADR 0010](0010-composites-nest-hierarchical-scale-by-path.md) for none of its decisions, but shares its source analysis

## Context

This repository now contains two simulation engines that share a build, a JSON
parser and nothing else.

- **`src/gssk.c`** — an Odum simulator. [`docs/odum_1972_conformance.md`](../odum_1972_conformance.md)
  scores it at 11 of Odum's 12 modules from the 1972 chapter, with §XIV (P-R)
  absent and §VII partial. It applies Giannantoni's IDC as its baseline solver
  wherever each module's law permits: exact for the linear core, exact as the
  Riccati duet for an isolated multiplicative junction, bounded-approximate for
  the cycling receptor, and piecewise across a switch.
- **`src/engine.c`** — a Maximum Ordinality engine. It implements persistence of
  form and the derivative drift ψ, the binary/duet/n-et branches, the harmony
  matrix over the (N−1) ordinal roots of unity, and ordinality-gated structural
  emergence. It does **not** implement Odum's language.

The obvious question is whether they should become one engine. The conformance
work says the two directions of that question have different answers, and only
one of them is an engineering problem.

### Direction one is settled and already shipped

Odum's module definitions feed Giannantoni's calculus. Every module's law is
either linear (exact under IDC), bilinear (exact as the duet, open as the n-et
for n ≥ 3), rational (bounded approximation), or discontinuous (piecewise
between events). That is a table with entries, not a blocker, and Phase 1 built
it.

### Direction two is not an engineering problem

Odum's language is **heterogeneous by construction**. Each of the dozen modules
carries its own mathematical definition, and the 1972 editor's description of
them as *"a mutually exclusive and exhaustive set"* is a claim that the typing is
load-bearing: a work gate is not a gain is not a switch, and none reduces to
another.

Giannantoni's MOP is **homogeneous**. Its state is Relational Space coordinates
and its structure is the N×N matrix of ordinal relationships `α_ij`, reduced by
the Harmony Relationships to one reference couple. Components enter as *relata*,
not as typed modules. The published framework does not address:

- a **switch** — ordinal coordinates are smooth, so a discontinuity has no
  direct form. It is handled piecewise, restarting at each event, which is what
  the kernel already does with Illinois crossing detection
- a **constant gain amplifier** — the split between a control signal and the
  separate power source that amplifies it is not given an ordinal form, though
  Odum's transformity hierarchy and Giannantoni's ordinality are plausibly the
  same axis, so this may be a gap in exposition rather than in the mathematics
- an **exchange / transactor** — two carriers coupled by price, against one
  relational space; a product space is not obviously forbidden, only undone
- the **conservative/non-conservative distinction between module kinds** — MOP
  asserts non-conservativeness globally, where Odum locates it precisely, in the
  multiplicative junction

**These are unaddressed, not proven impossible.** Only one item on the list is
genuinely closed off: the n-et for n ≥ 3, which Giannantoni's own papers leave
open. An earlier revision of this ADR called all four "no representation
exists", which overstated the case and is corrected here.

### Three ways to support both, not two

An earlier revision offered a false dichotomy — flatten Odum's typing into
ordinal relations, or extend MOP with typed relations. There is a third option,
and it is better than either:

- **(a) Flatten.** Express every Odum module as ordinal relations. Discards the
  distinctions Odum calls exhaustive, and the conformance score with them.
- **(b) Extend MOP with typed relations.** Genuine research; appears in neither
  the 2006 nor the 2023 paper.
- **(c) Make MOP a layer, not a rival simulator.** Let Odum's module laws supply
  the dynamics, and compute MOP's quantities — ordinality, the harmony
  relationships, the drift ψ, emergent quality — *over* the resulting system.

(c) is available now and requires nothing new from the mathematics. MOP's
claims are structural and analytical: they are about the relational organisation
of a system and how it originates new quality. **None of them requires MOP to own
the equations of motion.** The current `src/engine.c` owns them — assigning each
node a `φ` template instead of using Odum's law — and that was a design error on
this repository's part, not a limit of Giannantoni's framework.

The preference is **(c)**, with (b) as the long-term research direction that
would eventually subsume it. (a) is rejected: it destroys the thing being
conformed to.

### The concrete symptom, today

The gap is not abstract. The generative engine's own emitted vocabulary is
already incompatible with the kernel in two independent ways:

1. `src/engine.c` spawns nodes of type `"regulator"`. The schema's
   `PrimitiveNodeType` enum is
   `[storage, source, sink, constant, interaction, gain, loop_limited, exchange, switch]`,
   and the built-in composites are `producer`, `consumer`, `misc_box`,
   `system_frame`. **`regulator` is neither.** A `gia_generate` output cannot be
   loaded by `GSSK_Init`.
2. `GIA_NODE_CONSUMER` names a GSSK *composite*, and `GSSK_AddNode` accepts only
   fundamentals — it *"expands nothing at runtime and rejects composite and
   archetype names rather than mis-modelling them as single nodes."* So even the
   overlapping name would not survive runtime insertion.

Some of that divergence is gratuitous rather than principled, and can be closed
without touching the research question.

## Decision

**1. The two engines stay separate. Neither includes the other's header.**

`src/gssk.c` is the Odum conformance authority; `src/engine.c` is the MOP
research engine. `src/engine.c` continues not to include `gssk.h`, for the reason
given at the top of `include/engine.h`: the kernel carries state as a flat
`double *` of scalar storages, and Relational Space coordinates are complex.

**2. A one-way projection, GSSK model → MOP relational space, that declares what
it drops.**

The direction is forced. GSSK is the conformance authority and MOP cannot
express its vocabulary, so a MOP → GSSK direction would have to invent typing it
does not have. The projection therefore reads a GSSK model and produces a MOP
input, and for every module it cannot represent it emits a named finding rather
than a silent substitution — the same rule Phase 1 set for the solver, and ADR
0010 restated for folding: **approximate where you must, report what you did,
never do it silently.**

**3. The projection reports a coverage number, per model.**

For a given Odum model: what fraction of its nodes and edges have a MOP
representation? This turns the unquantified part of
[`docs/giannantoni_assessment.md`](../giannantoni_assessment.md) §5.3 — the Level 1
claim that "any system modelled through the MOP lens has an explicit solution" —
into a measurement against real models rather than an assertion. A model at 100%
coverage is one where the MOP claim is testable end to end; a model at 60% names
the four modules that stop it.

This is the deliverable that makes the disagreement between the two frameworks
*productive* instead of philosophical.

A caution on reading the number: a model scores below 100% because *that model*
uses modules the projection cannot yet carry — overwhelmingly the switch and the
multiplicative chain of length ≥ 3. It does not mean MOP is that fraction
correct. Most Odum models are expected to score high, and the report must name
the blocking modules rather than publish a bare percentage.

**4. The MOP engine's node vocabulary aligns to the GSSK primitive enum.**

`GIA_NODE_REGULATOR` is retired and `GIA_NODE_CONSUMER` stops naming a composite.
Emitted node types are drawn from `PrimitiveNodeType`, so a `gia_generate` output
is loadable by `GSSK_Init` and round-trips through the projection.

*Which* primitive an emergent component should be is deliberately not decided
here — that is the emergent-quality problem (MOP holds that the new component's
quality follows from ordinal structure, where the current code uses a hardcoded
string), and it is tracked in TODO Phase 10.3. This decision only requires that
whatever is emitted be a type that exists.

**5. Network coupling in the MOP engine is ordinary work and is not blocked.**

`src/engine.c` today computes each component's trajectory from its own `φ` alone;
deleting every edge from a model leaves the output CSV byte-identical. That is a
gap in this repository, not in MOP, and nothing in Giannantoni's framework
prevents Odum's flow laws from supplying the coupling. Requirements are in
[`docs/mop_network_coupling_requirements.md`](../mop_network_coupling_requirements.md).

Separating this from decision 6 matters: without it, "the engines cannot
converge" reads as a statement about the mathematics when most of the distance
is unbuilt code.

**6. Convergence into a single engine is not attempted, and the criteria that
would unblock it are recorded rather than left to taste.**

Revisit this ADR when *all four* hold:

- MOP has a representation for a **discontinuity**, so §XI switch is expressible
- MOP distinguishes a **control signal from the power source amplifying it**, so
  §IX gain is expressible
- MOP carries **more than one carrier**, so §XV exchange is expressible
- the **n-et for n ≥ 3** is solved, or the projection carries a stated error
  bound for the multiplicative junction

The first three are absences in the published framework, not in this code. The
fourth is open in Giannantoni's own papers.

## Consequences

**The trajectory comparison becomes possible, and bounded.** TODO Phase 10.4
carries a deferred bridge — run one model through both engines and diff the
trajectories, so the drift critique is measured against a real integrator rather
than only against a closed form. Decision 2 is what that bridge needs, and
decision 3 says what it may legitimately claim: the comparison is valid over the
projectable subset of a model and silent about the rest. Without the coverage
report, such a comparison would quietly compare a full Odum model against a
partial one and call the difference drift.

**`bench-giannantoni` keeps its current scope.** It compares mathematics against
mathematics with no graph machinery in the timing loop, which is the right way to
test the calculus claim. An end-to-end comparison against `bin/gssk` measures
infrastructure and answers a different question; it is not a replacement.

**Two engines is now a recorded position, not an accident.** It began as an
expedient — the fastest way to get MOP working without paying an integration
tax. Direction two says it is also correct, and this ADR is where that stops
being a coincidence.

**Risk: the projection becomes a de facto second model format.** It must remain a
read-only derivation from a GSSK model, with no authoring surface of its own.
If models start being written against the projection, the conformance authority
has silently moved.

**Risk: coverage is reported per model and read as a property of the framework.**
A 60% score means that model uses four modules MOP cannot express — not that MOP
is 60% correct. The report should name the modules, not just the fraction.

**No solver behaviour changes.** Nothing here alters a derivative, a logic type,
or an integration path in either engine. Golden CSVs are unaffected.
