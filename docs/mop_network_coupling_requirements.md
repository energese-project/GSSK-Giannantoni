# Requirements: network coupling in the MOP engine

## Why this exists

`src/engine.c` does not simulate a network. It computes each component's
trajectory from that component's own `φ` and nothing else:

```c
const gia_phi *p   = &m->nodes[i].phi;
double         idc = gia_idc_derivative(p, m->order, t);
```

The demonstration is blunt: **delete every edge from a model and the output CSV
is byte-identical.** Verified against `examples/giannantoni/input.json` — four
edges removed, `diff` reports no change.

So `flow_type` is parsed and written back out without ever being branched on;
edge `weight` reaches only the mean that sets an emergent edge's weight and the
modulus of the reference couple `α₁₂`; and there is no `F = k·Q`, no
`dQ/dt = ΣJ_in − ΣJ_out`, and no conservation of anything.

**This is a gap in this repository, not in Giannantoni's framework.** The `φ`
templates were chosen so the drift theorem could be checked exactly against
hand-derived closed forms — polynomial `φ` gives exact `φ⁽ⁿ⁾` — and that was the
right call for demonstrating persistence of form. It was the wrong thing to
leave standing as if it were a simulator. Nothing in MOP requires it: MOP's
claims are structural, and none of them requires MOP to own the equations of
motion (ADR 0011, option (c)).

## The design, which is not a fork

An earlier revision of this document posed a choice: let `engine.c` grow its own
coupling, or compute MOP over the kernel's stepped trajectories and recover
`φ = ln Q`. That was a badly posed question, and the second half of it was
actively wrong — recovering `φ` from a stepped trajectory reintroduces exactly
the integration error IDC exists to avoid, contaminating ψ with the drift of the
integrator used to measure it.

The correct treatment is already stated in both the source papers and this
repository's own assessment. For a network, the exponential form is a **matrix**,
not a scalar per node. Giannantoni 2023 gives Relational Space as

    {r}_s = e^{α(t)},    α an N×N matrix of ordinal coordinates

and `docs/giannantoni_assessment.md` §2 says the operational consequence
directly: *"the incipient solution reduces the simulation to a single matrix
exponentiation rather than thousands of Euler/RK4 steps."*

So:

1. Build the flow matrix **A** from Odum's edge laws — `F = k·Q_origin` for a
   barbed pathway, `F = k·Q_origin·Q_control` for the work gate, and so on. Each
   edge contributes to the entries for its endpoints, which is what makes a
   component's trajectory depend on its neighbours'.
2. The trajectory is `Q(t) = exp(A·t)·Q(0)`, evaluated by Padé (3,3). For
   constant **A** this is **exact and closed-form** — one exponentiation, no
   stepping, no accumulated error.
3. `α = A·t` *is* the matrix of ordinal coordinates, so the harmony
   relationships operate on the same object the dynamics do, rather than on a
   couple read off the first edge.

### What this does to the drift

It makes ψ mean something. For constant **A**, `d/dt exp(At) = A·exp(At)`, so the
traditional and incipient derivatives agree identically and ψ = 0 — the correct
answer for a constant-coefficient linear network, and the same result the scalar
theory gives for affine `φ`. ψ becomes non-zero exactly when **A** varies, with
`t` (forcing) or with `Q` (a work gate). The drift then reports a property of the
model rather than of a stipulated exponent.

The scalar `gia_phi` machinery is not discarded. It is the analytic core — the
only place where ψ is checked against hand-derived closed forms — and R5.3 keeps
those tests. The matrix layer is added beside it.

### What remains genuinely open

Only the case where **A** depends on **Q**: the multiplicative junction. Exact as
the duet, open as the n-et for n ≥ 3, Padé linearisation with a reported error
bound otherwise. That is Giannantoni's own open problem, not a choice this
repository has to make.

## Functional requirements

Numbered so they can be cited from tasks and tests.

### R1 — Edges participate in dynamics

**R1.1** A component's trajectory must depend on its inbound and outbound edges.
The falsifying test is the one above: removing an edge must change the CSV.

**R1.2** Flow must be computed per edge, from a named logic, not implied.

**R1.3** The node balance `dQ/dt = ΣJ_in − ΣJ_out` must hold for storage-like
components.

### R2 — Edge logic vocabulary

**R2.1** At minimum the logics that Odum 1972 defines mathematically:
`linear` (§III barbed, `F = k·Q_origin`), `interaction` (§X work gate,
`F = k·Q_origin·Q_control`), `reversible` (§III barb-less,
`F = k·(Q_a − Q_b)`), `limit` (§XIII cycling receptor, Michaelis-Menten) and
`threshold` (§XI switch).

**R2.2** `constant`, `ratio` and `subtract` to reach parity with the kernel.

**R2.3** `flow_type` in the seed format must either map onto a logic or be
rejected. Silently ignoring it, as today, is not acceptable — an inert
vocabulary that looks meaningful is worse than no vocabulary.

### R3 — Node vocabulary

**R3.1** Emitted and accepted node types must come from the schema's
`PrimitiveNodeType` enum, so a `gia_generate` output loads under `GSSK_Init`
(ADR 0011 decision 4). Today the engine emits `"regulator"`, which is neither a
primitive nor a built-in composite.

**R3.2** `sink`, `constant`, `gain`, `loop_limited`, `exchange` and `switch`
must be representable — the six primitives currently absent.

**R3.3** An unrecognised type must be rejected with a named error. Today it
becomes `GIA_NODE_UNKNOWN` and silently takes `φ = initial_value·t`.

### R4 — Conservation

**R4.1** A closed system's conserved carrier must balance to a stated tolerance,
and the residual must be reported, in the manner of
`GSSK_GetConservationError`.

**R4.2** Emergy is deliberately **not** conserved. Where emergy accounting
arrives it must be a separate pass with its own algebra, not folded into R4.1 —
that non-conservativeness is Giannantoni's motivation for IDC, and collapsing
the two accountings would erase the thing being demonstrated.

### R5 — psi must not be contaminated by the integrator

**R5.1** ψ is computed from a constructed `α`, never from a `φ` recovered by
taking logs of a stepped trajectory. There is no numerical integration in the
path that produces it.

**R5.2** For a constant flow matrix the engine must report ψ = 0 exactly, not
"small". A constant-coefficient linear network is a case where the two calculi
provably agree, and reporting a residue there would mean the implementation, not
the mathematics, is producing the drift.

**R5.3** The existing closed-form tests must keep passing unchanged. They are
the only *provable* results in the engine.

### R6 — What must not regress

**R6.1** Ordinality, cycle detection and the generative step.

**R6.2** The harmony matrix and both invariants.

**R6.3** The binary/duet/n-et branches.

**R6.4** Mode validation by structural diff, and the byte-identical seed
property on a functional run.

## Explicitly out of scope

- **Typed ordinal relations** — giving MOP a native form for the switch, the
  gain and the exchange. That is ADR 0011 decision 6, and research.
- **The n-et for n ≥ 3.** Open in Giannantoni's papers. Multiplicative chains
  get the kernel's Padé linearisation with a reported error bound, as today.
- **Emergy and transformity accounting.** Wanted (R4.2 reserves the shape) but a
  separate piece of work.
- **Hierarchy and composites.** ADR 0010.
