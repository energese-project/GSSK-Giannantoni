# ADR 0014 — Ordinality counts components on pathways that carry quantity

- **Status**: proposed
- **Date**: 2026-09-13
- **Task**: `adr-0014-ordinality-over-quantity`
- **Supersedes**: nothing
- **Depends on**: [ADR 0012](0012-where-a-law-lives.md) — where a law lives;
  [ADR 0013](0013-module-roles-not-position.md) — a module's inputs are named by role

## Context

Ordinality is the fraction of a system's components that lie on a closed
pathway. It is not a diagnostic in this engine: it **decides whether emergence
happens**. `gia_generate` scans for cycles, and if every component is on one it
returns the seed unchanged and the run is reported as functional. A second,
related verdict is closedness: `gia_system_is_closed` decides whether the
`--print` report shows a conservation residual that must be zero, or a net
boundary inflow that is expected.

ADR 0012 already recognised how much rides on that ratio. Its section *Emergence
is not desugaring* rejected having the loader invent gate nodes, because "a
bookkeeping node would count toward that ratio" and ordinality is "the one place
in the engine where bookkeeping must have no say."

That argument was made about **hidden** nodes. It applies just as much to the
explicit modules ADRs 0012 and 0013 then introduced, and nothing has applied it
to them yet.

### What changed underneath the cycle scan

Before module-hosted laws, a work gate's control was a *field* on a pathway —
`control_node` — and not an edge. Neither graph walk could see it, so it could
not close a loop and could not cross a boundary.

ADR 0013 made the control an **edge**, with `role: "control"`, and defined it as
*read and never consumed*. The emergy pass was brought into line in PR #14: a
control leg carries no flow. The two graph walks were not. `reaches()`, which
backs `gia_mark_cycles`, `gia_ordinality`, `gia_at_maximum_ordinality` and
`gia_generate`, follows every edge whatever its role, and `gia_system_is_closed`
inspects every edge's source whatever its role.

Separately, a module is counted as a component, even though PR #14 established
that it holds no carrier because it holds nothing: it is a hyperedge drawn as a
symbol, not a stock.

### Reproduction 1: a pure accumulator reported at maximum ordinality

`a` feeds a work gate `g` as its energy input; `g` delivers to both `a` and `b`;
`b` meters the gate as its control.

```
a --energy--> g --> a
b --control-> g --> b
```

`a` is on a real closed pathway. `b` is not: it receives from the gate, and
nothing it holds ever leaves. The trajectory confirms it — `b` only grows,
2.0 → 2.5347 over `t ∈ [0, 1]`, with conservation holding to 1.8e-15.

| | engine today | what the system is |
|---|---|---|
| ordinality | **1.0000** | 0.5 (`a` of `{a, b}`) |
| at maximum ordinality | **yes** | no |
| `gia_generate` | output **==** input — functional | an open pathway to close |

The loop `b → g → b` exists only through a leg that carries nothing, and it is
enough to tell the MOP step there is nothing left to originate.

### Reproduction 2: a conserving system reported open

A materially closed loop `a → g → b → a`. The only difference between the two
models is where the gate's control is read from.

| control read from | `gia_system_is_closed` | conservation residual at t=1 |
|---|---|---|
| stock `b` | yes | 1.8e-15 |
| constant `c` | **no** | 1.4e-14 |

The engine's own conservation measurement says nothing crossed the boundary,
while its closedness verdict says something did. The `--print` report then
explains the difference with a sentence that is false for this model:

```
  system              open -- a pathway leaves a held
                      component, which delivers quantity
                      without being depleted (Odum SecII)
  net boundary inflow 1.42109e-14  (expected, not an error)
```

The constant delivers no quantity at all. It is read.

### Reproduction 3: ordinality depends on where a law is drawn

One system, written twice. `a ⇄ b` is a closed loop. A source feeds a work gate
metered by `a`, which delivers into `a`. `a` and `b` are on a closed pathway;
the source is not. Written first with the gate as a pathway law, then as a
module. Three candidate rules were computed against both:

| rule | pathway spelling | module spelling |
|---|---|---|
| every edge walked, every node counted (**today**) | 0.6667 | **0.7500** |
| controls skipped, modules counted | 0.6667 | **0.5000** |
| controls skipped, modules passed through and not counted | 0.6667 | **0.6667** |

Today's rule over-counts: the control `a → g` closes a false loop through the
gate, so the gate is counted as on a cycle. Merely skipping controls
under-counts, because the gate is still in the denominator despite holding
nothing. Only the third rule gives the same answer for both spellings.

ADR 0012 decision 5 migrates every model from the first spelling to the second.
Under today's rule that migration would change the emergence verdict of a model
by how it is drawn, not by what it is.

## Decision

1. **A pathway is closed by legs that carry quantity.** A `control` leg is read
   and never consumed (ADR 0013), so it neither closes a pathway nor completes a
   cycle. Every other leg carries quantity and is walked: a pathway with no role,
   a module's `energy` input, its outputs, and a transactor's four legs.

2. **A module is passed through, not counted.** Quantity enters a module on its
   energy leg and leaves on its outputs, so a cycle may run *through* one and is
   no less closed for it. But a module holds nothing — the reason it names no
   carrier — so it is not a component. It is excluded from both the numerator
   and the denominator of ordinality, its own `on_cycle` flag stays false, and
   it is never chosen as the open component `gia_generate` closes.

3. **Closedness asks the same question at the boundary.** A system is open when
   *quantity* enters from a held component. A control read from a source or a
   constant moves nothing across the boundary and leaves the system closed. A
   pathway leaving a module does not cross the boundary either: a module only
   passes on what its energy leg brought, and that leg is where the question is
   asked.

4. **Invariance under re-spelling is the test.** The same system written with a
   law on a pathway and written as a module must have the same ordinality, the
   same closedness, and therefore the same generative verdict. That property is
   what licenses the migration in ADR 0012 decision 5, and it is asserted
   directly rather than inferred from the rule.

5. **One cycle scan.** `gia_generate` currently keeps its own inline copy of the
   ordinality computation beside `gia_ordinality`. There is one scan, used by
   both, so the reported ordinality and the one that decides emergence cannot
   disagree.

This completes ADR 0012's *Emergence is not desugaring* argument; it does not
reverse anything. 0012 kept bookkeeping from voting on emergence by refusing to
synthesise nodes. This keeps it from voting through the nodes and legs that were
then written explicitly.

## What this does not decide

- **An information or control ordinality.** A loop closed by controls is a real
  thing in Odum's diagrams — regulation, feedback of information — and it may
  deserve a measure. If it does, it is a *separate, named* measure. It is not
  folded into this one, whose job is to decide whether a system of flows has
  open relationships left to close.
- **Weighting by flow magnitude.** Ordinality stays a count of components. A
  pathway carrying a trickle closes a loop as fully as a torrent.
- **How `examples/giannantoni/input.json` is re-modelled.** Its `interaction_1`
  is written as an integrating stock with a pathway law on each side, and the
  consumer's feedback into it is *material*: the `ordinal_feedback` alias maps to
  the interaction law, which draws from the consumer. Rewriting it as a module
  with the consumer as a control is a different system, not a new spelling of
  the same one. The migration task decides which system the seed is meant to be
  and says so.
- **Whether a module keeps a CSV column.** A module is written today as a state
  column that is identically zero — `g_Q`, `g_Em` and the rest. That is the
  "state column nobody wrote" ADR 0012 cited against desugaring, and it follows
  from the same principle as decision 2, but it changes the output format rather
  than a verdict and is left to its own change.

## Consequences

- `gia_mark_cycles`, `gia_ordinality`, `gia_at_maximum_ordinality`,
  `gia_generate` and `gia_system_is_closed` change behaviour for every model
  containing a module. Models with no module and no role are unaffected: every
  leg carries quantity and every node is a component, which is today's rule.
- An existing test whose expected ordinality or closedness changes is updated
  with a comment saying which decision above moved it.
- The `--print` report stops printing one `components` count that includes
  modules, and reports components and modules separately.
- `migrate-models-to-modules` depends on this. Without it, migrating a model
  could change its generative verdict by spelling alone, and the migration's own
  acceptance test — that the two spellings agree — would fail for a reason that
  has nothing to do with the migration.
