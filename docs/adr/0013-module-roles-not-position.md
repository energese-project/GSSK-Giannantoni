# ADR 0013 — A module's inputs are named by role, never by position

- **Status**: proposed
- **Date**: 2026-09-11
- **Task**: `module-role-schema`
- **Supersedes**: nothing
- **Depends on**: [ADR 0012](0012-where-a-law-lives.md) — module laws leave pathways

## Context

ADR 0012 decided that Odum's work gate, amplifier, switch, cycling receptor and
transactor are *modules* — hyperedges over their neighbourhood — and left the
concrete schema for declaring one open. The schema has to say which of a
module's inputs does what: in a work gate, one input supplies the energy and is
consumed, while the others are controls that are read and not depleted.

GSSK already answers that question, and how it answers it is the reason this
ADR exists.

### GSSK decides a module's inputs by the order of lines in a file

GSSK's `Edge` schema has no role field. Its fields are `id`, `origin`,
`target`, `carrier`, `logic`, `params`, `forcing`, `output_mode`,
`coupled_edge`, `active` and `visual` — nothing that says what a pathway *is*
to the processing node it enters. So each evaluator in `src/gssk.c` assigns
roles by scanning the edge array and counting:

| node | first incoming edge | second incoming edge | further edges |
|---|---|---|---|
| `interaction` | energy — **consumed** | control — read | controls |
| `gain` | control — read | energy — **consumed** | ignored |
| `switch` | switched flow — **consumed** | sensor — read | ignored |
| `loop_limited` | the input | **silently ignored** | ignored |
| `exchange` | legs found by **carrier**, not position | | |

Two things are wrong with this, and the second is worse.

The roles are **positional**, so they are decided by the order in which
pathways happen to be listed in a JSON document. JSON does preserve the order of
an array, so the defect is not that the order gets lost — it is that the order
is given meaning at all. Nobody reading a model treats the sequence of its
pathways as significant, and the tools that touch a model change it freely
without anyone intending to change the model: a formatter sorting by id, a
merge, a generator emitting edges in whatever order it walks its graph.

And the convention is **not even consistent across node types.** The first
incoming edge is the consumed input of a work gate and of a switch, and the
*read-only* input of an amplifier. There is no rule to learn, only a table to
memorise.

### Demonstrated, not inferred

Two reproductions, each a single model run twice with nothing changed but the
order of two lines in `edges`.

**Work gate** — `interaction`, `k = 0.001`, inputs `grass` (100) and `sun` (50),
output `cow`, RK4 to `t = 5`:

| order | drained | grass | sun | cow |
|---|---|---|---|---|
| grass first | grass | **77.88** | 50.00 | 22.12 |
| sun first | sun | 100.00 | **30.33** | 19.67 |

**Amplifier** — `gain`, `k = 0.05`, inputs `signal` (10) and `power` (100),
output `load`, RK4 to `t = 5`:

| order | read as control | drained | load |
|---|---|---|---|
| signal first | signal | power 100 → 97.50 | **2.50** |
| power first | power | signal 10 → **0.00** | **25.00** |

The amplifier case changes the answer **tenfold** and exhausts a different
stock. Neither run reports an error or a warning, and no schema check could
catch either, because there is no field to validate. This is the same shape as
every defect this repository has found in itself: output that is plausible,
specific, and wrong.

The models:

```json
{"metadata":{"schema_version":4,"name":"gate order"},
 "nodes":[{"id":"grass","type":"storage","value":100.0},
          {"id":"sun","type":"storage","value":50.0},
          {"id":"gate","type":"interaction","value":0.0,"params":{"k":0.001}},
          {"id":"cow","type":"storage","value":0.0}],
 "edges":[{"id":"e_grass","origin":"grass","target":"gate"},
          {"id":"e_sun","origin":"sun","target":"gate"},
          {"id":"out","origin":"gate","target":"cow"}],
 "config":{"t_start":0,"t_end":5,"dt":0.01,"method":"rk4"}}
```

and for the amplifier, the same with `signal`/`power`, `"type":"gain"` and
`"k":0.05`. Swap the first two entries of `edges` to reproduce.

## Decision

**1. Every pathway into a module declares a `role`. Position is never
consulted.**

A module reads the role off each pathway that enters it, and nothing else. The
order of pathways in a document carries no meaning, so reordering them cannot
change a result — not by convention, but because the code that could let it
does not exist.

**2. Two roles cover every processing module: `energy` and `control`.**

Odum draws the same distinction in every module of §IX–§XIII: one input supplies
the energy and is used up; the others signal, and are read but not depleted.
The role vocabulary follows that distinction rather than each module's
particular law:

- **`energy`** — consumed. The module drains it.
- **`control`** — read. The module's law depends on it; it is not depleted.

What differs between a work gate, an amplifier and a switch is the law that
combines those inputs, not the kinds of input it has. Keeping the vocabulary to
two names shared by every module is what makes it closed, which ADR 0012 set as
the condition for adding module laws at all.

**3. Each module type states the roles it requires, and a document that does not
supply exactly those is rejected by name.**

| module | `energy` inputs | `control` inputs | law |
|---|---|---|---|
| `interaction` (§X) | exactly 1 | 0 or more | `F = k · Q_energy · Π Q_control` |
| `gain` (§IX) | exactly 1 | exactly 1 | `F = k · Q_control`, drawn from energy |
| `switch` (§XI) | exactly 1 | exactly 1 | flow on while `Q_control` exceeds threshold |
| `loop_limited` (§XIII) | exactly 1 | none | `F = k · Q · C / (C + Q)` |

A missing role, a surplus input, or an input without a role is a load-time
error naming the module and the offending pathway. `loop_limited` in particular
stops silently discarding a second input: a second input is either given a role
the law uses, or it is an error.

**4. Fan-out is in proportion to each output pathway's weight.**

A module's output flow `F` is divided among the pathways leaving it in
proportion to their `weight`, normalised to sum to one. Equal weights give an
even split, which is GSSK's rule, so a model written with GSSK's assumption
behaves identically. Unequal weights give a declared split, without inventing a
new field: the pathway already carries a weight, and giving it this meaning on a
module's output is what a modeller would expect it to mean.

This quantity split is separate from the emergy split, and they compose without
conflict. `output_mode` still governs emergy at a bifurcation — `partition`
shares it in proportion to each branch's share of the flow, and that share is
exactly the one decided here.

**5. `exchange` becomes a module whose four legs are four roles.**

ADR 0012 noted that as a module, `exchange` could read its whole neighbourhood,
which is what makes leg discovery possible. Roles are how that is made
unambiguous: `goods_in`, `goods_out`, `counter_in`, `counter_out`. The carrier
work declined to *infer* payer and receiver because carrier identity alone could
not settle it; a role on each leg settles it by saying so. Exchange keeps its
own role set rather than sharing `energy`/`control`, because a transaction
couples two flows rather than consuming one and reading another.

## A finding for GSSK, not a decision for this repository

The positional convention, its inconsistency across node types, and the two
reproductions above are recorded here because this repository found them. They
are GSSK's to act on. The fix there is the same in kind — a role on the edge,
validated against what each node type requires — but whether and how to migrate
existing GSSK models is GSSK's decision, with its own users to consider.

## What this does not decide

- **The spelling of `role` in a document** — a field on the pathway, or keyed
  under `params`. Either satisfies decision 1; it is settled when the loader is
  written, and must be one spelling, not both.
- **How the projection supplies roles** when reading a GSSK model, which has
  none. A GSSK model's roles are recoverable only by applying GSSK's own
  positional table — the very thing rejected above. The projection can do that
  as an explicit, reported translation step, naming in the coverage output that
  roles were assigned by position because the source document gave no other
  information. Whether that is acceptable, or whether such modules should be
  reported as uncarriable, is left to the projection work.

## Consequences

**Reordering pathways can no longer change a result.** That is the property
this ADR exists for, and it should be tested directly — the same model loaded
with its pathways in two orders must produce byte-identical output.

**Every module in a migrated model states what each of its inputs is for.** That
is more to write, and it is also exactly the information a reader of an Odum
diagram gets from the drawing itself: the energy input and the control arrive
at different sides of the symbol.

**Risk: the projection's GSSK translation reintroduces position.** It is the one
place position has to be consulted, because the source has nothing else. Kept
honest by being explicit and reported, and by never being the path a native
model takes.
