# ADR 0016 — Divide and subtract are actions of the interaction module

- **Status**: proposed
- **Date**: 2026-09-16
- **Task**: `adr-0016-interaction-actions`
- **Supersedes**: nothing
- **Depends on**: [ADR 0008](0008-nary-interaction-and-subtracting-action.md) — the subtracting action;
  [ADR 0012](0012-where-a-law-lives.md) — module laws leave pathways;
  [ADR 0013](0013-module-roles-not-position.md) — roles, never position

## Context

ADR 0012 decision 5 moves seven laws off pathways and onto modules:
`interaction`, `limit`, `ratio`, `subtract`, `gain`, `threshold` and `exchange`.
Its migration note counts the models that use each one, including one each for
`ratio` and `subtract`. Once nothing writes the pathway forms, they are removed.

Five of the seven have a module to move to. The loader accepts a `module` block
on exactly these node kinds:

| pathway law | module |
|---|---|
| `interaction` | `interaction` |
| `gain` | `gain` |
| `threshold` | `switch` |
| `limit` | `loop_limited` |
| `exchange` | `exchange` |
| `ratio` | — |
| `subtract` | — |

Neither ADR 0012 nor ADR 0013 says how `ratio` and `subtract` are written as
modules. ADR 0013 names roles for interaction, gain, switch and loop_limited,
and the transactor's four legs; it names none for a divisor or a subtracting
action. So the migration, which must rewrite every model that uses one of the
seven, cannot finish, and removing the seven pathway laws would take two of
Odum's arithmetic actions out of the language altogether.

### What Odum draws

ADR 0008 already records the figure. *Modeling for All Scales* Fig. 2.6 shows
**one** interaction glyph computing several ways: (a) a product of two inputs,
(c) a product of three, (d) a divisor action, (e) a subtracting action. "The
overloading is deliberate: one shape keeps the diagram legible while the
arithmetic inside it varies."

`ratio` is (d), and `subtract` is (e). They are not separate symbols. They are
what the work gate computes.

### Why the pathway spellings are separate today

ADR 0008 rejected the obvious spelling, a `params.op` string on `interaction`,
for two reasons, both about GSSK's edge laws:

1. **Two spellings of one thing.** `ratio` had already shipped as its own logic
   type, so `op: "div"` would have been a second way to write it.
2. **An unchecked string.** `logic` is an enum, and `-Wall -Wextra -Werror`
   turns a missed case at any of its switch sites into a build failure. A string
   compared at each site is checked by nothing, and a site that forgot one would
   return 0.0 and look like a modelling result.

## Decision

1. **`ratio` and `subtract` become actions of the `interaction` module.** A work
   gate's `module` block takes an `action`:

   ```json
   { "id": "gate", "type": "interaction",
     "module": { "k": 0.3, "action": "subtract" } }
   ```

   | `action` | Odum Fig. 2.6 | flow | controls |
   |---|---|---|---|
   | `multiply` (the default) | (a), (c) | `F = k · Q_energy · Π Q_control` | any number, as today |
   | `divide` | (d) | `F = k · Q_energy / max(Q_control, ε)` | exactly one |
   | `subtract` | (e) | `F = max(0, k · (Q_energy − Q_control))` | exactly one |

   Absent, `action` is `multiply`, so every interaction module written so far
   means what it meant before. The roles are ADR 0013's: one `energy` input,
   which is consumed, and controls, which are read and never consumed. Output
   fans out to the outgoing pathways by weight, like every other module's.

2. **Both of ADR 0008's objections are answered, not set aside.**

   - *Two spellings.* That objection was to an `op` field existing **beside** a
     `ratio` logic type. ADR 0012 removes the pathway form, so once the migration
     is done `action` is the **only** spelling. No second one survives to drift.
   - *An unchecked string.* `action` is parsed once, at load, into an enum. An
     unknown value is a load error that names the three valid ones. Everything
     past the loader switches on the enum, so `-Werror` keeps the check ADR 0008
     wanted.

3. **The arity rules are enforced at load.** `divide` and `subtract` need exactly
   one control. A quotient or difference of three things is not defined by the
   figure, and ADR 0008 already declined to invent an associativity. Violating
   the rule is a named error, the same as a cycling receptor with a surplus input.

4. **Each action keeps the numerical character its law has now.**

   - `divide` is linear in `Q_energy`, with a conductance `k / Q_control` that
     changes with state. The flow matrix is therefore not constant and `ψ` is
     generally non-zero, the same as for `multiply`. The `ε` floor is ADR 0002's:
     it keeps the flow finite as the control goes to zero.
   - `subtract` is **linear** in the state while the difference is positive, so
     between crossings the matrix is constant and the solution is exact. The
     clamp at zero is the semantics, not a guard (ADR 0008: a barbed pathway
     cannot carry negative flow). The instant `Q_energy = Q_control` is located
     as an **event**, the same way a switch module's threshold crossing is, so
     the piecewise solution has no seam. This is stronger than GSSK, where the
     crossing is bounded by `dt`.

## Considered and not chosen

**New node kinds, `ratio` and `subtract`.** This needs no `action` field. But it
creates two symbols Odum does not draw, next to the one he does. A reader
translating a Fig. 2.6 diagram would have to know the glyph maps to three kinds
depending on what is written inside it. The module form should follow the
figure, and the figure has one glyph.

**Keeping `ratio` and `subtract` as pathway laws while the other five move.**
This is the option ADR 0012 decided against for all seven: two places a law can
live, and a reader who cannot tell from a pathway whether it is a module law.
Exempting the two laws that happen to have no module yet would reintroduce that
ambiguity, only for these two.

## What this does not decide

- **How the seeds are re-modelled.** `examples/giannantoni/input.json` writes its
  interactions through the `generative_production` and `ordinal_feedback`
  labels, which load as the `interaction` pathway law between integrating
  stocks. Rewriting it as a module is a modelling change, as ADR 0014 noted. The
  migration PR decides it and says what it chose.
- **Actions beyond Fig. 2.6.** `min`, `max` and other combinations are not
  added. The table is the figure.
- **The projection from GSSK.** GSSK keeps `ratio` and `subtract` as edge logic
  types. Translating them into interaction modules with an action is the
  projection's job, and the coverage report states it.

## Consequences

- The loader accepts `action` on an `interaction` module and rejects it on any
  other kind. An unknown value, or the wrong number of controls for `divide` or
  `subtract`, is a named load error.
- The module flow, the emergy pass's module flow, and event location each gain
  the two actions. Event location is needed for `subtract` only.
- The migration can then rewrite the one `ratio` model and the one `subtract`
  model. Until the pathway laws are removed, a test asserts that each module
  form reproduces its pathway form's trajectory. After the removal, the pathway
  form no longer loads, so that test asserts the closed form instead.
- ADR 0012's removal becomes possible for all seven laws, with nothing from Odum
  Fig. 2.6 left without a spelling.
