# GSSK Schema v5: the General Systems Graph

- **Status**: draft, proposed. Nothing here is implemented.
- **Date**: 2026-10-05
- **Replaces, once accepted**: [`SPECIFICATION.md`](SPECIFICATION.md) (v2.0),
  [`gssk-schema.md`](gssk-schema.md) (v2 prose) and `gssk.schema.json` v4, and
  psde's `topology.json` v0.1 (`energese-project/psde`, `POC_SPECIFICATION.md` §2).
- **Builds on**: [ADR 0006](adr/0006-forcing-one-vocabulary-two-attachments.md) (forcing),
  [ADR 0012](adr/0012-where-a-law-lives.md) (a law lives where Odum draws it),
  [ADR 0013](adr/0013-module-roles-not-position.md) (roles, not position),
  [ADR 0016](adr/0016-interaction-actions.md) (interaction actions).

## 0. Why there is one schema

GSSK is the *General* Systems Simulation Kernel, so its model format should be the
most general description of a system it can execute. Until now there have been two
formats for the same thing: GSSK's JSON, where a flow's law is one of a closed list of
`logic` names, and psde's `topology.json`, where a law is a graph of arithmetic
operators. Each could express things the other could not. This schema is the
replacement for both. GSSK executes it, psde extracts into it and emits from it, and
there is no third format.

The test of "general" is the published record. Odum's BASIC listings were written for
real models, and they use the following:

- conditionals (`IF D > 30 THEN IV = 0`, `IF A < .001 THEN A = .001`, `IF X = 0 GOTO 348`)
- tabulated forcing (`READ`/`DATA`)
- randomness (`RND`)
- auxiliary variables, integer powers (`N*N`) and division (`M/P3`)
- reruns under changed conditions (`CONT`)

A schema that refuses these cannot hold the archive. Each one has a construct here
(§6–§9). Some are recognisably the same "IF", but they are different things, and each
gets its own construct rather than one general `if`.

### 0.1 Principles

1. **No mathematical strings.** A law is a tree of typed operator objects, never an
   expression to parse. (psde's IR rule, now the kernel's.)
2. **A law lives where Odum draws it** (ADR 0012). Pathways carry Odum §III laws only.
   Every other law is hosted by a module component.
3. **Inputs are named by role, never by position** (ADR 0013). The same applies inside a
   law. Every non-commutative operator names its operands (`numerator`, `minuend`,
   `condition`), so reordering a document can never change a result.
4. **No desugaring** (ADR 0012, "Emergence is not desugaring"). The loader never
   invents components. A law's operator tree lives *inside* the module that hosts it.
   It is not added to the system graph, so it cannot vote on ordinality and never
   appears as a state column.
5. **One spelling per thing.** A general law that is identical to a named module's law
   is rejected with a message naming the module to use (§5.4).
6. **Load or refuse by name.** Anything the kernel cannot execute faithfully is a load
   error naming the node, the field and the reason. Nothing falls back silently.

## 1. Document structure

```json
{
  "metadata":   { "schema_version": 5, "name": "…", "source": { … } },
  "carriers":   [ … ],
  "parameters": { "K1": 0.01, "IV": 0.3 },
  "signals":    [ … ],
  "nodes":      [ … ],
  "pathways":   [ … ],
  "events":     [ … ],
  "config":     { … },
  "scenarios":  [ … ],
  "archetypes": { … },
  "snapshot":   { … },
  "mutation_log": [ … ]
}
```

| Block | New in v5 | Purpose |
|---|---|---|
| `metadata` | `source` | v4 fields, plus optional provenance: citation key, locator, fidelity (`verbatim`, `corrected`, `adapted`, `original`), as the Odum archive's sidecars already record |
| `carriers` | — | unchanged from v4 |
| `parameters` | **yes** | named numbers. A law refers to them by name, and events and scenarios change them |
| `signals` | **yes** | named computed values that are not components (§5.5) |
| `nodes` | changed | components: storages, sources, sinks, modules (§2, §3) |
| `pathways` | renamed from `edges` | Odum §III pathways, with roles into modules (§4) |
| `events` | **yes** | latched discontinuities (§6.3) |
| `config` | extended | `seed`, plus v4 fields (§10) |
| `scenarios` | **yes** | reruns under changed parameters (§6.4) |
| `archetypes`, `snapshot`, `mutation_log` | — | unchanged from v4, except §11 |

`schema_version: 5` is required. v4 documents load through the migration in §12.1 or
are refused. They are never interpreted as v5.

## 2. Components

A component is a node in the system graph. Components are what ordinality counts
(ADR 0014), what carry quality, and what may appear as columns in the state vector.

| `type` | State? | Notes |
|---|---|---|
| `storage` | integrates | `value` = initial Q. Optional `bounds` (§6.2) |
| `source` | held | `value` is a number, or `forcing` (§7), or `law` (§5): the held value may be computed. This covers Odum's flow-limited source, whose remainder `R = J/(1 + Σ kQ)` is a law over the stores it feeds |
| `sink` | accumulates outflow | unchanged |
| `interaction`, `gain`, `switch`, `loop_limited`, `exchange` | none | the named modules (§3) |
| `process` | none | the general module: a law written as an operator tree (§3.2) |
| composites and archetypes | — | unchanged; they expand to the types above at init |

`constant` (v4) is removed. A fixed number with no flow belongs in `parameters`. A
fixed number that supplies flow is a `source`.

## 3. Modules

A module is a hyperedge over its neighbourhood (ADR 0012, decision 2):

- it reads every pathway that enters it, by role
- it computes one flow `F` from its law
- it drains `F` from its `energy` input
- it divides `F` over its outgoing pathways in proportion to their `weight` (ADR 0013,
  decision 4)

### 3.1 Named modules

These are Odum's 1972 modules, unchanged from ADRs 0013 and 0016. Each is *defined* by
the law tree given beside it. That definition is normative: it is what the module
means, and a kernel fast path must agree with it to `1e-12` relative.

| Module | Roles | Law, as a v5 tree (§5) |
|---|---|---|
| `interaction`, `action: multiply` | 1 energy, ≥0 control | `product(k, energy, control…)` |
| `interaction`, `action: divide` | 1 energy, 1 control | `divide(numerator: product(k, energy), denominator: max(control, ε))` |
| `interaction`, `action: subtract` | 1 energy, 1 control | `max(0, product(k, subtract(minuend: energy, subtrahend: control)))` |
| `gain` | 1 energy, 1 control | `product(k, control)` |
| `switch` | 1 energy, 1 control | `select(condition: compare(gt, control, threshold), then: k, else: 0)`, crossing located as an event |
| `loop_limited` | 1 energy | `divide(numerator: product(k, energy, C), denominator: sum(C, energy))` |
| `exchange` | `goods_in`, `goods_out`, `counter_in`, `counter_out` | ADR 0013 decision 5 |

Module parameters (`k`, `C`, `threshold`, `price`) are each either a number or a
reference `{"param": "K4"}`. A referenced parameter is what lets an event or a scenario
change a module.

### 3.2 The `process` module

`process` is the general module: Odum's rectangle with the function written inside
it. Its `law` is any tree from §5. It has exactly one `energy` input, which it drains,
and any number of `control` inputs, which it reads.

```json
{ "id": "residual_use", "type": "process",
  "law": { "op": "product",
           "factors": [ { "param": "K2" }, { "signal": "R" },
                        { "role": "energy" }, { "read": "N" } ] } }
```

`process` is what makes the schema general. It is also what ADR 0016 declined to add
under `interaction` ("the table is the figure"). That decision stands: `interaction`
stays exactly Fig. 2.6, and anything beyond it is a `process`, which is visibly a
different glyph. A reader can always tell whether a law is one of Odum's modules or
something written for this model.

## 4. Pathways

```json
{ "id": "p7", "origin": "D", "target": "drain_D", "law": "linear",
  "role": "energy", "carrier": "money", "weight": 1.0 }
```

| Field | Notes |
|---|---|
| `law` | `constant`, `linear` or `reversible`: Odum §III and nothing else (ADR 0012, decision 5). A pathway *into* a module has no `law`, because the module computes the flow |
| `k` | number or `{"param": …}`. Required for `constant`, `linear` and `reversible` |
| `role` | required on every pathway into a module: `energy` or `control`, or an exchange role. Forbidden elsewhere |
| `weight` | share of a module's output (ADR 0013, decision 4). Default 1 |
| `carrier`, `output_mode` | unchanged from v4 |
| `forcing` | drives `k` (ADR 0006), unchanged |

A `control` pathway carries information, not quantity. It drains nothing and is not a
quantity leg for ordinality (ADR 0014). This is psde's `state_vector` edge and Odum's
signal line.

`coupled_edge` is removed. Its job (pairing a goods flow with its money counter-flow)
is done by the `exchange` module's four roles.

## 5. Laws: the operator tree

A law is a JSON object. Leaves are values. Interior nodes are operators, and each
operator names its operands by role.

### 5.1 Leaves

| Leaf | Value |
|---|---|
| `{"const": 2.0}` | a literal number |
| `{"param": "K1"}` | a parameter's current value |
| `{"read": "N"}` | the current value of component `N` |
| `{"role": "energy"}` | the value at the module's energy input (inside a module only) |
| `{"signal": "R"}` | a signal's current value (§5.5) |
| `{"op": "time"}` | the solver's current time, at stage times |

**Reads must be drawn.** A module's law may `read` only components connected to it by
a pathway. The loader checks every `read` against the module's inputs and rejects an
undrawn dependency by name. This keeps every influence visible in the diagram, as
Odum's notation requires. Signals (§5.5) are the one exception, and they are reported.

### 5.2 Operators

Commutative operators take an array. Non-commutative operators take named operands.

| `op` | Operands | Value | Notes |
|---|---|---|---|
| `sum` | `terms: [≥2]` | Σ | |
| `product` | `factors: [≥2]` | Π | |
| `subtract` | `minuend`, `subtrahend` | a − b | replaces psde's `polarity` |
| `divide` | `numerator`, `denominator` | a / b | denominator 0 is `ERR_DIVERGENCE`, not ε. Use `max` to floor it |
| `power` | `base`, `exponent` (integer literal) | bⁿ | integer powers only. Real exponents use `exp`/`log` |
| `min`, `max` | `of: [≥2]` | | clamps |
| `negate`, `abs`, `floor` | `of` | | `floor` is BASIC `INT` |
| `exp`, `log`, `sqrt`, `sin`, `cos` | `of` | | out of domain is `ERR_DIVERGENCE` |
| `compare` | `relation` (`gt`, `ge`, `lt`, `le`), `left`, `right` | 1 or 0 | relation is an enum, not a string to parse |
| `select` | `condition`, `then`, `else` | | memoryless conditional (§6.1) |
| `previous` | `of`, `initial` | value of `of` at the last accepted step | §9 |
| `time` | — | t | |

`eq` is deliberately missing from `compare`. Floating-point equality between states is
a modelling error. BASIC's `IF X = 0` compares a *flag* parameter, and a flag is a
scenario (§6.4), not a condition inside a law.

### 5.3 Evaluation

- A law is evaluated at every solver stage.
- `previous` is latched per accepted step (§9).
- Evaluation is pure: a law has no side effects. Only events change parameters or
  state.

### 5.4 One spelling

At load, each `process` law is normalised and compared with the named modules'
defining trees (§3.1). An exact match is rejected:

> `process 'gate': this law is interaction(action: multiply); write it as that module`

This is ADR 0008's objection to two spellings, enforced by the loader rather than by
review.

### 5.5 Signals

```json
"signals": [
  { "id": "R",
    "law": { "op": "divide",
             "numerator": { "param": "I0" },
             "denominator": { "op": "sum", "terms": [
                 { "const": 1 },
                 { "op": "product", "factors": [ { "param": "K0" }, { "read": "N" }, { "read": "A" } ] },
                 { "param": "K1" } ] } } }
]
```

A signal is a named, computed value, such as a BASIC auxiliary (`R = I0/(1 + K0*N*A + K1)`).

- **It is not a component.** It holds no state, carries no quality and does not count
  for ordinality.
- **Ordering:** signals may read components, parameters and other signals. The loader
  orders them topologically and rejects a cycle by naming it. An algebraic loop is a
  modelling question for the author, not something the solver guesses at.
- **Outputs:** signals may be requested as output columns (§10).

Signals are the deliberate exception to "reads must be drawn". A listing's auxiliaries
often have no glyph, and refusing them would make most of the archive inexpressible.
So the kernel's dependency report lists every read through a signal, and the
undrawn dependency stays visible as a finding rather than being hidden.

Where Odum *does* draw the auxiliary (the remainder at a flow-limited source), it
belongs on the `source` as its `law` (§2), not as a signal.

## 6. Discontinuities: the four "IF"s

These are the conditionals in Odum's listings, by what each one does:

| Listing | What it is | v5 construct |
|---|---|---|
| flow on only while a store exceeds a level | memoryless switch | `switch` module, or `select` in a law (§6.1) |
| `IF A < .001 THEN A = .001` | floor on a stock | storage `bounds` (§6.2) |
| `IF D > 30 THEN IV = 0` | one-way change of a parameter | `event` (§6.3) |
| `IF X = 0 GOTO 348`, then `CONT` with `X = 1` | the same model under two conditions | `scenarios` (§6.4) |
| `IF T/T0 < 320 GOTO 200` | end of run | `config.t_end` |

These are deliberately not one construct. `IF D > 30 THEN IV = 0` is not a switch:
nothing ever sets `IV` back, so after D crosses 30, the investment stays off even if D
falls again. A switch would turn it back on, and the model would be wrong without any
error being raised.

### 6.1 Memoryless conditionals

`select(condition: compare(…), then: …, else: …)`. The value depends only on the
current state. Every `compare` in a law defines a crossing surface, which the solver
locates as an event (Illinois, as the switch module does now), so the piecewise
solution has no seam at the step boundary.

### 6.2 Storage bounds

```json
{ "id": "A", "type": "storage", "value": 0.1, "bounds": { "min": 0.001 } }
```

- **Applied:** after each accepted step, the state is projected into `[min, max]`.
- **Default:** `min: 0`, which is GSSK's current non-negativity clamp, now declared
  rather than global. `"min": null` allows negative values (a debt).
- **Ledger:** a projection creates or destroys quantity, so the kernel records it per
  storage as a **bound ledger**: the cumulative quantity added or removed by the bound.
  The conservation check reports this ledger instead of treating it as an error. The
  quality pass treats it as a boundary flow of *undeclared* transformity and says so,
  rather than assigning it one.

A bound is not a law. It does not change `dQ/dt`, only the state after the step. This
matches the BASIC exactly: the listing integrates and then clamps.

### 6.3 Events

```json
{ "id": "investment_stops",
  "when": { "op": "compare", "relation": "gt",
            "left": { "read": "D" }, "right": { "const": 30 } },
  "direction": "rising",
  "repeat": false,
  "actions": [ { "set_parameter": "IV", "to": { "const": 0 } } ] }
```

| Field | Meaning |
|---|---|
| `when` | a `compare` tree (§5.2) |
| `direction` | `rising` (false to true), `falling` (true to false), or `either` |
| `repeat` | `false`: the event fires once, then is spent (BASIC's latch). `true`: it fires on every crossing |
| `actions` | `set_parameter` (to a law), or `set_state` (a storage to a law, recorded in the bound ledger §6.2 because it creates quantity) |

Crossings are located within the step, and the actions take effect at the located
time. Events are evaluated in document order when several cross in one step.
`snapshot` records which events are spent, so a restored run resumes correctly.

**BASIC detail.** The listing tests its condition once per step, after the update. So
its event fires up to one `dt` later than the located crossing. The extractor records
this, and the round-trip tolerance for event listings is widened by one step's change.

### 6.4 Scenarios

```json
"scenarios": [
  { "id": "unlimited",  "parameters": { "X": 0 } },
  { "id": "renewable",  "parameters": { "X": 1 } }
]
```

A scenario is a complete run with some parameters overridden.

- **Output:** each scenario's output is labelled with its id.
- **Start state:** `start` is `initial` (the default: the document's initial values, as
  MACROEC's `CONT` path does by re-running lines 100–160) or `{"from": "<scenario id>"}`
  (continue from that scenario's final state).
- **Not a model change:** scenarios change no structure. A scenario that would need a
  different graph is a different model.

Inside a law, a scenario flag is read as `{"param": "X"}`. MACROEC's two production
functions are therefore one `process` law with
`select(condition: compare(lt, X, 0.5), then: …, else: …)`. Because a parameter is
constant within a run, this `compare` never crosses, and the solver knows it.

## 7. Forcing

ADR 0006's vocabulary is unchanged, with two additions: `random` (§8), and `table`,
for `READ`/`DATA` series and observed driving data.

```json
"forcing": { "waveform": "table",
             "times":  [0, 1, 2, 3],
             "values": [4.1, 5.0, 6.2, 5.8],
             "interpolation": "step",
             "extrapolation": "cycle" }
```

| Field | Values |
|---|---|
| `times` | strictly increasing |
| `values` | same length as `times` |
| `interpolation` | `step` (hold the last value: what a `READ` once per pass does) or `linear` |
| `extrapolation` | `hold` (keep the last value) or `cycle` (repeat with period `times[-1] − times[0]`: a seasonal year of monthly DATA) |
| `min`, `max` | as for every waveform |

Table forcing attaches to a source's value or to a pathway's `k`, as every waveform
does, and nowhere else. A data series comes from outside the system, as noise does
(§8), so a law that needs it reads a source over a `control` pathway.

## 8. Randomness

Randomness is a **source**, not an operator. In Odum's terms, noise reaches a system
from outside it, so it enters the diagram the way any outside driver does: as a source
whose held value is a forcing waveform. A law that needs it reads the source over a
`control` pathway, so the noise input is drawn like every other influence.

```json
{ "id": "weather", "type": "source", "value": 0,
  "forcing": { "waveform": "random", "distribution": "uniform", "low": 0, "high": 1 } }
```

- **`random` waveform:** `distribution` is `uniform` (`low`, `high`) or `normal`
  (`mean`, `sd`). `min`/`max` clamp it, as they do every waveform.
- **`jitter`** is unchanged, and is now defined as `random` `uniform` on
  `[mean − amplitude, mean + amplitude)`.
- **Draws:** from the instance PRNG (SplitMix64, as now). Each random source has its
  own stream position.
- **Latching:** each source draws **once per accepted step** and holds the value across
  that step's solver stages. A draw at every stage would make RK4 integrate a
  different function at each stage, and the result would depend on the solver.
- **Seed:** `config.seed` (a hex string, the snapshot's `rng_state` convention) seeds
  the run. Absent, `GSSK_DEFAULT_SEED` is used. Scenarios each start from the seed, so
  a scenario's result does not depend on which scenarios ran before it.
- **Solver:** a model with a random source runs fixed-step at `config.dt`. The step
  *is* the noise's correlation time, so changing `dt` changes the model, and the kernel
  says so (§10).
- **On a pathway:** a random waveform on a pathway's `k` is the same thing attached to
  a flow rather than a force (ADR 0006).

**Extraction.** Each `RND` occurrence in a listing becomes its own random source, with
a `control` pathway to whatever reads it. Two `RND` calls in one expression are two
independent draws in BASIC, so they are two sources here.

**Fidelity.** A BASIC `RND` stream cannot be reproduced bit for bit: each dialect has
its own generator. So a stochastic listing and its extraction are compared as
**ensembles** (the mean and quantiles of N runs of each), never trajectory by
trajectory. A trajectory match would only be evidence of having copied a PRNG.

## 9. Discrete-time listings

Many listings update variables in sequence within one pass. In MACROEC, for example:

```basic
367 I = S + DI * DT
368 Y = C + I
370 S = KS * Y
380 C = Y - S
```

Here `I` reads the `S` from the previous pass. The listing is then a difference
equation whose behaviour depends on `DT`, not an ODE that happens to be integrated by
Euler. `previous(of, initial)` expresses this exactly: the value of `of` at the last
accepted step.

A model that uses `previous`:

- is **step-dependent**. It runs fixed-step Euler at `config.dt` only, `method` must be
  `euler` or absent, and the dual-solver check is off with the reason reported (§10).
- is extracted this way only when the source's fidelity is `verbatim`. Otherwise the
  extractor writes the continuous model the listing approximates, and records
  `update_order: sequential` in its report (psde Module 0 §4.3).

Both are legitimate models, and they are different models. This schema keeps them
apart rather than pretending one is the other.

## 10. Config and the solver

`config` keeps v4's fields and adds:

| Field | Meaning |
|---|---|
| `seed` | §8 |
| `outputs` | ids of components and signals to report. Default: every storage |

**Method selection.** The solver is chosen by **analysing the laws**, not by the names
of the modules:

| What the laws contain | Integration | Dual check |
|---|---|---|
| only `constant`, `linear` and `reversible` pathways, and linear trees | IDC, exact | IDC vs RK4 |
| bilinear products of stores, isolated | IDC duet | IDC vs RK4 |
| rational (`divide`, `loop_limited`) | RK4 / DOPRI5 | bounded-approximate, reported |
| `compare`, `select`, events, `subtract`/`switch` modules | piecewise between located events | per segment |
| `exp`, `log`, `sin`, `cos`, `power` (non-polynomial) | RK4 / DOPRI5 | RK4 vs DOPRI5 |
| a random source, or `previous` | fixed-step at `dt` | **off**, reported as `GSSK_CONFIDENCE_STEP_DEPENDENT` |

Two consequences:

- A hand-written `process` law that happens to be linear gets the exact solver, just as
  a named module would.
- `GSSK_SolverConfidence` gains `STEP_DEPENDENT`, because "high" and "degraded" would
  both be false statements about a model whose answer is defined by its step.

## 11. Runtime mutation

| Call | v5 |
|---|---|
| `GSSK_AddNode` | accepts every §2 type except composites (unchanged rule), including `process` with its law |
| `GSSK_AddEdge` | renamed `GSSK_AddPathway`. Roles are validated against the target module at the call |
| `GSSK_SetEdgeK` | generalised to `GSSK_SetParameter(inst, name, value)`. A literal `k` on a pathway is set through the pathway's id as before |
| `GSSK_Deactivate*` | unchanged |
| events, signals | not mutable at runtime in v5.0. Adding one reorders evaluation, which is a structural change |

`mutation_log` gains `add_pathway` and `set_parameter`.

## 12. Migration

### 12.1 From GSSK v4

| v4 | v5 |
|---|---|
| node `storage`, `source`, `sink` | same |
| node `constant` | a `parameters` entry, or a `source` if a pathway draws from it |
| node `interaction`, `gain`, `loop_limited`, `switch`, `exchange` | same module, with roles on its input pathways. Positional roles are assigned by GSSK v4's own table and the migration report says so (ADR 0013, "What this does not decide") |
| edge logic `constant`, `linear`, `reversible` | pathway `law` |
| edge logic `interaction` | `interaction` module, `action: multiply` |
| edge logic `ratio` | `interaction` module, `action: divide` (ADR 0016) |
| edge logic `subtract` | `interaction` module, `action: subtract` (ADR 0016) |
| edge logic `limit` | `process`: `divide(product(k, energy), sum(1, divide(energy, C)))` |
| edge logic `threshold` | `switch` module |
| `params.control_node(s)` | `control` pathways into the module |
| `coupled_edge` | `exchange` module |
| global non-negative clamp | `bounds.min: 0` on every storage (the default) |
| `rng_state` seed | `config.seed` |

The migration creates modules, which are components. ADR 0012 forbids doing this
silently at load. So it is an explicit tool run (`gssk migrate v4 v5`), whose output
is a new document and a report of every module it introduced. The kernel itself never
loads v4.

### 12.2 From psde `topology.json` v0.1

| psde | v5 |
|---|---|
| `kind: stock` | `storage` |
| `kind: parameter` | `parameters` entry, or `source` when it feeds a flow directly |
| `kind: sink` | `sink` |
| `kind: operator` with a donor | a module: a named module if its tree matches one (§5.4), else `process` |
| operator without a donor | the tree of whichever law reads it; shared subtrees become signals |
| `mass_flow`, `entropy_loss` | pathways (`entropy_loss` targets a sink) |
| `state_vector`, `parameter_vector` | `control` pathways, or leaves in the law |
| `polarity` | operand roles (`minuend`/`subtrahend`, `numerator`/`denominator`) |

psde's symbolic derivation (`symbolic_solver.py`) reads v5 directly: a law tree maps
one-to-one onto a SymPy expression. psde's IR rule (no mathematical strings) is
satisfied unchanged.

### 12.3 From BASIC listings (the extractor's contract)

| Listing | v5 |
|---|---|
| assignment before the loop, constant | `parameters` |
| `X = X + (expr) * DT` in the loop | `storage` X, with its rate split into pathways and modules by psde Module 0 §3 |
| auxiliary in the loop | `signal`, or the `law` of a drawn source |
| `IF stock < c THEN stock = c` | `bounds` |
| `IF cond THEN param = v`, never reset | `event`, `repeat: false` |
| flag tested with `GOTO`, rerun with `CONT` | `scenarios` |
| `READ` inside the loop over `DATA` | `table` forcing, `interpolation: step` |
| `RND` | a `source` with `random` forcing, `uniform`, low 0, high 1, read by `control` pathway (§8) |
| `INT`, `ABS`, `EXP`, `LOG`, `SQR`, `SIN`, `COS` | the matching operator |
| sequential update, verbatim fidelity | `previous` |
| `PRINT`, `PSET`, `LINE`, `SCREEN`, `COLOR`, `CLS` | `config.outputs`. Graphics are dropped and listed in the report |
| `DIM` arrays | refused in v5.0 (§14) |

## 13. Worked example: DEVELOP

This fragment is from H.T. Odum's DEVELOP listing (lines 155–285; the run ends at line 400, `IF T / T0 < 320 GOTO 200`). It shows four
constructs together: the residual `R` as a signal, a floor on `A`, the latched
investment cut-off, and an interaction drawn from a drawn store.

```basic
155 A = .1
170 M = .01
180 D = 0
200 R = I0 / ( 1 + K0 * N*A + K1 )
210 DD = IV -N5*D*M
250 A = A + DA * DT
255 IF A <.001 THEN A= .001
280 D = D +DD*DT
285 IF D >30 THEN IV = 0
```

```json
{
  "metadata": { "schema_version": 5, "name": "DEVELOP (fragment)",
                "source": { "citation": "odum_develop", "fidelity": "adapted" } },
  "parameters": { "I0": 1, "K0": 2.3, "K1": 3, "N5": 0.01, "IV": 0.3 },
  "signals": [
    { "id": "R",
      "law": { "op": "divide",
               "numerator": { "param": "I0" },
               "denominator": { "op": "sum", "terms": [
                   { "const": 1 },
                   { "op": "product", "factors": [ { "param": "K0" }, { "read": "N" }, { "read": "A" } ] },
                   { "param": "K1" } ] } } }
  ],
  "nodes": [
    { "id": "investment", "type": "source",  "value": 1 },
    { "id": "D",          "type": "storage", "value": 0 },
    { "id": "M",          "type": "storage", "value": 0.01, "bounds": { "min": 0.0001 } },
    { "id": "A",          "type": "storage", "value": 0.1,  "bounds": { "min": 0.001 } },
    { "id": "N",          "type": "storage", "value": 0.3 },
    { "id": "debt_service", "type": "interaction", "module": { "k": { "param": "N5" } } },
    { "id": "heat",       "type": "sink",    "value": 0 }
  ],
  "pathways": [
    { "id": "invest",   "origin": "investment", "target": "D", "law": "constant", "k": { "param": "IV" } },
    { "id": "d_energy", "origin": "D", "target": "debt_service", "role": "energy" },
    { "id": "m_ctrl",   "origin": "M", "target": "debt_service", "role": "control" },
    { "id": "d_out",    "origin": "debt_service", "target": "heat" }
  ],
  "events": [
    { "id": "investment_stops",
      "when": { "op": "compare", "relation": "gt", "left": { "read": "D" }, "right": { "const": 30 } },
      "direction": "rising", "repeat": false,
      "actions": [ { "set_parameter": "IV", "to": { "const": 0 } } ] }
  ],
  "config": { "t_start": 0, "t_end": 320, "dt": 0.5, "outputs": ["N", "A", "M", "D", "R"] }
}
```

`N5*D*M` is a work gate: D supplies the energy and is drained, M controls and is read.
So it is the named `interaction` module, not a `process`. The rates of A, M and N
(lines 220–240) are omitted from the fragment. In full they are built the same way.

## 14. Not in v5.0

| What | Why | Path |
|---|---|---|
| `DIM` arrays (spatial or cohort models) | they need indexed components, which is a structural feature | ADR 0010's nested composites, with an index |
| events and signals added at runtime | they reorder evaluation | v5.1 |
| real exponents in `power` | `exp(e·log(b))` covers them, with log's domain check | none needed |
| `compare` `eq` | §5.2 | none intended |

## 15. Acceptance

v5 is accepted when, with tests written first:

1. Every v4 model in `examples/` migrates (§12.1) and reproduces its v4 trajectory to
   `1e-9` relative (or to the solver-change tolerance where the migration reports a
   change of method).
2. Every named module's kernel path matches its defining tree (§3.1) to `1e-12`
   relative.
3. The Odum archive's listings (MACROEC and DEVELOP first) extract to v5 through psde,
   run in GSSK, and match PC-BASIC in CI:
   - to the §6.3 tolerance for event listings
   - by ensemble for listings with random sources (§8)
   - exactly for `previous`-form verbatim extractions
4. `gssk.schema.json` v5 validates every document in 1–3, and the loader rejects each
   refusal case in this spec with the named error.

## 16. Open questions

- **Whether a control is consumed.** ADR 0016 leaves this open, and it decides a
  `process` module's output transformity as much as an interaction's. This spec
  follows current behaviour: energy supplies quality, controls do not.
- **Where the spec lives.** This draft is in GSSK-Giannantoni. It targets the GSSK
  kernel (`energese-project/GSSK`), whose published schema has external users. Its
  migration policy (ADR 0001's "deprecate neither" balance) has to be decided there.
- **Tree or flat graph.** Laws are written as nested trees here because they read like
  the expression they encode, and shared subexpressions have a home (signals). psde's
  v0.1 used flat operator nodes and edges. The two are equivalent, and psde converts
  losslessly. The question is only which one people read and diff.
