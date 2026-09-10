# Conformance: Odum's 1972 Energy Circuit Language

The reference is Odum, H.T. (1972), *An Energy Circuit Language for Ecological
and Social Systems: Its Physical Basis*, in Patten (ed.), *Systems Analysis and
Simulation in Ecology* Vol. II, ch. 4 — the earliest complete statement of the
language, and the one that gives each module a mathematical definition rather
than only a symbol. The PDF is in `docs/`.

Two reasons to use 1972 as the conformance target rather than *Modeling for All
Scales* (2000). It defines the modules mathematically, so conformance is
checkable rather than pictorial. And Odum presents them as a closed set — the
volume editor's introduction calls them *"a mutually exclusive and exhaustive
set sufficient for the functional representation of macroscale systems"* — which
means an omission is a real gap in expressiveness, not a stylistic choice.

This document answers two questions that are usually run together and should
not be:

1. **Does GSSK implement Odum's language?** (§1)
2. **Do Odum's module definitions and Giannantoni's IDC/MOP support each
   other?** (§2, §3) — and that question has two directions with different
   answers.

---

## 1. Module-by-module implementation status

| Odum 1972 | Definition | GSSK | Status |
|---|---|---|---|
| §II Energy Source | Forcing function; outside source of inflows | `source` (+ forcing, ADR 0006) | ✅ |
| §III Energy Pathway, barbed | Flow driven by upstream force alone | `linear` logic, `F = k·Q_origin` | ✅ |
| §III Energy Pathway, barb-less | Flow driven by the difference of forces at both ends | `reversible` logic (ADR 0007) | ✅ |
| §V Heat Sink | Pathway of used energy; absorbs, not depleted | `sink` | ✅ |
| §VI Passive Energy Storage | Tank; `dQ/dt = ΣJ_in − ΣJ_out`, force from `Q` | `storage` | ✅ |
| §VII Potential-Generating Storage | Storage whose filling does work against its own backforce, with heat dispersal; max power when backforce is half the input force | `consumer` composite (`body` + `heat`) | ⚠️ partial — the storage and its dissipation are modelled; the **maximum-power backforce adjustment is not** |
| §IX Constant Gain Amplifier | `J = k·N`; control signal, separate energy source supplies the power | `gain` | ✅ |
| §X Work Gate | `J = k·N₁·N₂`, multiplicative junction | `interaction` (n-ary and subtracting, ADR 0008) | ✅ |
| §XI Switch Module | Discontinuous / on-off action | `switch` node, `threshold` logic | ✅ |
| §XII Self-Maintaining Module | Multiplicative feedback onto the upstream flow; logistic growth; *"autocatalytic by chemists and logistic by population biologists"* | `producer` composite (`body` + `gate` + `heat`) | ✅ — see naming note below |
| §XIII Cycling Receptor Module | Multiplicative feedback constrained by a recycling material constant over short periods; hyperbolic saturation; photocell and photosynthesis | `loop_limited` node, `limit` logic | ✅ |
| §XIV Production and Regeneration (P-R) | *"formed by combining a cycling receptor module, a self-maintaining module which it feeds, and a work feedback loop"*; **the green plant** | — | ❌ **absent** |
| §XV Economic Transactor | `J_energy = k·J_currency`, counter-flowing money and energy | `exchange` (transaction diamond, ADR 0001) | ✅ |

**Score: 11 of Odum's 12 modules implemented; 1 absent (§XIV P-R); 1 partial
(§VII).**

Beyond the 1972 set, GSSK also ships `constant` (no Odum symbol), and `ratio`
logic for the divisor action of *Modeling for All Scales* Fig. 2.6d (ADR 0002).

### The `producer` naming note

What GSSK calls `producer` is Odum's §XII **Self-Maintaining Module**. Odum's own
producer — the green plant — is §XIV **P-R**, which is not implemented.

The name is deliberately kept (ADR 0010, decision 7). `producer` is what
modellers call a self-maintaining autotrophic unit, and the 1972 module names
are a physical taxonomy rather than a modelling vocabulary. The consequence to
remember: **if P-R is ever added it must take a different type name**, because
the obvious one is taken.

### Why §XIV is not merely a missing feature

P-R cannot currently be *written*, not just not-shipped. Odum defines it as a
combination of two other modules, and a GSSK archetype template may only name
fundamentals — `ArchetypeDefn.nodes.items.type` is
`{"$ref": "#/$defs/PrimitiveNodeType"}`. Since `loop_limited` is a fundamental
but the self-maintaining half is the `producer` *composite*, the definition has
no expression. This is the motivating case for ADR 0010.

---

## 2. Direction one — does Odum's mathematics support Giannantoni's IDC?

IDC's persistence-of-form theorem is stated for the exponential form
`f(t) = e^φ(t)`, where `(d̃/d̃t)ⁿ e^φ = (φ°)ⁿ e^φ`. So the question per module is:
**is that module's law expressible as `e^φ` with `φ` constructively available?**

| Odum module | Law | Exponential form | IDC treatment in GSSK |
|---|---|---|---|
| §II Source | forced / fixed | `φ` affine when constant | **exact** |
| §III barbed pathway | `F = k·Q` | `φ` affine | **exact** |
| §III barb-less | `F = k·(Q_a − Q_b)` | linear, exactly integrable | **exact** |
| §V Heat Sink | accumulation | linear | **exact** |
| §VI Passive Storage | linear balance | `exp(A·t)` | **exact** |
| §VII Potential-Generating Storage | linear + dissipation | linear | **exact** for the part implemented |
| §IX Gain | `J = k·N` | linear in control | **exact** |
| §X Work Gate | `J = k·N₁·N₂` | **bilinear** | **exact only as the duet** — closed Riccati form for an isolated 2-node pair; Padé linearisation otherwise. The general case is the **n-et**, open for n ≥ 3 |
| §XII Self-Maintaining | `dQ/dt = aQ − bQ²` | logistic / Riccati | **exact when isolated** — this is the duet |
| §XIII Cycling Receptor | `J = k·N·I/(C+N)` | **rational** | **approximate, bounded** — effective conductance `g = k·C/(C+Q)`, with per-edge error reported |
| §XI Switch | discontinuous | **no smooth φ** | **piecewise only** — `φ` has a kink at the event, so persistence-of-form holds *between* events; GSSK locates each crossing by Illinois iteration and restarts |
| §XV Economic Transactor | `J_e = k·J_$` | linear if price constant | **exact** for constant price; time-varying under endogenous price |

### What the table supports

Three distinct outcomes, and they should not be blurred:

- **The linear core is fully supported and genuinely benefits.** Source, pathway,
  storage, sink, gain and transactor are linear or linear-time-varying, which is
  exactly the class where IDC replaces stepping with one matrix exponential.
- **The multiplicative junction is the boundary.** §X and §XII are bilinear. IDC
  solves the two-body case exactly and leaves n ≥ 3 open — and that is open **in
  Giannantoni's own papers**, not merely unimplemented here.
- **Two modules resist the form.** §XIII is rational: formally any positive `f` is
  `e^{ln f}`, but `φ` is only *constructively* available once you have already
  solved the ODE, so universality here is vacuous. §XI is discontinuous, so no
  `φ` exists across the event at all.

### One deep convergence, from 1972

The volume editor's introduction to the 1972 chapter says the modules are joined
by *"additive (corresponding to conservative) and multiplicative (not
corresponding) interactions."*

That is the same boundary Giannantoni draws thirty-four years later. His argument
for IDC is that emergy algebra is **non-conservative** — the output of a
generative process is irreducible to its inputs — and that TDC cannot express it
because classical derivatives impose conservative relationships.

**Both frameworks locate the non-conservativeness in the multiplicative
junction.** Odum's work gate is Giannantoni's generative step. That is not a
coincidence to be smoothed over; it is the strongest evidence that the two are
talking about the same thing, and it was written down in 1972.

It is also why §X and §XII are exactly where the mathematics gets hard in both.

---

## 3. Direction two — does Giannantoni's framework support Odum's language?

**As published: no.** This direction is much weaker than the first, and it is
worth being blunt about.

Odum's language is **heterogeneous by construction**. Each of the dozen modules
has *its own mathematical definition*, and the editor's claim that they are
"mutually exclusive and exhaustive" is a claim that the typing is load-bearing —
a work gate is not a gain is not a switch, and none reduces to another.

Giannantoni's MOP is **homogeneous**. Its state is Relational Space coordinates
`{r}_s = e^{σi ⊕ φj ⊕ ϑk}`, and its structure is the N×N matrix of ordinal
relationships `α_ij` reduced by the Harmony Relationships to a single reference
couple. Components enter as *relata*, not as typed modules.

There is therefore **no MOP representation of**:

- a **switch** — MOP has no discontinuity; ordinal relations are smooth
- a **constant gain amplifier** — the distinction between a control signal and
  the power source that amplifies it has no ordinal counterpart
- an **exchange / transactor** — counter-flowing carriers coupled by price is a
  two-carrier construct; MOP has one relational space
- the **conservative/non-conservative distinction between module kinds** — MOP
  asserts non-conservativeness globally rather than locating it in particular
  junctions, which is precisely what Odum's additive-vs-multiplicative typing
  does

Mapping Odum onto MOP therefore requires one of:

1. **Flatten the typing** — express every module as ordinal relations and lose
   the distinctions Odum calls exhaustive. This is what the current
   `src/engine.c` effectively does: it reads five node kinds and reduces them to
   a `φ` polynomial each.
2. **Extend MOP with typed relations** — not in the 2006 or 2023 papers, and a
   research undertaking rather than an implementation task.

### Consequence for this repository

The asymmetry explains the current shape of the code, and it is not an accident:

- **`src/gssk.c` is a faithful Odum simulator** — 11 of 12 modules, typed,
  with the IDC solver applied where each module's law permits.
- **`src/engine.c` is a Giannantoni MOP engine** — it implements the calculus and
  the ordinality machinery, and it does *not* implement Odum's language. Its
  five node kinds (`source`, `interaction`, `storage`, `consumer`, `regulator`)
  are not the 1972 modules, and it has no `gain`, `switch`, `exchange`,
  `loop_limited`, `sink` or `constant`.

### What `src/engine.c` implements, and what it still lacks

An earlier revision of this section recorded a much larger gap: no network at
all, node kinds that were only labels selecting a `phi` template, no edge logic,
and a trajectory that was byte-identical with every edge deleted. That has been
closed. The record of it is kept in
`docs/mop_network_coupling_requirements.md`, because the shape of the mistake is
worth remembering: the engine demonstrated the drift theorem correctly and was
presented as a simulator, which it was not.

**Now implemented.** The exponential form for a network is the matrix one --
`{r}_s = e^{α(t)}`, so `Q(t) = exp(A·t)·Q(0)` where `A` is assembled from the
pathway laws below. The matrix is augmented to `(n+1)×(n+1)` so that pathways
delivering a rate rather than a conductance stay exact rather than being
approximated away.

| Odum law | GSSK logic | In `engine.c` |
|---|---|---|
| §III barbed, `F = k·Q_origin` | `linear` | ✅ |
| §III barb-less, `F = k·(Q_a − Q_b)` | `reversible` | ✅ four matrix entries |
| §X work gate, `F = k·Q_a·Q_ctl` | `interaction` | ✅ linearised; duet exact, n-et open |
| §IX amplifier, `F = k·Q_control` | `gain` | ✅ entry in the control's column |
| §XIII receptor, Michaelis-Menten | `limit` | ✅ `g = k·C/(C+Q)` |
| §XI switch, `F = k` above threshold | `threshold` | ✅ **piecewise, Illinois event location** |
| §XV transactor, `J_e = P·J_$` | `exchange` | ✅ counter-flow at `F/P` |
| Fig 2.6d divisor (ADR 0002) | `ratio` | ✅ `g = k / max(Q_c, ε)` |
| Fig 2.6e subtracting (ADR 0008) | `subtract` | ✅ clamped, clamp is an event |
| fixed rate | `constant` | ✅ affine, via the augmented column |

All nine GSSK primitive node types are accepted: `source` and `constant` are
held rather than integrated (§II), a `sink` is never depleted (§V), and
`storage`, `interaction`, `gain`, `loop_limited`, `switch` and `exchange`
integrate.

Two structural results follow, and both are tested against closed forms:
conservation holds exactly for a closed system, and `ψ` is **exactly zero** when
the flow matrix is constant — the matrix analogue of "affine φ is drift-free" --
and non-zero precisely at the multiplicative junction.

**Emergy and transformity are implemented**, as a second accounting pass over
the same topology under a different algebra. This was the most consequential
gap, because emergy's non-conservativeness is Giannantoni's *stated motivation*
for IDC — the engine had his calculus without the quantity that motivated it.

Odum's rules, using the neutral vocabulary of `docs/emergy_synthesis.md`:

- a source injects quality at its declared `quality_input`;
- emergy travels a pathway as `F x Tr_origin`;
- `output_mode: "partition"` (the default) is a **split** — each branch takes
  emergy in proportion to its share of the energy, so emergy out equals emergy
  in and the accounting is conservative;
- `output_mode: "replicate"` is a **co-production** — each product takes the
  **whole** emergy, because each required all of it. Emergy out exceeds emergy
  in, and the difference is the irreducible excess;
- an edge closing a cycle carries quantity but re-injects no emergy, so feedback
  is not counted twice (Odum's fourth rule).

`gia_emergy_excess()` reports the emergy created across the network. It is
exactly zero when every bifurcation is a partition and positive as soon as one
is a replication — the non-conservativeness reported as a number rather than
argued about. Odum's own wind-and-rain example gives 10 units in and 10 units
created:

```
  component      Em        Tr
  sun            10         1
  atm            10    1.2254
  wind           10   2.92423
  rain           10   2.92423
  emergy created 10
```

Quantity remains conserved over the same run. The two accountings share the
topology and disagree about the algebra, which is the point.

**Still absent.**

1. **Leg discovery, money-stock gating, and a price resolved from a node**
   (ADR 0001). Carriers now exist: a component declares what it holds,
   conservation is checked per carrier, a pathway may not cross carriers, and
   an `exchange` must couple two different ones with both currency legs on the
   counter-carrier. What remains is the rest of the kernel's diamond.

   An exchange is **not** required to couple two different carriers. Barter is
   a real process — grain for sheep, or one commodity traded between two
   markets at a ratio — and Odum's §XV transactor is written for money only
   because that is the case he was modelling. The structure it describes, two
   counter-flowing quantities coupled by a ratio, does not depend on either
   side being money. An earlier revision rejected a same-carrier counter-flow;
   worse, it did so only when the model happened to contain a second carrier
   somewhere else, so one edge's validity depended on unrelated parts of the
   graph. Both are gone. The field names follow suit: `counter_origin`,
   `counter_target` and `exchange_ratio` are the neutral spellings, with
   `currency_origin`, `currency_target` and `price` accepted as aliases.

   What is still rejected: a non-exchange pathway crossing carriers, counter-
   flow legs holding different carriers from each other, and both legs landing
   on the same component, which pays nothing to nobody.

   Leg discovery is deliberately **not** inferred. Given a goods flow a -> b and
   a money carrier, which money component is the payer and which the receiver
   is not recoverable from carrier identity alone — it needs the ownership the
   diamond's shape encodes. Guessing would be inventing semantics Odum's figure
   does not show, so the legs stay named and are validated by carrier instead.

   The change is worth stating plainly, because the old behaviour was not
   merely incomplete. `gia_conservation_residual` summed every integrating
   component into one total, so a model losing 2 units of grain while gaining
   2 units of money reported a residual of **exactly zero**. A surplus in one
   carrier cancelled a deficit in another. The residual is now the worst
   carrier rather than the sum, and a regression test pins that the old summed
   figure would have hidden it.
3. **Forcing** (ADR 0006) — both attachment points are now implemented, and
   the difference between them is the sharpest statement of what drift is.

   A driven **value** is additive, `dQ/dt = A Q + k S(t)`. If the waveform
   generates itself it is absorbed into `A`, the augmented system stays
   time-invariant, and `psi` is exactly zero.

   A driven **rate** is multiplicative, `dQ/dt = k(t) Q`. That is bilinear in
   driver and state, and no augmentation linearises it. So it is the one case
   in the suite where `psi` is non-zero with no interaction, limit, ratio,
   threshold or subtract edge anywhere in the model — drift arising from time
   dependence alone, which nothing else here can show.

   It is also where the closed form ends. The solution is composed over
   subintervals with a midpoint step, and `gia_integration_error()` reports
   what that cost, Richardson-estimated and scaled by 4/3 because the raw
   N-versus-2N difference is only three quarters of a second-order method's
   error. The test checks the estimator against the error it can actually
   measure against the closed form — same size, not merely small — because a
   bound that understates is worse than none when `psi` is read against it. In
   the test model `psi` is some four orders of magnitude above the integration
   error, so it is measuring the calculi rather than the solver.

   Node attachment: a source or constant may declare a `forcing`
   waveform, and `sine`, `ramp` and `exponential` are carried exactly.

   The reason those three and not the other five is worth stating, because it
   is the same reason the engine is exact at all. Each of them **generates
   itself** — `d/dt [s;c] = [[0,w],[-w,0]][s;c]`, `d/dt r = 1`,
   `d/dt e = lambda e` — so the driver can be carried as extra state and the
   augmented system stays linear and *time-invariant*. `Q(t) = exp(A t) Q(0)`
   therefore remains a closed form under forcing, checked against
   hand-integrated solutions (`A(1 - cos wt)/w` and so on). A square wave,
   sawtooth or jitter is not its own generator, so it is refused rather than
   quietly turned into a step-and-hope.

   **A driven model still has psi exactly zero, and that is a result rather
   than an oversight.** Drift is about the coefficient depending on the
   *state*, not on time: a driver absorbable into `A` leaves the two calculi in
   agreement. The case that does drift without any multiplicative junction is
   Odum's *other* attachment point — forcing an edge RATE, where the flow is
   `k(t)·Q`, bilinear in driver and state and therefore not absorbable. That is
   what remains, and it is the interesting half for the drift claim.
4. **Composites and archetypes.** `producer`, `consumer`, `misc_box` and
   `system_frame` are refused with a message pointing at ADR 0010, rather than
   silently given storage semantics.
5. **The n-et for n ≥ 3.** Open in Giannantoni's own papers; multiplicative
   chains take the Padé linearisation with `ψ` reported.

### The Level 1 claim, measured

`docs/giannantoni_assessment.md` §5.3 makes a claim it flags as unquantified:
that any system modelled through the MOP lens has an explicit solution. The
projection of ADR 0011 turns that into a number. Over the 24 GSSK models in
`examples/`:

    mean coverage 92.3%     none below 71%

What stops the remaining 7.7% is not scattered. It is four named things, each
already recorded:

| blocker | models | recorded in |
|---|---|---|
| a processing node's law lives in its own `params` | 3 | below |
| endogenous price resolved from a node | 4 | ADR 0001 |
| composites and user archetypes | 2 | ADR 0010 |
| n-ary interaction, more than one control node | 1 | ADR 0008 |

The first is an architectural difference worth stating on its own, because it
was not visible until the projection measured it. **GSSK configures a
processing node through the node's own `params` block; this engine puts laws on
pathways.** A GSSK `interaction`, `gain`, `loop_limited`, `exchange` or
`switch` node carrying parameters has nowhere to put them here, and a node
emitted without them is a storage wearing the name. The projection drops it and
says so, rather than producing a model that runs and is not the one anybody
wrote.

That rule matters more than the percentage. An earlier revision of the
projection dropped node `params` silently and scored these same models at
**97.2%** — a figure that was higher, and wrong.

**Odum's mathematics can feed Giannantoni's calculus, module by module, with the
limits tabulated in §2. Giannantoni's framework cannot yet carry Odum's
vocabulary.** Any convergence of the two engines has to solve direction two
first, and direction two is open research, not integration work.

---

## 4. Summary

| Question | Answer |
|---|---|
| Does GSSK implement Odum's 1972 language? | **11 of 12 modules.** §XIV (P-R) absent and currently inexpressible; §VII partial (no max-power backforce). |
| Does Odum's mathematics support IDC? | **Yes for the linear core** (exact, and this is where the O(steps) → O(1) gain is real). **Duet-exact** for the multiplicative junction, open for n ≥ 3. **Bounded-approximate** for the cycling receptor. **Piecewise only** across a switch. |
| Does Giannantoni's framework support Odum's language? | **No, as published.** MOP has relational coordinates, not typed modules; there is no MOP form for a switch, a gain, or an exchange. |
| Are the two frameworks aligned in principle? | **Yes, and demonstrably since 1972** — both locate non-conservativeness in the multiplicative junction. |

### Open items this document identifies

- Implement §XIV P-R (blocked on ADR 0010 — archetype nesting)
- Model §VII's maximum-power backforce adjustment, or document its omission as
  deliberate
- Whether the two engines converge is settled by
  [ADR 0011](adr/0011-two-engines-declared-lossy-projection.md): they stay
  separate, joined by a one-way declared-lossy projection that reports a MOP
  coverage number per model. §3 above is that ADR's source analysis
