# ADR 0017 — A control input is a flow of energy, and its emergy reaches the product

- **Status**: proposed
- **Date**: 2026-09-17
- **Task**: `control-input-consumed`
- **Supersedes**: the definition of `control` in [ADR 0013](0013-module-roles-not-position.md) decision 2
- **Depends on**: [ADR 0013](0013-module-roles-not-position.md) — roles;
  [ADR 0014](0014-ordinality-over-quantity-legs.md) — ordinality over legs that carry quantity;
  [ADR 0016](0016-interaction-actions.md) — interaction actions

## Context

ADR 0013 gave every module two input roles. `energy` is consumed. `control` is
"read. The module's law depends on it; it is not depleted." It attributed that
to Odum: "Odum draws the same distinction in every module of §IX–§XIII: one
input supplies the energy and is used up; the others signal, and are read but
not depleted."

ADR 0014 then built on it: "A `control` leg is read and never consumed," so it
closes no pathway and crosses no boundary. The emergy pass followed too: a
control leg's flow is zero.

Review of ADR 0016 brought Odum's own figure and definition to the question, and
they do not say what ADR 0013 says.

### What the symbol says

The definition of the interaction symbol in *Modeling for All Scales* opens:
"Two or more flows that are different and both required for a process are
connected to an interaction symbol." Fig. 2.6 (p. 29) labels the inputs of panel
(a) "Use of Ingredient A" and "Use of Ingredient B", with a "Used Energy Flow"
to the heat sink. Panel (b) shows a "Low Transformity Input" and a "High
Transformity, Control" producing "**Moderate** Transformity Products".

### What Odum 1972 says

The modules ADR 0013 cites, read in the 1972 paper:

- **§X, work gate.** "The work of one **flow of energy** controlling and
  facilitating conductivity of a second." "The limiting and controlling **lesser
  energy flow** is shown passing in from above. The heat sink represents the
  usual spontaneous heat increases associated with any process **including those
  from the coupling of the two forces**. The smaller flow is the signal in an
  amplifier effect … The smaller control flow may come from another part of the
  network or from outside energy sources."
- **§X.C, multiplicative junctions.** Both inflows are drawn in proportion to
  the reaction, `dN₁/dt = J₁ − sJ₃`, "where s is the stoichiometric ratio … of
  necessary use". And: "There are, however, many double-flux processes where
  only one of the flows contributes material to the outflow, **the other inflow
  contributing only energy**."
- **§IX, amplifier.** The combination follows the second energy principle "by
  dispersing more potential energy into heat in the amplifier flow than is
  stored in the flow amplified". "If the gain is high, the down-circuit process
  may be guided by the up-circuit process **without drawing much power** from
  the input" — much, not none.
- **§XI, switch.** "Work done by one or more energy circuits controls another."
  "The switch requires maintenance energy by at least one of the inflowing
  processes and … has a heat drain."
- **§VIII, feedback.** Odum does have an influence that carries no flow, and he
  distinguishes it by exactly that: "The feedback may be a pathway of flow of
  force and energy … Also there are feedbacks of effect **without a special
  pathway**."

Odum does call the control the signal. But in the same sentence it is "the
smaller flow". ADR 0013 took "signal" to mean read and not depleted, and nothing
in §IX–§XI says that. Every one of those sections describes the control as a
lesser flow of energy that does work and ends up as heat.

### What emergy algebra says

The rules as given by Brown and Herendeen (1996), quoted in Bastianoni et al.
(2011):

> Rule 1. All source emergy of a process is assigned to process output.
>
> Rule 4. Emergy cannot be counted twice within a system: (a) emergy in feedback
> cannot be double counted; (b) co-products, when reunited, cannot add up to a
> sum greater than the source emergy from which they were derived.

A control is a source to the process, so its emergy is assigned to the product.
Panel (b)'s "moderate" is Rule 1: a small flow of high transformity added to a
large flow of low transformity gives a product between the two.

### Measured on `main`

A work gate fed by `sun` (transformity 1) on its energy leg, with a source `H`
on its control leg:

| `H`'s transformity | empower at product | product transformity |
|---|---|---|
| 1 | 6.0000 | 1.0000 |
| 1000 | 6.0000 | **1.0000** |

The control's transformity makes no difference to the product. Because a control
leg's flow is zero, Rule 1 never sees it. Whatever the control is, the product
gets exactly the energy input's transformity, where panel (b) says it should be
moderate.

## Decision

1. **A control pathway is a flow. It draws `J_c = s · F` from its source.**
   `F` is the module's flow, and `s` is the pathway's `use_ratio` — §X.C's
   "stoichiometric ratio of necessary use". The form follows Eqs. (63)–(65):
   every flow at the junction is proportional to the same product of forces, with
   its own coefficient. This applies to every control leg: a work gate's under
   each of ADR 0016's actions, an amplifier's, and a switch's sensor.

   ```json
   { "source": "H", "target": "gate", "role": "control", "use_ratio": 0.01 }
   ```

2. **The control's quantity is dissipated, not added to the product.** This is
   §X.C's second case: the control is "the other inflow contributing only
   energy", and §X, §IX and §XI all send it to the heat sink. The energy input's
   quantity still becomes the product. The control's leaves on the module's
   `used` leg (decision 5).

3. **The control's emergy reaches the product.** This is Rule 1. The product's
   emergy includes `J_c · Tr_c`, combined with the energy input's under the
   engine's existing Rule 4 handling, so a control fed back from downstream is
   not counted twice. The `used` leg carries no emergy onward.

   Worked through for the measurement above, with `s = 0.01` and no feedback:

   ```
   Em_product = F · Tr_energy + s · F · Tr_control
   Tr_product = Tr_energy + s · Tr_control = 1 + 0.01 · 1000 = 11
   ```

   11 lies between 1 and 1000, which is panel (b)'s "moderate". In general the
   product stays below the control while `s < 1 − Tr_energy / Tr_control`, and
   Odum's control is "the lesser energy flow". At `s = 0` the product is exactly
   1: the read-only limit is precisely the case the figure contradicts.

4. **`use_ratio` is required on every control leg, and `0` is allowed.** A zero
   means the control's energy is neglected, and the report names every such leg.
   The value is required rather than defaulted because the measured defect *is*
   a silent zero, and ADR 0013's own rule is that nothing about a module's inputs
   is silently defaulted. Odum's name for an influence carrying no flow is a
   "feedback of effect without a special pathway". A control leg is a drawn
   pathway, so if it carries nothing the model should say so.

5. **The used-energy flow belongs to the module.** Every panel of Fig. 2.6 draws
   it. A module with any control whose `use_ratio > 0` names a `used` output leg,
   ending at a sink that holds that control's carrier. Its flow is `Σ s_i · F`
   over the controls of that carrier. A missing `used` leg is a load error naming
   the module and the carrier.

   A price read as a control is the case this catches. Money is not dissipated,
   so a price control has `use_ratio: 0`. Writing a nonzero one forces the
   modeller to draw a money sink, which is the prompt to reconsider.

## Consequences for merged decisions

- **ADR 0013.** The roles keep their names, their number, and the arity each
  module requires. What changes is the definition of `control`, from "read, not
  depleted" to "drawn at `use_ratio · F` and dissipated". Order-invariance is
  unaffected. ADR 0013's attribution of the old definition to §IX–§XIII is
  superseded by the readings above.
- **ADR 0014.** At `use_ratio: 0`, nothing changes. Above zero, a control leg
  carries quantity, but only through the module to its `used` leg and so to a
  sink, never into the product's outputs. The cycle scan therefore passes a
  control through to `used` and nowhere else. Sinks are dead ends, so a control
  still closes no pathway through the product, and every reproduction in ADR 0014
  keeps its verdict. Decision 3 there is amended: a control with `use_ratio > 0`
  read from a source or a constant draws quantity across the boundary, so the
  system is open.
- **ADR 0015.** Unaffected. `used` legs end at sinks, which the emergent quality
  never closes.
- **ADR 0016.** `divide` and `subtract` controls draw `s · F` like `multiply`'s.
  When `subtract`'s clamp holds `F` at zero, the control draws nothing.
- **The emergy pass.** A control leg's flow is `s · F` rather than 0, the
  product's incoming emergy includes it, and the `used` leg carries none.
- **Models.** Every control leg in the suite gains a `use_ratio`, and each model
  that sets one above zero gains a `used` leg and a sink. This is the breaking
  change ADR 0012 accepted for a repository with no external users. The
  migration writes them.
- **Pathway laws.** The pathway `interaction` law's `control_node` has no drain
  and no `used` leg. It is removed under ADR 0012, so no pathway spelling of
  `use_ratio` is added.

## What this does not decide

- **Dissipation on the energy input.** Odum's heat sink also takes heat from the
  main flow. The module still delivers all of `F` to its outputs. An efficiency
  on the energy leg would change every trajectory, and is its own change.
- **More than one input contributing material to the product.** This covers
  Fig. 2.6(a)'s two ingredients and §X.C's stoichiometric limiting factors,
  Eqs. (66)–(69). ADR 0013's single `energy` input stays.
- **Deriving `use_ratio`.** It is a stated parameter, like `k`, not computed
  from transformities or anything else.
- **A separate construct for influence without flow.** A control with
  `use_ratio: 0` is that influence, stated on the leg. Whether it also deserves
  an ordinality of its own is the question ADR 0014 already left open.
