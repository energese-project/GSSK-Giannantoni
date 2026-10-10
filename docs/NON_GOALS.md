# Non-goals

Things this project deliberately does not build. Each entry says why, so that a request for one can be
answered by pointing here rather than by re-arguing it. TODO.md schedules a quarterly review of this
list. If a non-goal keeps being requested, that review is where it is reconsidered, with an ADR, not by
quietly building it.

## The kernel (`gssk.h`)

- **No stiff-solver zoo** (BDF, Rosenbrock, SDIRK). A problem that is stiff and linear with constant
  coefficients is already handled by the matrix exponential (`"method": "expm"`). A problem that is
  stiff *and* not of that form is a modelling problem, not a solver gap.

  This entry used to say "IDC-eligible", which equated the incipient calculus with the matrix
  exponential. They are different things (PLAN.md §1 B5, ADR 0022).
- **No DAE or index-reduction support.** Models are ODEs over storages. An algebraic constraint
  belongs in the model's structure, not in the solver.
- **No symbolic Jacobian engine.**
- **No automatic-differentiation framework.** Adjoint sensitivity is hand-coded per flow primitive.
- **No general scientific-computing parity** with SciPy or DifferentialEquations.jl. The kernel is an
  Odum energy-systems simulator, not a general ODE library.
- **No GUI builder.** The schema is the UI contract (`gssk.schema.json`). Visualisation is a
  downstream concern.

## The Giannantoni engine (`engine.h`, `idc.h`, `mop.h`)

- **Nothing the sources do not define.** PLAN.md §6 lists what is excluded and why. Examples: the
  "emergence of Harmony" (never written as an equation), the solar-system and N-body closed forms
  (their parameters are not given), the Abel n-et, and the relational-valued First Equation for
  k > 1. The engine refuses each by name (FR-OUT-002). A request to "just approximate it" is a
  request to invent a result and attribute it to the sources; ADR 0018 forbids that.
- **No second model format.** The projection (`gia_project`) is a read-only derivation from a GSSK
  model, with no authoring surface of its own (ADR 0011). Models are written against the GSSK schema,
  or as Giannantoni seeds in `examples/giannantoni/`.
- **No merged engine.** The two engines share a build and a JSON parser, nothing else (ADR 0011,
  NFR-SEP-001). `bin/gia_bridge` compares them. It does not join them.
