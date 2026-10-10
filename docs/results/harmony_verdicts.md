# Harmony verdicts, per seed

PLAN.md §7 (`mop-network-beta`, `mop-harmony-detector`); srs FR-HAR-002, FR-OUT-003. These are the
`harmony.<construction>: <verdict>` lines `giannantoni_sim` prints for each example seed.
`tests/mop_cli.sh` reads this table and fails if a run disagrees with a row, so the table cannot
silently go stale.

The verdicts come from the detector (`gia_harmony_verdict`, vv-plan.md §6). They are not assumed:

- **The constructions are evaluated at N = the seed's number of components.** The residual
  [23 Eq 5.6.5] needs N ≥ 3, so a seed with fewer components gets no verdict, and its report says
  "not evaluated".
- **The three constructions' verdicts are PLAN R7's.**
  - The First Equation carries its input's harmony through: `transported`.
  - The Second Equation's row and the EQS roots are roots of unity whatever the input: `imposed`.
- **The network's own Matrioska** (FR-MOP-008) is observed on the reference couple's row. Its values
  are real emergies, so for N ≥ 4 it cannot match the non-real roots (T-HAR-04, PLAN R8): `absent`.
  A seed whose pathways carry no emergy (no `quality_input` on a source) has no logarithm, and is
  not evaluated.

| Seed | Construction | Verdict |
|---|---|---|
| `examples/giannantoni/input.json` | first_equation | not evaluated |
| `examples/giannantoni/input.json` | second_equation | not evaluated |
| `examples/giannantoni/input.json` | eqs | not evaluated |
| `examples/giannantoni/input.json` | network | not evaluated |
| `examples/giannantoni/closed_loop.json` | first_equation | transported |
| `examples/giannantoni/closed_loop.json` | second_equation | imposed |
| `examples/giannantoni/closed_loop.json` | eqs | imposed |
| `examples/giannantoni/closed_loop.json` | network | not evaluated |
| `examples/giannantoni/trophic_chain.json` | first_equation | transported |
| `examples/giannantoni/trophic_chain.json` | second_equation | imposed |
| `examples/giannantoni/trophic_chain.json` | eqs | imposed |
| `examples/giannantoni/trophic_chain.json` | network | absent |
| `examples/giannantoni/coproduction.json` | first_equation | transported |
| `examples/giannantoni/coproduction.json` | second_equation | imposed |
| `examples/giannantoni/coproduction.json` | eqs | imposed |
| `examples/giannantoni/coproduction.json` | network | absent |
| `examples/giannantoni/harmonic_couples.json` | first_equation | transported |
| `examples/giannantoni/harmonic_couples.json` | second_equation | imposed |
| `examples/giannantoni/harmonic_couples.json` | eqs | imposed |
| `examples/giannantoni/harmonic_couples.json` | network | absent |

Why each "not evaluated" row:

- **`input.json`** has two components, `store_1` and `consumer_1`, so N = 2 < 3.
- **`closed_loop.json`'s network** has no residual: `source_1` has no `quality_input`, so every
  pathway carries zero emergy, and `e^α = 0` has no logarithm.

Why each network row is `absent`:

- **`trophic_chain.json`** (reference `producer → herbivore`, set in its `mop` block).
  - Its row holds the two direct pathways `producer → herbivore` and `producer → decomposer`.
  - Both carry positive real emergy, so their ratio is a positive real number.
  - The one root for N = 3 is −1, so R_H = |ratio + 1| > 1.
  - This is the only example where the network residual is computed. The run's MOP CSV shows it at
    every step.
- **`coproduction.json`** (reference `nutrients`, `seeds`: the first two by id). No pathway runs
  directly between two components; every one passes through the `growth` module. The row is
  therefore unrelated, and the residual is undefined.
- **`harmonic_couples.json`** (reference `a`, `b`). `a` has a direct pathway to `b` but not to `c`,
  so the row is not fully related.

  Its *seeded* couples are a different matter. They make the MOP CSV's `R_H` exactly 0 (to
  rounding): it is the `mop` block's First Equation row, not the network's, that is harmonic.
