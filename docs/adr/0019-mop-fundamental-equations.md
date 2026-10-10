# ADR 0019 — The MOP's two Fundamental Equations, as the engine solves them

- **Status**: accepted
- **Date**: 2026-10-09
- **Task**: `adr-mop-equations`
- **Depends on**: [ADR 0018](0018-giannantoni-test-protocol.md) — the test protocol (G4 errata, G5 definitions)
- **Records**: PLAN.md §4 R1, R4, R7 and erratum X1; srs.md FR-MOP-001…005, FR-HAR-001…003
- **Precedes**: `mop-first-equation`, `mop-second-equation`, `mop-harmony-detector`

## Context

`[22 Eq 11–12]` and `[23 Eq 4.1–4.2]` state the Maximum Ordinality Principle as two Fundamental
Equations in Relational Space. The engine has neither (PLAN §2 E8). Implementing them raises three
questions the sources do not settle by citation alone:

1. The printed solution of the First Equation, `[23 Eq 5.5.7–5.5.8]`, does not solve the equation.
2. No source writes the Second Equation, `[23 Eq 4.2]`, as an explicit differential equation. Only
   its solution is printed, `[23 Eq 6.1–6.3]`.
3. The sources say Harmony Relationships *emerge*. The engine must report whether its constructions
   produce harmony, and on what evidence.

## Decision

### 1. First Fundamental Equation: per couple, with the derived solution (R1, X1)

For each couple `(i, j)` of components, the engine solves `(d̃/dt)^k α_ij = β_ij` `[23 Eq 5.4.2]`.
The incipient derivative is `(α̇/α)^k α` `[23 Eq 5.5.2]`, `[02 Eq 14.9.5]`. Couples are independent.
The Matrioska is their collection, with a zero diagonal `[23 Eq 5.6.1]`.

The solution used is the one that follows from `[23 Eq 5.5.6]`:

```
α(t) = { (1/k) ∫₀ᵗ β^{1/k} }^k
```

For `β = (a + bt)^p`:

- **b ≠ 0:** `α = [((a+bt)^{(p+k)/k} − a^{(p+k)/k}) / (b(p+k))]^k`;
- **b = 0:** `α = (β^{1/k} t / k)^k`.

The engine evaluates both cases by numerics.md N1's unified form, which has no cancellation near
`t = 0` and is continuous through `b = 0`.

- The lower limit 0 of `[23 Eq 5.5.5]` fixes `α(0) = 0`. This is not configurable.
- Fractional k uses the principal branch, and only for real `a > 0`, `b ≥ 0` (FR-MOP-002).
- A sampled β is integrated by adaptive Gauss–Kronrod quadrature (N2), which reports its error.

**Erratum X1.** The printed `[23 Eq 5.5.7]` puts the prefactor `1/k` outside the k-th power.
`[23 Eq 5.5.8]` also drops the lower limit and misplaces `1/b`. The probe
`docs/sources/probes/eq_5_5_8_residual.py` evaluates the residual of `[23 Eq 5.5.2]` at
`a = 1, b = 0.25, p = 1, k = 2`, over `t = 0.5, 1, 2, 4`:

| Form | Residual |
|---|---|
| derived | `0, 0, 0, 0` |
| printed 5.5.7 | `1.125, 1.25, 1.5, 2.0` |
| printed 5.5.8 | `−0.5625, −0.625, −0.75, −1.0` |

Per ADR 0018 G4, T-MOP-02 asserts that both printed forms fail, at `t = 1`: residual `1.25` and
`−0.625`.

### 2. Second Fundamental Equation: the printed solution, with a reconstructed oracle (R4)

The engine implements the solution exactly as `[23 Eq 6.1–6.3]` prints it:

- `A(t) = (α₁₂(0) + λ₁₂(0)) ∘ (ᴺ⁻¹√1)^{↑N/N} + ln(c₁ + c₂t)`, with λ null `[23 §8 ii]`;
- `B(t) = [[A, −A], [−A, A]]`, the duet-binary specular form of `[22 Eq 8]` (`[23 Eq 6.2]`);
- `{r} = e^{B(t)} ∘ (roots 13 … 1N)` `[23 Eq 6.1]`.

The only time dependence is `ln(c₁ ⊕ {c₂, t})`. The duet `{c₂, t}` reduces cardinally to `c₂·t`
(R12, ADR 0020).

Its test oracle is the equation that this time dependence satisfies:
`d/dt ln(c₁ + c₂t) = c₂/(c₁+c₂t)` is the general solution of the Riccati equation `u' + u² = 0`. That
fits `[23 §6]`'s description of the result as the solution of "a typical Riccati's Equation of
Ordinal Nature". **This oracle is reconstructed from the printed solution**, and the code and docs
say so wherever it appears. The engine refuses `c₁ + c₂t ≤ 0` on `[0, t]`, where the logarithm is
undefined.

### 3. Harmony: transported, imposed, or absent, and never assumed by the solver (R7)

The engine never asserts harmony. It measures it, and it reports one verdict per construction.

**The residual.** Over ℂ:

```
R_H = max_{j=2..N−1} |α_{1,j+1}/α_{12} − e^{2πi(j−1)/(N−1)}|
```

It uses ratios to the reference couple. `[23 Eq 5.6.5]` read literally at `j = 1` gives
`α₁₂ = ω₁ α₁₂`, which collides with the reference couple itself (X9); `[23 App. A1]` "updates … the
same reference couple".

**The verdict** comes from the procedure in vv-plan.md §6. The construction is evaluated on one
harmonic input and on 8 fixed perturbations of it:

| Verdict | Condition |
|---|---|
| `imposed` | every output has `R_H` within tolerance |
| `transported` | every output's `R_H` equals its input's |
| `absent` | otherwise |

An observed Matrioska is `present` or `absent`.

**The three verdicts the sources determine.** These are not hard-coded. The detector computes them,
and VAL-07 checks that it reproduces them:

| Construction | Verdict | Why |
|---|---|---|
| First Equation | `transported` | Decoupled per couple and homogeneous of degree 1, so `α_1j = ω α_12 ⟺ β_1j = ω β_12`. The solution has exactly the harmony supplied in β |
| Second Equation | `imposed` | `[23 Eq 6.1]` writes the roots `(ᴺ⁻¹√1)_{13…1N}` into the formula. Its free inputs are `α₁₂(0)`, `λ₁₂(0)`, `c₁`, `c₂`, so couples `1j` have no input through which they could differ |
| EQS | `imposed` | By its own statement: it "originate[s] from the assumption that the Harmony Relationships…" `[23 §8 ii]` |

**Emergence is not reproducible.** The "Diffusive Generativity" from which `[22 §12]` and
`[23 §5.6, App. A1]` say harmony emerges is never written as an equation. `[23 App. A1]` itself calls
it "not a necessary consequence". The documentation says so (PLAN §6), and no verdict is `emerged`.

**No harmony in the solver** (FR-HAR-003). `src/mop.c` does not call the harmony constructor.
`bin/test_mop_emergence` links `mop.o` and `relational.o` *without* `harmony.o`, so a solver that
reached for the constructor would fail to link (T-HAR-01).

## Consequences

- The First Equation's numbers will not match `[23 Eq 5.5.8]` evaluated as printed. That is the
  point, and T-MOP-02 records it.
- The Second Equation's oracle is weaker than a printed equation: it is the equation that the printed
  solution satisfies. Wherever it appears it is labelled reconstructed, so no one mistakes it for
  Giannantoni's own.
- A network's Matrioska from real empower (FR-MOP-008, R8) is real-valued. For N ≥ 4 its ratios
  cannot equal non-real roots of unity, so the detector necessarily reports `absent` (T-HAR-04). That
  is a consequence of the sources, which assign no phase to a network couple. It is not a defect of
  the detector.

## Alternatives considered

- **Implement `[23 Eq 5.5.8]` as printed, since it is the source.** Rejected. It does not solve the
  equation it claims to solve, and G4 exists for exactly this case.
- **Report harmony as satisfied wherever the constructor built it** (today's behaviour, PLAN E5).
  Rejected. A construction checked against itself is a tautology (B11).
- **Leave the Second Equation out until a source writes it explicitly.** Rejected. Its solution is
  printed and usable, and the reconstructed oracle is stated as such rather than hidden.
