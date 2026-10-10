# ADR 0020 — The Relational Space algebra: the printed table, De Moivre exponentials, angle powers

- **Status**: accepted
- **Date**: 2026-10-09
- **Task**: `adr-relational-algebra`
- **Depends on**: [ADR 0018](0018-giannantoni-test-protocol.md); [ADR 0019](0019-mop-fundamental-equations.md)
- **Records**: PLAN.md §4 R2, R3, R12 and errata X10, X11; srs.md FR-REL-001…004, FR-MOP-006, FR-EM-006
- **Precedes**: `mop-relational-algebra`, `mop-eqs`, `emergy-ordinal-forms`

## Context

The MOP's coordinates live in Relational Space: elements `x = x_i·i ⊕ x_j·j ⊕ x_k·k` over the three
spinors of `[23 Eq 5.1.1]`, where i is the real unit. The engine today carries the harmony matrix
over ℂ, which is not Giannantoni's algebra (PLAN E5). Building the algebra raises three questions:

1. **The product table.** `[23 Eq 5.1.3–5.1.5]` (page image p. 3183) prints
   `i∘i = +1, i∘j = j, i∘k = k, j∘i = j, j∘j = −1, j∘k = k, k∘i = k, k∘j = k, k∘k = −1`. Read literally
   and bilinearly, this product is commutative and **non-associative**. Is the literal reading the
   intended one, or a misprint of an associative (quaternion-like) product?
2. **Roots of unity.** `[23 Eq A2.5–A2.6]` call the EQS roots `(ᴺ⁻¹√1)_l` "roots of unity". Under
   which product is that true?
3. **The circle product of the earlier papers** `[02 Eq 14.11.5]`, `[06b Eq 2]`. How does it relate to
   ordinary multiplication?

## Decision

### 1. The table is adopted as printed, and products are never reassociated (R2)

**The evidence is a test of the reading against the source's own use of it.** Take the EQS root in
De Moivre form, `r_l = cos(√2ψ_l) + (j + k)·sin(√2ψ_l)/√2`, which is `B_l + C_l j + C_l k` of
`[23 Eq 7.4]`. The probe `docs/sources/probes/eqs_relational_product.py` computes
`r_l ∘ {Σ₀, Φ₀, Θ₀}` under the literal table over 1000 random draws:

- It reproduces the three brackets of `[23 Eq 7.1.1, 7.2, 7.3]` to **8.9 × 10⁻¹⁶**.
- The associative alternative (`j∘k + k∘j = 0`) does **not**. For one draw it gives `−1.2116` where
  Eq 7.3 requires `−0.2686`, because it loses Eq 7.3's `+C_l(Φ₀+Θ₀)` term.

So the source's own operative equations were computed with the literal table.

The engine therefore implements:

- **`rel_mul`:** the nine-term bilinear sum of the table, in a fixed order (numerics N8).
- **`rel_mul3` and longer products:** evaluated strictly left to right in the order written, never
  reassociated. The probe's associativity witness shows the order is observable:
  `(j∘j)∘k = −k`, `j∘(j∘k) = +k`.

### 2. The exponential is De Moivre's, never a power series (R2)

```
Exp{a·i ⊕ b·j ⊕ c·k} := eᵃ [cos ρ + (b·j + c·k) sin ρ / ρ],    ρ = √(b² + c²)
```

This is the generalised De Moivre form of `[23 Eq 5.1.2]`, with the limit `eᵃ` at `ρ = 0`.

A power series in a non-associative algebra depends on how its powers are bracketed. The source
defines the exponential by the closed form, so the engine uses that. Near `ρ = 0`, `sin ρ/ρ` is
evaluated by its series to stay accurate (N8).

### 3. A root's power multiplies its angle; the table power is a different function (R3, X10)

`[23 Eq A2.5–A2.6]` define `(ᴺ⁻¹√1)_l = Exp{(α i ⊕ β j ⊕ γ k)/(N−1)}`, with periods given per spinor.
Unity there comes from **angle periodicity**.

Under repeated *table* products, the root is not a root of unity. `ordinal_root_power.py` gives, for
N = 4, l = 1, `(r∘r)∘r = (0.540721, 0, −0.665721)`, not `(1, 0, 0)`. It reaches unity only where the
sine vanishes (N = 3, or l = (N−1)/2).

The engine exposes both functions, under distinct names:

| Function | Definition | Guarantee |
|---|---|---|
| `rel_root_pow(N, l, m)` | the root at angle `m·√2ψ_l` | `rel_root_pow(N, l, N−1) = 1` |
| `rel_mul_pow(x, m)` | m-fold left-to-right `rel_mul` | (none; T-REL-04 pins the X10 value above) |

Neither may be substituted for the other; T-REL-04 pins both. This records erratum **X10**: "roots of
unity" in `[23 A2.5]` holds under angle multiplication, not under the table product.

### 4. EQS equal-angle restriction (X11)

`[23 Eq A2.6, 7.4]` give the spinor i a period of 4π. But i is the real unit `[23 Eq 5.1.3]`, and a
real exponential has no period. The EQS uses `E_{l,1}` only as a scale factor of `S_l`
`[23 Eq 7.1.1]`, and Eq 7.4's `√2ψ` collapses two angles into one, which is valid only for
`ε₂ = ε₃`. The engine therefore implements the EQS as printed and **refuses `ε₂ ≠ ε₃`** with a reason
naming X11 (FR-MOP-006).

### 5. The circle product is structural; its cardinal reduction is multiplication (R12)

- `[02 Eq 14.11.4–14.11.5]`: `l∘l ≠ l²`, and `l∘l = [l, l]`, "a du-et of real numbers".
- `[06b Eq 2]`: `(a₁; a₂) ∘ [b₁, b₂] = [(a₁b₁; a₂b₁), (a₁b₂; a₂b₂)]`.

These are one rule, written once unreduced and once reduced. The engine's circle product keeps every
pair of factors in the outer arrangement. Its **cardinal reduction** maps each pair to its product, so
`l∘l` is stored as `[l, l]` and reduces to `l²`. By the same reduction, the duet `{c₂, t}` in
`[23 Eq 6.3]` is `c₂·t` (ADR 0019 §2).

## Consequences

- Relational elements are three doubles, `(i, j, k)`. No quaternion type is introduced: the source's
  algebra is not the quaternions, and the probe shows that adopting them would break Eq 7.3.
- Code that writes `x∘y∘z` must choose `rel_mul3` (left to right) explicitly. There is no
  variadic product that could hide the order.
- Division and non-integer powers are not defined by the sources in this algebra. Relational-valued
  First Equation couples are therefore solved only for k = 1, componentwise; k > 1 is refused
  (PLAN §6, FR-MOP-007).

## Alternatives considered

- **Treat the table as a misprint of the quaternion product.** Rejected. It fails the source's own
  Eq 7.3, by the probe.
- **Define powers by repeated `rel_mul` and accept that the roots are not roots of unity.**
  Rejected as the *only* definition: `[23 A2.5–A2.6]` define the roots by angle. It is kept as a
  separately named function, because the X10 behaviour is a fact about the table worth pinning.
