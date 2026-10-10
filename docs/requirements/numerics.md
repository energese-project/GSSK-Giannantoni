# Numerical design — Giannantoni kernel

For each algorithm: the formulation the code shall use, its accuracy, and how it fails. A formulation
here is normative — a requirement in [srs.md](srs.md) that cites `Nn` is met only by this method or
one shown (by a probe and a revision here) to be at least as accurate. Probes are in
`docs/sources/probes/`.

## 1. Branch and representation policy

| Quantity | Branch / representation | Why |
|---|---|---|
| `log`, `pow` of a complex number at a single point | Principal (C99 `clog`, `cpow`) | One stated choice (NFR-NUM-004) |
| `(1+x)^q` along `t ∈ [0, T]` in N1 | `exp(q·log1p(x))` with principal `log1p`; continuous on the segment because the segment from 1 to `1+x` meets the cut only if `1+x ≤ 0` is real, which FR-MOP-002 refuses | Continuity without phase tracking |
| `β(t)^{1/k}` for sampled β (N2) | Principal at `t = 0`, continued: `|β|^{1/k} e^{iθ/k}` with θ unwrapped across samples | A straight segment not through 0 sweeps less than π |
| Fractional-order branches (FR-IDC-005) | Branch m has argument `(arg φ' + 2πm)·p/q`, `arg ∈ (−π, π]` | Matches the existing implementation and [06 Eq 3.9] |
| Relational elements | Three doubles `(i, j, k)`; i is the real unit | [23 Eq 5.1.1] |
| Integer powers | Binary exponentiation (multiplication), never `cpow` | Exact sign and branch, determinism |

## 2. Algorithms

### N1 — First Equation, affine-power boundary condition (FR-MOP-001, NFR-NUM-002)

Write `q = (p+k)/k`, `x = b t / a`. The printed closed form
`[((a+bt)^q − a^q)/(b(p+k))]^k` cancels catastrophically as `x → 0` (relative error 4.6e-8 at
`t = 1e-8`, 5.9e-5 at `t = 1e-11` for a=1, b=0.25, p=1, k=2) and needs a separate case at `b = 0`.
Use the unified form

```
S(t) = (a^{p/k} · t / k) · E(x),   E(x) = expm1(q·log1p(x)) / (q·x),   E(0) = 1
α(t) = S(t)^k
```

which is algebraically identical, keeps relative error ≤ 4e-16 across `x ∈ [1e-11, 1]`, and is
continuous through `b = 0` (probe `first_equation_numerics.py`).

- **Real a, b** use C99 `log1p`/`expm1`.
- **Complex** use `clog1p(z) = z` if `1+z == 1`, else `clog(1+z)·z/((1+z)−1)` (Kahan), and
  `cexpm1(x+iy) = (expm1(x)·cos y − 2 sin²(y/2)) + i·eˣ·sin y`.
- **`S^k`**: integer k by binary exponentiation. Rational k only for real `S > 0` (FR-MOP-002), as
  `exp(k·log S)`.
- **Overflow**: if `Re(k·log S) > log(DBL_MAX) ≈ 709.78`, return `GIA_E_RANGE` (NFR-NUM-003).
- **Accuracy**: relative error ≤ 1e-12 wherever `|E(x)|` is bounded away from 0, i.e. `1+x` not near
  0 (refused at 0).

### N2 — First Equation, sampled boundary condition (FR-MOP-004)

β is the piecewise-linear interpolant of samples `(tᵢ, βᵢ)`, strictly increasing `tᵢ`, `t₀ = 0`. A
sample or segment through `β = 0` is refused (`GIA_E_DOMAIN`). Unwrap: `θ₀ = arg β₀`,
`θᵢ₊₁ = θᵢ + arg(βᵢ₊₁/βᵢ)`; within a segment `θ(t) = θᵢ + arg(β(t)/βᵢ)`. The integrand is
`g(t) = |β(t)|^{1/k} e^{iθ(t)/k}`.

Integrate `∫₀ᵗ g` segment by segment with adaptive Gauss–Kronrod 7–15, absolute-plus-relative
tolerance `1e-10·(1 + |I|)`, at most 50 bisections per segment. Return the summed error estimate; if
any segment fails, return `GIA_E_CONVERGENCE` (NFR-NUM-005). Then `α = ((1/k)·I)^k` as in N1.

### N3 — Second-order LDE: roots, labels, integrals (FR-IDC-006)

**Roots** of `α̃² + a₁α̃ + a₀ = 0` by the stable form: `s = csqrt(a₁² − 4a₀)` (principal),
`σ = ±1` maximising `|a₁ + σs|`, `r₁ = −(a₁ + σs)/2`, `r₂ = a₀/r₁` (`r₁ = 0` only if
`a₁ = a₀ = 0`, a double root at 0).

**Double root:** declared when `|a₁² − 4a₀| ≤ 1e-12·(|a₁|² + |a₀|)` at every accepted node of the
sweep; then the solution is `c·e^{∫α̃}` (FR-IDC-006, X12).

**Sweep:** step from 0 to `t_max`. Each step `[tₙ, tₙ₊₁]` evaluates the coefficients at 5
Gauss–Legendre nodes in order, labelling the two roots at each node by the assignment minimising
`Σ|rᵢ(node) − rᵢ(prev)|`. A step is accepted when:

- its integral of each `rᵢ` agrees with two half-steps to `1e-12·(1 + |∫|)`, and
- no root moves by more than `½|r₁ − r₂|` evaluated at `tₙ`.

Otherwise halve the step, with a minimum of `1e-12·max(1, t_max)`. Reaching the minimum means the
roots collide at an isolated time: `GIA_E_CONVERGENCE`, with `why` naming the time.

Store cumulative integrals at accepted step ends. `eval(t)` integrates from the nearest stored node with
the same rule. **Initial conditions:** solve `[[1, 1], [r₁(0), r₂(0)]] c = (f₀, f₁)`; if
`|r₁(0) − r₂(0)| ≤ 1e-12·max(1, |r₁(0)|)`, refuse.

### N4 — Binary function (FR-IDC-007)

`u₁, u₂` from N3's root formula with constant `a₁ = A`, `a₀ = B`. Constants per branch by Cramer's rule
on `[[1, 1], [u₁, u₂]]`; if `|u₁ − u₂| ≤ 1e-12·max(1, |u₁|)`, refuse. Exponents `uᵢ²`; evaluation
termwise.

### N5 — Riccati by linearisation (FR-IDC-008)

y from N3 applied to `y'' + p₁ y' + p₀ y = 0` with
`p₁ = −(R' − QR)/R`, `p₀ = −P R`, `y(0) = 1`, `ỹ'(0) = R(0) f₀`. R' is supplied by the caller (a
coefficient function), not differenced. Then `f = (Σ cᵢ rᵢ Eᵢ) / (R · Σ cᵢ Eᵢ)` with
`Eᵢ = e^{∫rᵢ}`. If `|Σ cᵢEᵢ| ≤ 1e-300` or `R = 0`, return `GIA_E_RANGE` (a pole of f). The traditional
Riccati residual is reported using `f'` from the same termwise expression, not by differencing.

### N6 — Exponentials and overflow (NFR-NUM-003)

Every `e^z` is guarded by `Re z ≤ 709.78`, else `GIA_E_RANGE`. Every value returned or written is
checked with `isfinite` (real and imaginary parts) before it leaves the function; a non-finite value
is a `GIA_E_RANGE`, never an output.

### N7 — Output-projection drift (FR-IDC-014)

For each component with `Q_i(t) ≠ 0`, `Q_i' = F_i(Q(t), t)` (the flow balance). `Q_i''` is:

- **constant flow matrix:** `(A² Q)_i` exactly;
- **otherwise:** the Richardson-extrapolated central difference
  `D(h) = [F_i(Q(t+h), t+h) − F_i(Q(t−h), t−h)]/(2h)`, with `h₀ = 1e-3·max(1, |t|)` and the
  extrapolation `(4·D(h₀/2) − D(h₀))/3`. Q at `t ± h` comes from the engine's own solution, not a
  re-integration. At `t = 0` the one-sided second-order form is used.

The drift is `(Q_i'' − Q_i'²/Q_i)·Δ²/2`, with Δ the output step. If two successive extrapolations
differ by more than 1e-6 relative, return `GIA_E_CONVERGENCE`.

### N8 — Relational algebra (FR-REL)

- **`rel_mul`:** the nine-term bilinear sum of the table, in the fixed order i·i, i·j, …, k·k.
- **`rel_exp`:** `ρ = hypot(b, c)`; `sin ρ/ρ` by its series `1 − ρ²/6 + ρ⁴/120` when `ρ < 1e-4`.
- **`rel_root(N, l)`:** `ψ = 2πl/((N−1)√2)` and returns `(cos √2ψ, sin(√2ψ)/√2, sin(√2ψ)/√2)`.
- **`rel_root_pow(N, l, m)`:** returns `rel_root` at angle `m·√2ψ`, reduced mod 2π in long double
  before conversion.
- **`rel_mul_pow`:** left-to-right repeated `rel_mul`.
- **EQS (FR-MOP-006):**
  - bracket = `rel_mul(root_l, {Σ₀, Φ₀, Θ₀})` with `root_l = (B_l, C_l, C_l)` from [23 Eq 7.4–7.5];
  - `ρ₁ⱼ = A·exp(ψ₁,₁·E_{l,1}·bracket.i)` (N6-guarded);
  - `φ₁ⱼ = ψ₁,₂·E_{l,2}·bracket.j`;
  - `θ₁ⱼ = ψ₁,₃·E_{l,3}·bracket.k`.

### N9 — Incipient Taylor projection (FR-IDC-010)

`a = f₁/f₀` (refuse `f₀ = 0`); `Σ_{k=0}^{n} (aΔ)ᵏ/k!` by Horner in `aΔ`; `n ≤ 170` (factorial range),
else `GIA_E_LIMIT`.

### N10 — Harmony residual (FR-HAR-001)

Ratios `α_{1,j+1}/α₁₂`, computed with C99 complex division. If
`|α₁₂| ≤ 1e-12·max_j |α_{1j}|`, refuse (`GIA_E_DOMAIN`): the ratios are ill-conditioned. The roots
`e^{2πi(j−1)/(N−1)}` are computed from the angle in long double.

## 3. Determinism

- No algorithm depends on hash order, address order, or allocation order.
- Iteration over nodes, couples and candidates is in id order (FR-ORD-007, NFR-DET-001).
- `qsort` comparators take their context by argument, through a sort over an array of
  `(key, index)` pairs, not a file-scope pointer (NFR-REE-001).
- Cross-toolchain differences (libm, FMA contraction on arm64 clang) are absorbed by the tolerances in
  vv-plan.md §3–4. Bit-identity is required only for the same binary (NFR-DET-001).
