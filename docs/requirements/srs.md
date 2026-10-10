# Software requirements specification — Giannantoni kernel

Conventions are in [README.md](README.md). Equations are cited `[key Eq n]` with keys from
[docs/sources/README.md](../sources/README.md); `PLAN Rn` and `Xn` refer to PLAN.md §4 and §5. Numerical
methods (`N1`…) are in [numerics.md](numerics.md); interfaces (`IF-…`) in [icd.md](icd.md); test IDs in
[vv-plan.md](vv-plan.md).

Throughout, "the kernel" means the Giannantoni library units (`src/engine.c`, `src/idc.c`,
`src/mop.c`, `src/relational.c`, `src/harmony.c`) and `bin/giannantoni_sim`. "Shall refuse" means:
return `GIA_E_UNSUPPORTED` or `GIA_E_DOMAIN` (IF-API-001) with a reason naming the source, produce no
partial output, and never fall back silently to another method.

---

## 1. Definitions

| Term | Meaning here |
|---|---|
| Incipient derivative of cardinality n | Of a single function given pointwise: `(d̃/dt)ⁿ f = (f'/f)ⁿ · f`; for `f = e^φ`, `(φ')ⁿ e^φ` [02 Eq 14.8.2, 14.9.5], [06 Eq 3.2], [23 Eq 3.2.5] |
| Incipient derivative of a superposition | Of a finite sum of exponential terms: termwise, `(d̃/dt)ⁿ Σ cᵢ e^{φᵢ} = Σ cᵢ (φᵢ')ⁿ e^{φᵢ}`. The pointwise form is not additive, and only the termwise reading makes the superpositions the sources print solve their equations: [06 Eq 3.6] has residual −1.67 pointwise, 0 termwise (probe `incipient_superposition.py`). Every superposition here — FR-IDC-006, FR-IDC-007, FR-IDC-011 — is termwise |
| Traditional derivative | The ordinary derivative; for `e^φ`, Faà di Bruno: `e^φ · Bₙ(φ', …, φ⁽ⁿ⁾)` |
| Solution drift | Difference between the incipient and traditional solutions of a *known* equation [06 §4] |
| Output-projection drift | `[09 Eq 13]`: per-term difference between the traditional and incipient Taylor series of an *output* trajectory |
| Component | A non-module, non-boundary node; boundary = `source`, `sink`, `constant` (PLAN R5) |
| Couple | An ordered pair `(i, j)`, `i ≠ j`, of components |
| Matrioska | The N×N array of couple coordinates `α_ij(t)`, zero diagonal (internal representation [23 Eq 5.6.1]) |
| Relational element | `x = x_i·i ⊕ x_j·j ⊕ x_k·k` over the spinors of [23 Eq 5.1.1] |

---

## 2. Functional requirements

### 2.1 Incipient calculus (single variable)

### FR-IDC-001 — Incipient derivative of an exponential
For `φ(t)` a real polynomial of degree ≤ 8 and integer `0 ≤ n ≤ 8`, the kernel shall compute
`(d̃/dt)ⁿ e^φ = (φ'(t))ⁿ · e^{φ(t)}`.
- **Source:** [02 Eq 14.8.2, 14.9.4], [09 Eq 12, 25], [22 Eq 5.2], [23 Eq 3.2.5]
- **Verification:** T-IDC-01 · **Priority:** Must · **Status:** implemented · **Task:** —

### FR-IDC-002 — Incipient derivative of a general function
Given `f(t) ≠ 0` and `f'(t)` in ℂ and integer `n ≥ 0`, the kernel shall compute `(f'/f)ⁿ · f`. For
`f(t) = 0` it shall refuse (`GIA_E_DOMAIN`).
- **Source:** [02 Eq 14.9.5], [23 Eq 5.5.2]
- **Verification:** T-IDC-03 · **Priority:** Must · **Status:** implemented · **Task:** idc-general-f

### FR-IDC-003 — Traditional derivative of an exponential
For the same φ and n as FR-IDC-001, the kernel shall compute the traditional derivative as
`e^φ · Bₙ(φ', …, φ⁽ⁿ⁾)` using the complete Bell polynomial.
- **Source:** [02 p. 175], [06 Table 1], [10 Table 1], [22 Eq 5.1]
- **Verification:** T-IDC-01 · **Priority:** Must · **Status:** implemented · **Task:** —

### FR-IDC-004 — Drift identity
The kernel shall compute `ψₙ = Bₙ(φ', …, φ⁽ⁿ⁾) − (φ')ⁿ` and shall report a polynomial φ as drift-free
exactly when its degree is ≤ 1, decided from the coefficients without evaluation.
- **Source:** [02 p. 175], [09 Eq 11–13], [10 Eq 14–16]
- **Verification:** T-IDC-01 · **Priority:** Must · **Status:** implemented · **Task:** —

### FR-IDC-005 — Incipient derivative of fractional order
For `p/q` with integers `p ≥ 0`, `1 ≤ q ≤ GIA_MAX_BRANCHES`, the kernel shall return all q branches of
`(φ')^{p/q} e^φ` as one state, branch m having argument `(arg φ' + 2πm)·p/q`.
- **Source:** [06 Eq 3.8–3.9], [02 Eq 14.7.2]
- **Verification:** T-IDC-02, VAL-06 · **Priority:** Must · **Status:** implemented · **Task:** —

### FR-IDC-006 — Second-order incipient LDE with variable coefficients
Given coefficient functions `a₁(t), a₀(t)` and initial conditions `f(0) = f₀`, `f̃'(0) = f₁`, the
kernel shall solve `f̃'' + a₁ f̃' + a₀ f = 0` as `f(t) = Σ cᵢ exp(∫₀ᵗ α̃ᵢ)`, where `α̃ᵢ(t)` are the roots
of `α̃² + a₁(t) α̃ + a₀(t) = 0` labelled by continuity in t (numerics N3), and `c₁ + c₂ = f₀`,
`c₁α̃₁(0) + c₂α̃₂(0) = f₁`; the derivative is termwise (§1). Where the discriminant vanishes on the
whole interval (a double root α̃), the solution set is the one family `c·e^{∫α̃}`: the kernel shall
solve it when `f₁ = α̃(0) f₀` and refuse otherwise. It shall not use [06 Eq 3.7], which solves the
equation under no reading (X12). Where the roots collide at an isolated time, it shall refuse.
- **Source:** [06 Eq 3.3–3.6], [09 Eq 3, 7], [10 Eq 8.1, 10.1]; PLAN R13, X12
- **Verification:** T-IDC-04, T-IDC-05, VAL-05 · **Priority:** Must · **Status:** implemented · **Task:** idc-lde2

### FR-IDC-007 — The binary function
Given constants `A, B` and the four initial conditions `f_σ(0)`, `f_σ^{(½)}(0)` for `σ ∈ {1, 2}`, the
kernel shall solve `f' + A f^{(½)} + B f = 0` as `f_σ = c_σ1 e^{u₁²t} + c_σ2 e^{u₂²t}`, with `u₁, u₂`
the roots of `u² + A u + B = 0`, half-derivative `Σ c_σi uᵢ e^{uᵢ²t}`, and constants from
`[[1, 1], [u₁, u₂]] c_σ = (f_σ(0), f_σ^{(½)}(0))`. For `u₁ = u₂` it shall refuse (no source defines it).
- **Source:** [02 Eq 14.7.1–14.7.7], [06 Eq 3.10–3.15]; PLAN R9, X7
- **Verification:** T-IDC-06 · **Priority:** Must · **Status:** implemented · **Task:** idc-binary

### FR-IDC-008 — Riccati equation by linearisation
Given `Q(t), R(t), P(t)` with `R ≠ 0` and `f(0) = f₀`, the kernel shall solve
`f' + Q f + R f² = P` by `f = y'/(R y)`, where y solves `R y'' − (R' − Q R) y' − P R² y = 0` as an
incipient LDE (FR-IDC-006) with `y(0) = 1`, `ỹ'(0) = R(0) f₀`, and shall report the traditional
Riccati residual of the result. The printed substitution of [06 Eq 3.17] shall not be used.
- **Source:** [06 Eq 3.16–3.18]; PLAN R11, X2
- **Verification:** T-IDC-07 · **Priority:** Should · **Status:** implemented · **Task:** idc-riccati

### FR-IDC-009 — The nonlinear equation of [02 Eq 14.10.1]
For constants `A, B`, the kernel shall return the two roots of `4u² + A u + B = 0` as the solutions
`F = e^{ut}` of `F·(d̃²/dt²)F² + A F²·(d̃/dt)F + B F³ = 0`.
- **Source:** [02 Eq 14.10.1–14.10.2]; X8
- **Verification:** T-IDC-08 · **Priority:** Could · **Status:** implemented · **Task:** idc-nonlinear-14-10

### FR-IDC-010 — Incipient Taylor projection
Given `f(t₀) ≠ 0`, `f'(t₀)`, a horizon Δ and an order `n ≥ 0`, the kernel shall compute
`f*(t₀+Δ) = f(t₀) Σ_{k=0}^{n} (aΔ)ᵏ/k!` with `a = f'(t₀)/f(t₀)`.
- **Source:** [09 Eq 10, 16–22], [10 Eq 13]; PLAN R10
- **Verification:** T-IDC-09, VAL-01 · **Priority:** Must · **Status:** implemented · **Task:** idc-taylor

### FR-IDC-011 — Solution drift on network trajectories
For a network whose flow matrix is constant, the kernel shall report solution drift as exactly zero
for every component and order 1–4, without computing it numerically. For any other network it shall
not report a solution drift: the sources define none for coupled nonlinear systems.
- **Source:** [06 §4 (i)]; PLAN E4b
- **Verification:** T-IDC-11, T-IDC-13 · **Priority:** Must · **Status:** implemented · **Task:** idc-drift-coupled

### FR-IDC-012 — Linearity in initial conditions
The solutions of FR-IDC-006 and FR-IDC-007 shall be linear in their initial conditions.
- **Source:** [06b §Scientific challenges ii]
- **Verification:** T-IDC-10 · **Priority:** Should · **Status:** implemented · **Task:** idc-lde2

### FR-IDC-013 — Refuse what the sources do not define
The kernel shall refuse: the direct Riccati duet form [06 Eq 3.22] (X3); Abel's n-et [06 Eq 3.25–3.27];
solution drift on a network with a non-constant flow matrix (FR-IDC-011).
- **Source:** PLAN §6
- **Verification:** T-IDC-13 · **Priority:** Must · **Status:** implemented · **Task:** idc-riccati

### FR-IDC-014 — Output-projection drift on network trajectories
For each component with `Q_i(t) ≠ 0`, the kernel shall report the second-order output-projection drift
over the output step Δ, `[09 Eq 13]` at k = 2: `(Q_i'' − (Q_i')²/Q_i) · Δ²/2`, with `Q_i'` the flow
balance and `Q_i''` its time derivative along the solved trajectory (numerics N7). It shall be labelled
as output-projection drift, distinct from FR-IDC-011.
- **Source:** [09 Eq 9–13], [10 Eq 12–16]
- **Verification:** T-IDC-12 · **Priority:** Should · **Status:** implemented · **Task:** idc-drift-coupled

### 2.2 Emergy algebra

### FR-EM-001 — Source emergy to output, and the split
Emergy entering a process shall be assigned to its output (rule 1); at a partition each leg shall carry
the emergy in proportion to its share of the flow, creating none (rule 3).
- **Source:** [02 p. 23 rules 1, 3; Eq 3.16–3.17]
- **Verification:** T-EM-01 · **Priority:** Must · **Status:** implemented · **Task:** —

### FR-EM-002 — Co-production
At a replicating process with n products, each product shall carry the whole input emergy, and the
emergy created shall equal `(n − 1)·Em(u)`.
- **Source:** [02 p. 23 rule 2; Eq 3.6–3.8]
- **Verification:** T-EM-01, T-EM-03 · **Priority:** Must · **Status:** implemented · **Task:** —

### FR-EM-003 — No double counting
Co-products from one process, when reunited, shall contribute the maximum and not the sum; emergy
returning on a feedback pathway shall not be re-injected; independent inputs shall still sum.
- **Source:** [02 p. 23 rule 4a, 4b]
- **Verification:** T-EM-02 · **Priority:** Must · **Status:** implemented · **Task:** —

### FR-EM-004 — Interaction in steady state
For an interaction whose inputs are all drawn, the product's emergy shall equal the sum of the drawn
inputs' emergy, so the equivalent source term `Φ(u₁, u₂)` is zero.
- **Source:** [02 Eq 3.9, 3.12, 3.15]; ADR 0017
- **Verification:** T-EM-04 · **Priority:** Must · **Status:** implemented · **Task:** emergy-source-terms

### FR-EM-005 — Ordinal forms of the three processes
The kernel shall construct co-production as a binary (two branches, each `Em(u)`), interaction as a
duet `[Em(u₁), Em(u₂)]`, and feedback as the duet-binary `[[a₁, a₂], [a₂, a₁]]`.
- **Source:** [22 Eq 6–8], [06b Eq 6, 10], [10 §Incipient Derivative]
- **Verification:** T-EM-05 · **Priority:** Should · **Status:** planned · **Task:** emergy-ordinal-forms

### FR-EM-006 — Circle product
The circle product shall keep every pair of factors in the outer arrangement
(`(a₁; a₂) ∘ [b₁, b₂] = [(a₁·b₁; a₂·b₁), (a₁·b₂; a₂·b₂)]` before reduction), and its cardinal reduction
shall map each pair to its product; `l ∘ l` shall be stored as `[l, l]` and reduce to `l²`.
- **Source:** [02 Eq 14.11.4–14.11.5], [06b Eq 2–3]; PLAN R12
- **Verification:** T-EM-06 · **Priority:** Should · **Status:** planned · **Task:** emergy-ordinal-forms

### FR-EM-007 — Global emergy balance arithmetic
Given input emergies, co-injection, co-production and re-normalisation factors and source terms, the
kernel shall evaluate the global balance of [02 Eq 3.21] and solve it for one unknown source term.
- **Source:** [02 Eq 3.18–3.26]
- **Verification:** T-EM-07, VAL-04 · **Priority:** Could · **Status:** implemented · **Task:** emergy-source-terms

### 2.3 Maximum Ordinality Principle

### FR-MOP-001 — First Fundamental Equation, one couple, affine-power boundary condition
For `β(t) = (a + b t)^p` and cardinality k, the kernel shall compute the solution of
`(d̃/dt)^k α = β` with `α(0) = 0`:
`α(t) = [((a+bt)^{(p+k)/k} − a^{(p+k)/k}) / (b(p+k))]^k` for `b ≠ 0`, and `α(t) = (β^{1/k} t / k)^k`
for `b = 0`, computed per numerics N1. The printed forms [23 Eq 5.5.7–5.5.8] shall not be used.
- **Source:** [23 Eq 5.4.2, 5.5.2–5.5.6]; PLAN R1, X1
- **Verification:** T-MOP-01, T-MOP-02, T-MOP-03, T-NUM-01 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### FR-MOP-002 — Domain of the First Equation
The kernel shall accept integer `k ≥ 1` with real or complex a, b, and rational non-integer k only with
real `a > 0`, `b ≥ 0`. It shall refuse any other k, `a + bt = 0` on `[0, t]`, and `b(p + k) = 0`
with `b ≠ 0`.
- **Source:** PLAN R1; numerics N1
- **Verification:** T-MOP-05 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### FR-MOP-003 — The Matrioska
For N components and a boundary condition per couple, the kernel shall solve each couple independently
and return the internal-representation Matrioska with a zero diagonal and `N(N−1)` entries; couples
without a boundary condition shall be marked unrelated, not zero.
- **Source:** [23 Eq 5.4.1–5.4.2, 5.6.1]
- **Verification:** T-MOP-04 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### FR-MOP-004 — Sampled boundary conditions
For β given as samples `(tᵢ, βᵢ)`, linearly interpolated, the kernel shall compute
`α(t) = {(1/k) ∫₀ᵗ β^{1/k}}^k` by adaptive quadrature with the branch of `β^{1/k}` continued along t
from the principal branch at `t = 0` (numerics N2), and shall refuse t outside the sampled range.
- **Source:** [23 Eq 5.5.5–5.5.7]; PLAN R1
- **Verification:** T-MOP-06 · **Priority:** Should · **Status:** planned · **Task:** mop-first-equation

### FR-MOP-005 — Second Fundamental Equation (printed solution)
Given `α₁₂(0)`, `c₁`, `c₂` and N, the kernel shall evaluate `A(t) = α₁₂(0) ∘ r + ln(c₁ + c₂t)`, the
specular `B(t) = [[A, −A], [−A, A]]`, and `{r} = e^{B(t)} ∘ (roots 13…1N)` of [23 Eq 6.1–6.3], with λ
null; it shall refuse `c₁ + c₂t ≤ 0` on `[0, t]`.
- **Source:** [23 Eq 6.1–6.3, §8 ii]; PLAN R4, R12
- **Verification:** T-MOP-07 · **Priority:** Must · **Status:** planned · **Task:** mop-second-equation

### FR-MOP-006 — The EQS operative form
Given the reference couple's coordinates `Σ₀(t), Φ₀(t), Θ₀(t)`, the factors `ψ₁,ᵢ`, `ψ₂`, ε₁, ε₂ = ε₃,
A and N, the kernel shall compute for each `l = 1…N−1` the coordinates of [23 Eq 7.1–7.5], with the
brackets evaluated as the relational product of the De Moivre root `B_l + C_l j + C_l k` and
`{Σ₀, Φ₀, Θ₀}` (FR-REL-001). It shall refuse ε₂ ≠ ε₃: [23 Eq 7.4] is defined only for equal angles.
- **Source:** [23 Eq 7.1–7.5]; PLAN R2, X11
- **Verification:** T-MOP-08, VAL-02 · **Priority:** Must · **Status:** planned · **Task:** mop-eqs

### FR-MOP-007 — Relational-valued couples
For k = 1 the kernel shall solve the First Equation componentwise on relational elements; for k > 1 it
shall refuse, because the sources define no division or power in the relational algebra.
- **Source:** PLAN §6
- **Verification:** T-MOP-09 · **Priority:** Should · **Status:** planned · **Task:** mop-first-equation

### FR-MOP-008 — Boundary conditions from the network
On request, the kernel shall set `e^{α_ij(t)}` to the empower carried from component i to component j
by the emergy pass, derive `β_ij` by [23 Eq 5.5.2] with k = 1, and mark couples with no pathway
unrelated.
- **Source:** PLAN R8
- **Verification:** T-MOP-10 · **Priority:** Should · **Status:** planned · **Task:** mop-network-beta

### 2.4 Relational Space algebra

### FR-REL-001 — The relational product, as printed
The kernel shall implement the bilinear product of relational elements with the table of
[23 Eq 5.1.3–5.1.5]: `i∘i = 1, i∘j = j, i∘k = k, j∘i = j, j∘j = −1, j∘k = k, k∘i = k, k∘j = k,
k∘k = −1`.
- **Source:** [23 Eq 5.1.3–5.1.5]; PLAN R2
- **Verification:** T-REL-01, VAL-02 · **Priority:** Must · **Status:** planned · **Task:** mop-relational-algebra

### FR-REL-002 — The De Moivre exponential
The kernel shall define `Exp{a·i ⊕ b·j ⊕ c·k} = eᵃ [cos ρ + (b·j + c·k) sin ρ / ρ]` with `ρ = √(b²+c²)`,
and the limit `eᵃ` at `ρ = 0`; it shall not use a power series.
- **Source:** [23 Eq 5.1.2]; PLAN R2
- **Verification:** T-REL-03 · **Priority:** Must · **Status:** planned · **Task:** mop-relational-algebra

### FR-REL-003 — Ordinal roots and their powers
The kernel shall construct the canonical ordinal root `r_{N,l}` (ε = 0, `√2·ψ = 2πl/(N−1)`) and raise
a root to an integer power by multiplying its angle, so that `r_{N,l}^{N−1} = 1`. It shall expose the
table-product power separately and shall not substitute one for the other.
- **Source:** [23 Eq A2.5–A2.6]; PLAN R3, X10
- **Verification:** T-REL-04 · **Priority:** Must · **Status:** planned · **Task:** mop-relational-algebra

### FR-REL-004 — No reassociation
Products of three or more relational elements shall be evaluated left to right in the order written;
the kernel shall not reassociate them.
- **Source:** [23 Eq 5.1.3–5.1.5] (non-associative); PLAN R2
- **Verification:** T-REL-02 · **Priority:** Must · **Status:** planned · **Task:** mop-relational-algebra

### 2.5 Ordinality and the generative step

### FR-ORD-001 — Couple classification
For each unordered couple of components the kernel shall assign: 2/2 if each reaches the other along
quantity-carrying legs (the ADR 0014 walk); otherwise 2 if both feed one interaction module; otherwise
½ if both are products of one replicating process; otherwise unrelated.
- **Source:** [22 Eq 6–8, 11.1], [10 §MOP]; PLAN R5, X6
- **Verification:** T-ORD-01 · **Priority:** Must · **Status:** planned · **Task:** mop-ordinality

### FR-ORD-002 — The Ordinality record
The kernel shall report `{k, n₂₂, n₂, n½, n_unrelated}`, with k the number of components (boundary
nodes and modules excluded) and the counts from FR-ORD-001.
- **Source:** [23 Eq 3.2], [22 Eq 11.1]; PLAN R5
- **Verification:** T-ORD-01 · **Priority:** Must · **Status:** planned · **Task:** mop-ordinality

### FR-ORD-003 — Maximum Ordinality
A network shall be at Maximum Ordinality exactly when every couple is 2/2, i.e. when its component
graph is strongly connected.
- **Source:** [22 §12.1, Eq 11.1]; PLAN R5
- **Verification:** T-ORD-02 · **Priority:** Must · **Status:** planned · **Task:** mop-ordinality

### FR-ORD-004 — Closure, labelled as a proxy
The fraction of components on a closed pathway shall be reported only as `closure`, labelled a proxy,
and shall decide nothing.
- **Source:** PLAN R5, E6
- **Verification:** T-ORD-03 · **Priority:** Must · **Status:** planned · **Task:** mop-ordinality

### FR-ORD-005 — The generative step under the Maximum Em-Power Principle
Below Maximum Ordinality, the kernel shall repeatedly add one `linear` pathway, of the seed's mean
weight, from a component in a sink strongly-connected component of the condensation to one in a source
component, choosing the candidate that maximises total empower at `t_end` with ties broken by
lexicographic `(from, to)`, until the network is at Maximum Ordinality; at Maximum Ordinality it shall
add nothing.
- **Source:** [02 Eq 5.3], [22 Eq 2, §12.1]; PLAN R6
- **Verification:** T-ORD-04 · **Priority:** Must · **Status:** planned · **Task:** mop-generative-empower

### FR-ORD-006 — Mode decided structurally
Whether a run was generative or functional shall be decided by comparing the output graph with the
seed, not by the input flag.
- **Source:** ADR 0011; TODO 10.4
- **Verification:** T-ORD-05 · **Priority:** Must · **Status:** implemented · **Task:** —

### FR-ORD-007 — Invariance under re-spelling
Reordering a model's nodes, edges or a module's legs shall not change any result.
- **Source:** ADR 0013, ADR 0015
- **Verification:** T-ORD-06 · **Priority:** Must · **Status:** implemented · **Task:** —

### 2.6 Harmony

### FR-HAR-001 — The harmony residual
For a Matrioska α with N ≥ 3, the kernel shall compute
`R_H = max_{j=2..N−1} |α_{1,j+1}/α_{12} − e^{2πi(j−1)/(N−1)}|` in ℂ, refusing `α₁₂ = 0`.
- **Source:** [23 Eq 5.6.5]; PLAN R7, X9
- **Verification:** T-HAR-02, T-HAR-03, T-HAR-04 · **Priority:** Must · **Status:** planned · **Task:** mop-harmony-detector

### FR-HAR-002 — Verdict for a construction
Given a construction (a map from boundary conditions to a Matrioska), the kernel shall evaluate it on a
harmonic input and on the fixed perturbations of vv-plan.md §6, and report `imposed` if `R_H` of every
output is within tolerance, `transported` if every output's `R_H` equals its input's, and `absent`
otherwise. For an observed Matrioska it shall report `present` or `absent`.
- **Source:** PLAN R7
- **Verification:** T-HAR-05, T-HAR-06, T-HAR-07 · **Priority:** Must · **Status:** planned · **Task:** mop-harmony-detector

### FR-HAR-003 — No harmony in the solver
The First and Second Equation solvers and the EQS shall not call the harmony constructor; the
detector shall not be linked against it.
- **Source:** PLAN §1 G-rules, B11
- **Verification:** T-HAR-01 · **Priority:** Must · **Status:** planned · **Task:** mop-harmony-detector

### FR-HAR-004 — The constructor, labelled
The existing harmony constructor shall live in `src/harmony.c` as `gia_harmony_assume_*`, and every
output it produces shall be labelled `assumed`.
- **Source:** PLAN E5, W8
- **Verification:** T-HAR-08 · **Priority:** Must · **Status:** planned · **Task:** mop-harmony-detector

### 2.7 Labelling and reporting

### FR-OUT-001 — Every output says what produced it
The run report shall list every output column and quantity with one label from `implemented`,
`classical`, `assumed`, `proxy`, `illustrative`.
- **Source:** PLAN §2 (E1–E7), BR-009
- **Verification:** T-OUT-02 · **Priority:** Must · **Status:** implemented · **Task:** mop-claims-remediation

### FR-OUT-002 — Refusals name their reason
Every refusal shall name the feature and the source reason (equation or PLAN §6 row), on stderr for
the CLI and through the `why` out-parameter for the library.
- **Source:** PLAN §6
- **Verification:** T-OUT-03 · **Priority:** Must · **Status:** planned · **Task:** mop-claims-remediation

### FR-OUT-003 — Harmony verdicts reported verbatim
The report shall carry one line per construction, `harmony.<construction>: <verdict>`, with the verdict
from FR-HAR-002.
- **Source:** PLAN R7
- **Verification:** T-OUT-02 · **Priority:** Must · **Status:** planned · **Task:** mop-harmony-detector

### 2.8 GSSK kernel

### FR-KER-001 — The kernel's method label
The GSSK kernel shall accept `"method": "expm"`, and shall accept `"method": "incipient"` as a
deprecated alias that produces byte-identical output and prints one notice that it is the matrix
exponential.
- **Source:** PLAN R15
- **Verification:** T-KER-01 · **Priority:** Should · **Status:** planned · **Task:** kernel-method-label

---

## 3. Non-functional requirements

### NFR-NUM-001 — Accuracy
Closed-form results shall have relative error ≤ 1e-12 where the condition number of the formula
(numerics.md) is ≤ 10; residuals of defining equations shall be ≤ 1e-6 relative; quadrature shall meet
a requested relative tolerance of 1e-10 or report failure.
- **Source:** vv-plan.md §3
- **Verification:** T-MOP-01, T-IDC-04, T-NUM-03 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### NFR-NUM-002 — No cancellation near t = 0
FR-MOP-001 shall retain relative error ≤ 1e-12 for `|b t / a|` from 1e-12 to 1, and shall be continuous
as b → 0.
- **Source:** numerics N1
- **Verification:** T-NUM-01, T-MOP-03 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### NFR-NUM-003 — Never emit a non-finite value
No library function shall return, and the CLI shall not write, an infinity or NaN; a result that would
overflow shall be refused with `GIA_E_RANGE`.
- **Source:** `AGENTS.md` §Fail-Safe; numerics N1, N6
- **Verification:** T-NUM-02 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### NFR-NUM-004 — One branch policy
Every multivalued function shall use the branch stated in numerics.md for it, and no other.
- **Source:** numerics.md §1
- **Verification:** T-MOP-05, T-MOP-06 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### NFR-NUM-005 — Quadrature reports its error
Every adaptive quadrature shall return an error estimate, and shall return `GIA_E_CONVERGENCE` when it
cannot meet its tolerance within its subdivision limit.
- **Source:** numerics N2
- **Verification:** T-NUM-03 · **Priority:** Must · **Status:** planned · **Task:** mop-first-equation

### NFR-NUM-006 — No clamping in the Giannantoni units
The Giannantoni units shall not clamp, floor or take absolute values of signed or complex coordinates
to keep them "physical"; `AGENTS.md`'s clamp applies to the GSSK kernel only.
- **Source:** PLAN §1 B6, G6
- **Verification:** T-NUM-04 · **Priority:** Must · **Status:** planned · **Task:** guard-agents-wording

### NFR-DET-001 — Bit-identical repeats
The same binary on the same input shall produce byte-identical output.
- **Source:** ADR 0013; BR-011
- **Verification:** T-DET-01 · **Priority:** Must · **Status:** implemented · **Task:** guard-reentrancy

### NFR-DET-002 — Agreement across toolchains
Every catalogue test shall pass, at its stated tolerance, on every CI toolchain (macOS clang, Linux
clang, Linux GCC); no result may depend on a libm difference beyond those tolerances.
- **Source:** vv-plan.md §4
- **Verification:** DEM-POR-01 · **Priority:** Must · **Status:** planned · **Task:** guard-reentrancy

### NFR-REE-001 — Reentrant
The Giannantoni units shall hold no file-scope or static mutable state; two models processed on two
threads shall give the same results as processed sequentially. (Defect at baseline: `src/engine.c:2830`
`static const gia_model *sort_model` feeds a `qsort` comparator.)
- **Source:** `AGENTS.md` §Architecture ("No global or static variables")
- **Verification:** T-REE-01, T-REE-02 · **Priority:** Must · **Status:** implemented · **Task:** guard-reentrancy

### NFR-MEM-001 — Memory
Every allocation shall have a matching free reachable from the public API; the Giannantoni tests shall
run clean under AddressSanitizer and LeakSanitizer.
- **Source:** `AGENTS.md` §Architecture
- **Verification:** T-MEM-01 · **Priority:** Must · **Status:** implemented · **Task:** guard-reentrancy

### NFR-ERR-001 — No process exits in the library
No library unit shall call `exit`, `abort` or `assert`-to-abort; every failure shall be a returned
`gia_status`.
- **Source:** IF-API-001
- **Verification:** T-ERR-01 · **Priority:** Must · **Status:** implemented · **Task:** guard-reentrancy

### NFR-LIM-001 — No silent truncation
Exceeding any fixed limit shall be refused with `GIA_E_LIMIT`. (Defect at baseline: `combine_inflows`
and `gia_emergy_at` in `src/engine.c` silently ignore a component's inflows beyond 64, and co-production
masks beyond 64 components.)
- **Source:** BR-009
- **Verification:** T-LIM-01 · **Priority:** Must · **Status:** implemented · **Task:** emergy-source-terms

### NFR-ROB-001 — Hostile seeds
The `mop` block parser shall reject malformed input with an error and never crash; it shall be fuzzed
with a committed corpus.
- **Source:** `AGENTS.md` §Fail-Safe ("JSON parsing must be strict")
- **Verification:** T-ROB-01, T-IN-02 · **Priority:** Should · **Status:** planned · **Task:** mop-first-equation

### NFR-POR-001 — Builds clean everywhere
The Giannantoni units shall compile without warnings under `-std=c99 -Wall -Wextra -Werror` with GCC,
Linux clang and Apple clang.
- **Source:** `AGENTS.md` §Tech Stack
- **Verification:** DEM-POR-01 · **Priority:** Must · **Status:** planned · **Task:** guard-reentrancy

### NFR-PERF-001 — Scale
The kernel shall solve the First Equation for N = 64 (4,032 couples) at 1,000 output times, affine-power
β, in under 2 s on a CI runner; N ≤ 64 is the supported limit (NFR-LIM-001 above it).
- **Source:** `src/engine.c` 64-wide masks; BR-011
- **Verification:** T-PERF-01 · **Priority:** Could · **Status:** planned · **Task:** mop-first-equation

### NFR-SEP-001 — Two engines stay separate
No Giannantoni unit shall include `gssk.h`.
- **Source:** ADR 0011
- **Verification:** INS-SEP-01 · **Priority:** Must · **Status:** implemented · **Task:** guard-trace

### NFR-TRC-001 — Traceability is enforced
`make check-trace` shall run in CI and fail on any of the conditions in README.md §The trace check.
- **Source:** PLAN §1 G5, G7; BR-010
- **Verification:** T-TRC-01 · **Priority:** Must · **Status:** implemented · **Task:** —

### NFR-COV-001 — Coverage of the Giannantoni units
Line coverage of the Giannantoni units shall be ≥ 90%, measured in CI; an unparsable figure shall fail.
- **Source:** PLAN §1 B3, G2
- **Verification:** T-COV-01 · **Priority:** Must · **Status:** implemented · **Task:** guard-coverage-giannantoni

### NFR-API-001 — No untested or stubbed API
Every function declared in `idc.h`, `mop.h` and `relational.h` shall be called by a test, and none
shall return a not-implemented code.
- **Source:** PLAN §1 B4, G3
- **Verification:** T-API-02 · **Priority:** Must · **Status:** implemented · **Task:** guard-api-called
