# MosaiQ-Lindhard: C++20 Deterministic Physics Sandbox

**Status:** Architecture Blueprint (Phase 1 — two-fermion plasma)  
**Theory Reference:** `manuscript/two-fermi.tex` — *Linear response representation of two-fermion plasmas*  
**Docs map:** [`README.md`](README.md) · Phase~2 extension: [`BEREZINIAN_ARCHITECTURE.md`](BEREZINIAN_ARCHITECTURE.md)  
**Purpose:** Bare-metal validation of the finite-temperature Lindhard function and RPA structure factor using the Kramers–Krönig / sinc-quadrature pathway formalized in Section IV, bypassing the analytically defective Matsubara pole-summation routes documented in Appendices F and G.

---

## 1. Mission Statement

The MosaiQ-Lindhard engine is an **absolute zero-point physics sandbox**: a C++20 system that evaluates thermodynamic and linear-response quantities with deterministic, auditable numerics. It exists to prove that the manuscript's central claim — that Hilbert-transform evaluation via sinc-quadrature is the only rigorous route when contour integration and infinite Matsubara sums fail — survives contact with silicon, not just LaTeX.

Design priorities, in order:

1. **Determinism** — identical inputs yield bit-identical outputs across platforms (within documented floating-point policy).
2. **Type-safe physics** — wave vectors, frequencies, and temperatures are distinct types; dimensional aliasing is a compile error.
3. **Zero hidden state** — numerical pipelines are pure functions over immutable views; no silent mutation across integration steps.
4. **Traceability** — every published curve in the manuscript must be reproducible from a single `main` invocation with a pinned parameter file.

---

## 2. Directory Structure Blueprint

```
simulator/
├── CMakeLists.txt          # C++20 project root; strict warnings-as-errors
├── src/
│   ├── core/
│   │   ├── SincQuadrature.hpp / .cpp   # Cauchy principal-value & Hilbert kernels
│   │   ├── RootFinder.hpp              # Newton–Raphson (chemical potential only)
│   │   ├── BrentRootFinder.hpp         # Bracketed Brent root solver (Phase 6)
│   │   └── FermiDirac.hpp / .cpp       # Occupation functions & half-integer FD integrals
│   ├── physics/
│   │   ├── Constants.hpp               # constexpr atomic Hartree, ℏ, k_B, …
│   │   └── ChemicalPotential.hpp / .cpp
│   ├── engine/
│   │   ├── Lindhard.hpp / .cpp         # χ⁰(q, ω; T) real & imaginary parts
│   │   ├── RPA.hpp / .cpp              # Two-component coupling, ε(q, ω), S(q, ω)
│   │   └── PlasmonPoleExtractor.hpp / .cpp  # Re[ε]=0 root finding, ω_p(q) (Phase 6)
│   └── main.cpp                        # CLI entry: parameter ingest → pipeline → output
└── tests/
    ├── test_sinc_quadrature.cpp        # Error bound regression (target: 10⁻²⁵⁹ class)
    ├── test_fermi_dirac.cpp
    ├── test_lindhard.cpp               # Cross-check vs. known limits (T→0, q→0)
    ├── test_rpa.cpp
    └── test_plasmon_pole_extractor.cpp # Brent bracket + Re[ε]=0 regression (Phase 6)
```

### Build Contract (`CMakeLists.txt`)

| Requirement | Policy |
|-------------|--------|
| Standard | `CMAKE_CXX_STANDARD 20`, extensions off |
| Warnings | `-Wall -Wextra -Wpedantic -Werror` (MSVC: `/W4 /WX`) |
| Optimization | `-O3 -march=native` for release; `-O0 -g` for debug |
| Tests | CTest + GoogleTest or Catch2 (TBD at Phase 1 kickoff) |
| I/O | Header-only JSON or TOML for run manifests; no runtime network |
| Editor / clangd | Each of `src/`, `src/core/`, `src/physics/`, `src/engine/`, `tests/` must contain `compile_flags.txt` (`-std=c++20`, include path to `src/`). New subdirectories copy this file or Problems-tab diagnostics regress. |

---

## 3. C++20 Design Philosophy — Weapons of Deterministic Control

### 3.1 Topological Type Enforcement via `concept`

Physical quantities carry their **semantic dimension** at compile time. A `WaveVector` cannot be passed where a `Frequency` is expected, even if both are stored as `double` internally.

```cpp
template<typename T>
concept ScalarPhysical = std::is_floating_point_v<T>;

template<ScalarPhysical T>
struct WaveVector { T q; /* … */ };

template<ScalarPhysical T>
struct Frequency { T omega; /* … */ };

// Lindhard API — dimensionally typed surface
template<ScalarPhysical T>
struct LindhardResult {
    T real_part;
    T imag_part;
};

template<ScalarPhysical T>
[[nodiscard]] LindhardResult<T>
evaluate_lindhard(WaveVector<T> q, Frequency<T> omega, Temperature<T> tau);
```

**Goal:** block the class of bugs that produced unphysical aliasing in legacy Matsubara/Euler–Maclaurin implementations — wrong units, swapped arguments, and illegal order-of-operations interchanges — before link time.

### 3.2 Zero-Cost Abstraction via `constexpr`

Fundamental constants and closed-form limits are **burned into the translation unit** at compile time:

```cpp
namespace mosaiq::constants {
    inline constexpr double hartree_to_ev = /* … */;
    inline constexpr double pi            = 3.141592653589793238462643383279502884;
    // Atomic Hartree energy, ℏ, m_e, … — single source of truth
}
```

Thermodynamic identities that reduce to rational expressions (e.g., Sommerfeld expansions, T→0 limits of Fermi–Dirac step functions) live in `constexpr` functions so the optimizer can fold them away in hot loops while unit tests still evaluate them at compile time via `static_assert`.

### 3.3 Functional Purity via `<ranges>`

Numerical pipelines compose over **immutable range views**; accumulators never mutate shared buffers:

```cpp
auto sinc_weights(int N, double h) -> std::vector<double>;

template<std::ranges::range R>
[[nodiscard]] double integrate_sinc(R&& samples, /* kernel params */) {
    return std::ranges::fold_left(
        std::forward<R>(samples),
        0.0,
        [](double acc, auto term) { return acc + term; }  // pure fold
    );
}
```

Infinite series are expressed as lazy views with explicit truncation policies (`take_while`, `take(N)`) so truncation order and cutoff are visible in the call graph — a direct countermeasure to the illegal Fubini interchanges and divergent `∑ j^{-1/2}` tails identified in the manuscript.

---

## 4. Module Responsibilities

### 4.1 `core/`

| Module | Responsibility | Manuscript Anchor |
|--------|----------------|-------------------|
| **SincQuadrature** | Cauchy principal-value and Hilbert transform via modified sinc rule; pole correction terms | Eq. (sinc-quadrature-Cauchy-principal-value-integral), practical rule |
| **RootFinder** | Newton–Raphson for implicit equations with analytic derivatives | Chemical potential inversion only |
| **BrentRootFinder** | Derivative-free bracketed root finding (Brent) | Phase 6 plasmon poles — **no** numerical derivatives |
| **FermiDirac** | `n(ε)`, half-integer FD integrals, Sommerfeld helpers | Thermodynamic potential chapter |

**Acceptance criterion (Phase 1):** sinc quadrature of a benchmark integrand matches reference to ≤ 10⁻²⁵⁹ relative error at documented `(N, h)` — matching the error class quoted in Section IV.

#### Sinc Quadrature Singularity Handling (L'Hôpital Regularization)

To compute the real part of the Lindhard function, the engine employs a Sinc quadrature rule for the Cauchy Principal Value integral. At extremely high momenta ($q \ge 15.0$), kinematic variables like $\xi_-$ frequently align perfectly with the integration grid nodes (e.g., $s=0.0$). This alignment triggers a $0/0$ grid resonance trap in the integration kernel. To prevent IEEE 754 `NaN` propagation, the engine explicitly applies L'Hôpital's rule: $\lim_{\Delta \to 0} [1 - \cos(c\Delta)] / (c\Delta) = 0.0$. This architectural choice completely avoids technical debt such as arbitrary grid-offsetting, ensuring absolute determinism and numerical safety across all temperature and momentum scales.

### 4.2 `physics/`

| Module | Responsibility |
|--------|----------------|
| **Constants** | Unit system, conversions, constexpr limits |
| **ChemicalPotential** | Invert `n(μ)` constraint for target density; feeds Lindhard `γ`, `ξ±` |

### 4.3 `engine/`

| Module | Responsibility | Manuscript Anchor |
|--------|----------------|-------------------|
| **Lindhard** | `Re χ^L(q, ω; τ)` via Kramers–Krönig / Hilbert of `Im χ^L`; **no** Matsubara pole sum | Section IV, Eq. (real-Lindhard-Kramers-Kronig) |
| **RPA** | Two-species susceptibility coupling, dielectric function `ε(q, ω)`, structure factor `S(q, ω)` | Main text RPA sections |
| **PlasmonPoleExtractor** | Deterministic extraction of `ω_p(q)` from `Re ε = 0`; simultaneous Landau damping `Im ε(q, ω_p)` | Section IV causality constraint; collective-mode dispersion |

**Explicit non-goals:** Gouedard–Deutsch parametric sums, Kim–Lee Euler–Maclaurin hybrids — retained only as negative regression fixtures proving divergence, not production code paths.

---

## 5. Data Flow

```
┌─────────────┐     ┌──────────────────┐     ┌─────────────────┐
│  Constants  │────▶│ ChemicalPotential │────▶│    Lindhard     │
│  (constexpr)│     │  (RootFinder)     │     │ (SincQuadrature)│
└─────────────┘     └──────────────────┘     └────────┬────────┘
                                                       │
                                                       ▼
                                              ┌─────────────────┐
                                              │      RPA        │
                                              │ ε(q, ω), χ^RPA  │
                                              └────────┬────────┘
                                                       │
                              ┌────────────────────────┴────────────────────────┐
                              ▼                                                 ▼
                     ┌─────────────────┐                           ┌──────────────────────┐
                     │  main / output  │                           │ PlasmonPoleExtractor │
                     │  (CSV, JSON)    │                           │  Re[ε]=0 → ω_p(q)    │
                     └─────────────────┘                           └──────────────────────┘
```

All arrows denote **pure function calls**; no module holds process-global mutable state.

---

## 6. Phased Execution Roadmap

### Phase 1 — Math Core

**Deliverables**

- `SincQuadrature` with modified/practical sinc rules and pole correction
- `RootFinder` with bracket validation and monotonicity guards
- Unit tests: sinc error ≤ 10⁻²⁵⁹ on canonical integrands; root finder on known transcendental equations

**Exit gate:** CTest green; benchmark log archived in `docs/benchmarks/phase1.md`.

### Phase 2 — Thermodynamics

**Deliverables**

- `FermiDirac` half-integer integrals
- `ChemicalPotential` solver for target `n` at finite `τ`
- Tests: Sommerfeld limit, step-function limit at `τ → 0`

**Exit gate:** chemical potential inversion reproduces manuscript table values within tolerance.

### Phase 3 — The Lindhard

**Deliverables**

- `Lindhard::evaluate` — imaginary part (closed or quadrature) + real part via Hilbert/sinc
- Tests: compare `Re/Im χ^L` against Section IV limiting cases; **regression test** showing Matsubara-sum path blows up or fails tolerance

**Exit gate:** full `(q, ω, τ)` grid matches reference data generated from manuscript formulas (not legacy pole code).

### Phase 4 — RPA Coupling

**Deliverables**

- `RPA` two-component susceptibility matrix inversion
- Structure factor extraction `S(q, ω)` and static `S(q)`
- Integration test: two-fermion plasma scenario from `two-fermi.tex`

**Exit gate:** end-to-end `main` run produces publication-ready curves with pinned seed and manifest.

### Phase 5 — CLI Structure-Factor Export

**Deliverables**

- `mosaiq_simulator` CLI: `(q, ω)` grid sweep over two-component RPA
- Multi-channel output: `S_ee`, `S_ii`, `S_ei` with fluctuation–dissipation consistency
- Deterministic output manifest under `output/`

**Exit gate:** structure-factor grids reproducible from pinned `(r_s, T)` invocation.

### Phase 6 — Deterministic Plasmon Pole Extraction

**Deliverables**

- `PlasmonPoleExtractor` — Brent-based root finding on `Re ε(q, ω) = 0`
- `BrentRootFinder` — pure, stateless bracketed solver (derivative-free)
- `test_plasmon_pole_extractor` — bracket validation, root residual, Landau damping extraction

**Exit gate:** `|Re ε(q, ω_p)|` ≤ solver tolerance at documented reference `(r_s, T, q)`; dispersion curve `ω_p(q)` monotonic in the plasmon branch.

See **Section 10** for the full Phase 6 design specification.

---

## 7. Testing & Integrity Policy

| Layer | Scope |
|-------|-------|
| Unit | Each `core/` and `physics/` function in isolation |
| Integration | `Lindhard` fed by live `ChemicalPotential`; `PlasmonPoleExtractor` fed by live `RPA` |
| Regression | Known-bad Matsubara/Euler–Maclaurin paths must **fail** documented assertions |
| Phase 6 | `test_plasmon_pole_extractor`: Brent bracket + `Re ε = 0` residual at `ω_p` |
| Reproducibility | Every release tag stores `{compiler, flags, git SHA, manifest hash}` |

Tests live under `simulator/tests/` and are registered with CTest. No test may depend on network or wall-clock randomness.

---

## 8. Relationship to the Manuscript

| Manuscript Element | Simulator Module |
|--------------------|------------------|
| Section IV — Kramers–Krönig / Hilbert | `engine/Lindhard`, `core/SincQuadrature` |
| Sinc-quadrature error bound 10⁻²⁵⁹ | `tests/test_sinc_quadrature.cpp` |
| Appendix F — contour / Jordan breakdown | Documented **anti-pattern**; not implemented |
| Appendix G — Kim–Lee hybrid | Documented **anti-pattern**; not implemented |
| Two-component RPA | `engine/RPA` |
| Collective-mode dispersion `ω_p(q)` | `engine/PlasmonPoleExtractor` |
| Causality-constrained real/imaginary pair | Kramers–Krönig pipeline; pole extraction inherits `Im ε` at root |

The simulator is the executable counterpart to the manuscript's claim: **deterministic sinc-quadrature is the only production pathway**; pole-sum shortcuts are mathematically excluded by construction. Phase 6 extends this philosophy to collective-mode extraction: the plasmon dispersion is computed entirely within the **strictly constrained physical space** defined by exact causality, without derivative pollution or Matsubara detours.

---

## 10. Phase 6: Deterministic Plasmon Pole Extraction

Phase 6 closes the loop between the manuscript's causality-constrained Lindhard evaluation and the observable collective-mode spectrum of a partially degenerate two-component plasma. The `PlasmonPoleExtractor` module locates the plasmon dispersion trajectory `ω_p(q)` by root finding on the dielectric function assembled from the Phase 3–4 pipeline. Every operation is pure, stateless, and auditable: the extractor holds no hidden mutable state and produces `std::optional<PlasmonState>` so that failed brackets or non-convergent roots are reported explicitly rather than masked by silent fallbacks.

### 10.1 Design Philosophy — Strictly Constrained Physical Space

The manuscript establishes that the real and imaginary parts of the susceptibility must satisfy the Kramers–Krönig relation as a **hard constraint**, not a post-hoc consistency check. Phase 6 applies the same principle to collective-mode identification:

- The objective function `Re[ε(q, ω)]` is evaluated exclusively through the sinc-quadrature Hilbert transform of the convergent `Im χ^L`, followed by two-component RPA assembly of `ε(q, ω)`.
- No Matsubara pole summation, no contour-residue shortcuts, and no ensemble-mixed analytic continuations enter the evaluation path.
- A located root `ω_p` is accepted only if it lies within a **validated bracket** on the causal frequency axis; otherwise the extractor returns `std::nullopt`.

This maps plasmon dispersion extraction into the same strictly constrained physical space as the Lindhard engine itself: the root is a property of the causality-respecting dielectric function, not of an auxiliary numerical regularization.

### 10.2 The Objective — Roots of `Re[ε(q, ω)] = 0`

Collective excitations in the RPA dielectric formalism correspond to poles of the resolvent `(ε(q, ω))^{-1}`. In the weak-damping limit, the plasmon dispersion manifold is identified by the condition

\[
\mathrm{Re}\,\varepsilon(q, \omega_p) = 0 ,
\]

at fixed wave vector `q`. The Phase 6 engine therefore reduces dispersion extraction to a one-dimensional root-finding problem in `ω` at each `q`, with the dielectric function supplied by:

```
χ^L_e, χ^L_i  ←  evaluate_lindhard (Kramers–Krönig / sinc)
ε(q, ω)       ←  evaluate_dielectric (two-component RPA determinant)
```

The C++20 surface constrains the scalar objective through the `DielectricFunction` concept:

```cpp
template<typename F>
concept DielectricFunction = std::invocable<F, double> && requires(F f, double omega) {
    { f(omega) } -> std::convertible_to<double>;
};
```

Any callable satisfying this concept may be passed to the bracket and root stages; the production implementation binds `f(ω) = Re[ε(q, ω)]` via `PlasmonPoleExtractor::evaluate_epsilon`.

The output type `PlasmonState` records:

| Field | Meaning |
|-------|---------|
| `q` | Wave vector (strong-typed `WaveVector`) |
| `omega_p` | Root frequency `ω_p` (`Frequency`) |
| `landau_damping` | `Im[ε(q, ω_p)]` — dissipation at the pole |

### 10.3 Derivative-Free Determinism — Brent's Method

**Newton–Raphson and all derivative-based root finders are explicitly excluded from the plasmon extraction path.**

Rationale: the sinc-quadrature engine achieves a documented Cauchy principal-value error class of **O(10⁻²⁵⁹)** on the Hilbert transform that pins `Re χ^L`. Numerical differentiation of `Re ε(q, ω)` with respect to `ω` would superimpose finite-difference noise on this pristine objective, destroying the error budget that the manuscript treats as non-negotiable. The chemical-potential solver (`RootFinder.hpp`) retains Newton–Raphson because its target function admits stable analytic derivatives; the dielectric objective does not justify the same assumption at production quadrature order.

Instead, Phase 6 implements **Brent's method** in `core/BrentRootFinder.hpp`:

| Mechanism | Role |
|-----------|------|
| Bisection | Guaranteed bracket convergence |
| Secant step | Acceleration when interpolation is unreliable |
| Inverse quadratic interpolation | Superlinear convergence near simple roots |

Properties enforced by design:

1. **Pure** — `brent_root(f, a, b, policy)` is a free function; no object holds iteration state between calls.
2. **Deterministic** — fixed bracket, fixed iteration cap, fixed tolerance; identical inputs yield identical outputs.
3. **Derivative-free** — only evaluations of `f(ω)` are required, each one a full deterministic pass through the Lindhard → RPA stack.

The root finder accepts any callable satisfying `BracketedRootFunction` (`f(x) → double`) with a validated sign change on `[ω_low, ω_high]`.

### 10.4 Physics-Informed Bracketing — Bohm–Gross Anchor

Bracketed root finding is only as reliable as its initial interval. An unconstrained frequency scan risks locking onto spurious sign changes associated with non-plasmon structure, numerical noise, or — conceptually — the unphysical resonances that arise when Matsubara/contour shortcuts contaminate the dielectric function. Phase 6 therefore anchors the search space with the **Bohm–Gross long-wavelength dispersion relation**.

In electron reduced units (`q̄ = q/k_F`, `ω̄ = ω/E_F`, `τ = k_B T/E_F`), the anchor frequency is

\[
\bar{\omega}_\mathrm{BG}^2 \approx \bar{\omega}_p^2 + 6\,\tau\,\bar{q}^2 ,
\]

where `ω̄_p = ω_p / E_F` is the reduced plasma frequency derived from the target density (via `r_s`). The routine `PlasmonPoleExtractor::bohm_gross_frequency` implements this estimate; `compute_initial_bracket` then:

1. Performs a **deterministic forward scan** from `ω_floor`, with geometric step growth, until `Re ε` changes sign.
2. Caps the scan at `scan_ceiling_factor × ω̄_BG` to confine the search to the plasmon sector suggested by the long-wavelength limit.
3. Falls back to symmetric expansion about `ω̄_BG` if the forward scan does not encounter a sign change.

This two-stage policy prevents the root finder from wandering into unphysical regions of the `(q, ω)` plane while remaining fully deterministic: no random restarts, no adaptive heuristics that depend on wall-clock or platform.

### 10.5 Simultaneous Landau Damping Extraction

Landau damping is not computed in a separate post-processing pipeline. Once Brent's method converges to `ω_p`, the engine evaluates the **same** complex dielectric function `ε(q, ω_p)` that supplied the real-part objective:

\[
\gamma_L(q) \equiv \mathrm{Im}\,\varepsilon(q, \omega_p) .
\]

Because `Im χ^L` is obtained from the closed-form grand-canonical expression and `Re χ^L` from the Hilbert transform of that imaginary part, the pair `(Re ε, Im ε)` at the root inherits the manuscript's causality constraint by construction. The Landau damping stored in `PlasmonState::landau_damping` is therefore the exact imaginary response at the plasmon pole — not a finite-difference estimate, not a separate Matsubara channel, and not an ensemble-mixed analytic continuation.

Operationally:

```cpp
const auto state = PlasmonPoleExtractor::extract(inputs);
// state->omega_p          : Re[ε(q, ω_p)] ≈ 0
// state->landau_damping   : Im[ε(q, ω_p)]
```

### 10.6 Module Map and Acceptance Criteria

| File | Responsibility |
|------|----------------|
| `core/BrentRootFinder.hpp` | Generic Brent solver; `BracketedRootFunction` concept |
| `engine/PlasmonPoleExtractor.hpp` | `DielectricFunction` concept, bracket policy, `extract()` API |
| `engine/PlasmonPoleExtractor.cpp` | Bohm–Gross anchor (`bohm_gross_frequency`) |
| `tests/test_plasmon_pole_extractor.cpp` | Brent sanity check, end-to-end `Re ε → 0` regression |

**Acceptance criteria**

- `|Re ε(q, ω_p)| ≤ 10⁻⁸` at reference thermodynamic state (documented `(r_s, T, q)`).
- Bracket construction succeeds without manual `ω` tuning for the electron plasmon branch at `q ≥ 0.1`.
- Repeated invocations with identical `PlasmonPoleInputs` yield bit-identical `PlasmonState` (Release build, pinned compiler flags).

**Explicit non-goals:** Newton–Raphson on `Re ε`; finite-difference damping rates; Matsubara-pole or contour-based dielectric evaluations; stochastic global optimization.

---

## 9. Open Decisions (to resolve at Phase 1 kickoff)

1. Floating-point policy: `double` everywhere vs. `long double` in sinc accumulation only.
2. Test framework: GoogleTest vs. Catch2.
3. Output format: CSV vs. JSON for figure reproduction scripts.
4. Parallelism: single-threaded default; optional `std::execution` for grid sweeps in Phase 4.

---

*Document version: 1.1 — Phase 6 plasmon pole extraction (`PlasmonPoleExtractor`, `BrentRootFinder`) aligned with `manuscript/two-fermi.tex` causality constraint.*
