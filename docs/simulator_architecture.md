# SOPHUS-Lindhard: C++20 Deterministic Physics Sandbox

**Status:** Architecture Blueprint  
**Theory Reference:** `manuscript/two-fermi.tex` — *Classical linear response representation of two-fermion plasmas*  
**Purpose:** Bare-metal validation of the finite-temperature Lindhard function and RPA structure factor using the Kramers–Krönig / sinc-quadrature pathway formalized in Section IV, bypassing the analytically defective Matsubara pole-summation routes documented in Appendices F and G.

---

## 1. Mission Statement

The SOPHUS-Lindhard engine is an **absolute zero-point physics sandbox**: a C++20 system that evaluates thermodynamic and linear-response quantities with deterministic, auditable numerics. It exists to prove that the manuscript's central claim — that Hilbert-transform evaluation via sinc-quadrature is the only rigorous route when contour integration and infinite Matsubara sums fail — survives contact with silicon, not just LaTeX.

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
│   │   ├── RootFinder.hpp / .cpp       # Bracketed root solvers (chemical potential)
│   │   └── FermiDirac.hpp / .cpp       # Occupation functions & half-integer FD integrals
│   ├── physics/
│   │   ├── Constants.hpp               # constexpr atomic Hartree, ℏ, k_B, …
│   │   └── ChemicalPotential.hpp / .cpp
│   ├── engine/
│   │   ├── Lindhard.hpp / .cpp         # χ⁰(q, ω; T) real & imaginary parts
│   │   └── RPA.hpp / .cpp              # Two-component coupling, S(q, ω), HNC/RPA hooks
│   └── main.cpp                        # CLI entry: parameter ingest → pipeline → output
└── tests/
    ├── test_sinc_quadrature.cpp        # Error bound regression (target: 10⁻²⁵⁹ class)
    ├── test_fermi_dirac.cpp
    ├── test_chemical_potential.cpp
    ├── test_lindhard.cpp               # Cross-check vs. known limits (T→0, q→0)
    └── test_rpa.cpp
```

### Build Contract (`CMakeLists.txt`)

| Requirement | Policy |
|-------------|--------|
| Standard | `CMAKE_CXX_STANDARD 20`, extensions off |
| Warnings | `-Wall -Wextra -Wpedantic -Werror` (MSVC: `/W4 /WX`) |
| Optimization | `-O3 -march=native` for release; `-O0 -g` for debug |
| Tests | CTest + GoogleTest or Catch2 (TBD at Phase 1 kickoff) |
| I/O | Header-only JSON or TOML for run manifests; no runtime network |

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
namespace sophus::constants {
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
| **RootFinder** | Robust bracketing (Brent/ bisection) for implicit equations | Chemical potential inversion |
| **FermiDirac** | `n(ε)`, half-integer FD integrals, Sommerfeld helpers | Thermodynamic potential chapter |

**Acceptance criterion (Phase 1):** sinc quadrature of a benchmark integrand matches reference to ≤ 10⁻²⁵⁹ relative error at documented `(N, h)` — matching the error class quoted in Section IV.

### 4.2 `physics/`

| Module | Responsibility |
|--------|----------------|
| **Constants** | Unit system, conversions, constexpr limits |
| **ChemicalPotential** | Invert `n(μ)` constraint for target density; feeds Lindhard `γ`, `ξ±` |

### 4.3 `engine/`

| Module | Responsibility | Manuscript Anchor |
|--------|----------------|-------------------|
| **Lindhard** | `Re χ^L(q, ω; τ)` via Kramers–Krönig / Hilbert of `Im χ^L`; **no** Matsubara pole sum | Section IV, Eq. (real-Lindhard-Kramers-Kronig) |
| **RPA** | Two-species susceptibility coupling, structure factor `S(q, ω)`, static limit checks | Main text RPA sections |

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
                                              │ (2-component)   │
                                              └────────┬────────┘
                                                       │
                                                       ▼
                                              ┌─────────────────┐
                                              │  main / output  │
                                              │  (CSV, JSON)    │
                                              └─────────────────┘
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

---

## 7. Testing & Integrity Policy

| Layer | Scope |
|-------|-------|
| Unit | Each `core/` and `physics/` function in isolation |
| Integration | `Lindhard` fed by live `ChemicalPotential` |
| Regression | Known-bad Matsubara/Euler–Maclaurin paths must **fail** documented assertions |
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

The simulator is the executable counterpart to the manuscript's claim: **deterministic sinc-quadrature is the only production pathway**; pole-sum shortcuts are mathematically excluded by construction.

---

## 9. Open Decisions (to resolve at Phase 1 kickoff)

1. Floating-point policy: `double` everywhere vs. `long double` in sinc accumulation only.
2. Test framework: GoogleTest vs. Catch2.
3. Output format: CSV vs. JSON for figure reproduction scripts.
4. Parallelism: single-threaded default; optional `std::execution` for grid sweeps in Phase 4.

---

*Document version: 1.0 — initial architecture blueprint aligned with `manuscript/two-fermi.tex`.*
