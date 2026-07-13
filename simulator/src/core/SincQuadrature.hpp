/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#pragma once

#include "core/Concepts.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <functional>
#include <vector>

namespace mosaiq {

/// Parameters for the modified sinc rule (Davis–Rabinowitz Eq. 3.4.6.21; legacy N=256, h = sqrt(7 pi / N)).
template<ScalarPhysical T = double>
struct SincQuadratureParams {
    std::size_t half_nodes{constants::default_gv_sinc_nodes};
    T strip_factor{constants::sinc_strip_factor};

    [[nodiscard]] constexpr T step() const noexcept
    {
        return static_cast<T>(std::sqrt(strip_factor * constants::pi /
                                        static_cast<T>(half_nodes)));
    }

    [[nodiscard]] constexpr std::size_t node_count() const noexcept
    {
        return 2 * half_nodes + 1;
    }
};

/// Arguments for the GV Eq. (4.44) grand-canonical real Lindhard sinc integrand.
template<ScalarPhysical T = double>
struct GvLindhardSincArgs {
    WaveVector<T> q;
    Frequency<T> omega;
    ReducedTemperature<T> tau;
    ReducedChemicalPotential<T> gamma;
};

/// One sinc abscissa sample: index k, mapped x = asinh(e^{kh}), and quadrature weight factor.
template<ScalarPhysical T = double>
struct SincSample {
    long k{};
    T kh{};
    T x{};
    T weight{};  // 1 / sqrt(1 + exp(-2 kh)) in legacy get_GVrealLindhard
};

/// Callable integrand at a single sinc node (pure; no mutation).
template<ScalarPhysical T = double>
using GvLindhardIntegrand = std::function<T(const SincSample<T>&, const GvLindhardSincArgs<T>&)>;

/// Declarations only — implementations deferred to SincQuadrature.cpp (Phase 1).
template<ScalarPhysical T = double>
class SincQuadrature {
public:
    explicit SincQuadrature(SincQuadratureParams<T> params = {});

    [[nodiscard]] SincQuadratureParams<T> params() const noexcept;

    /// Lazy view of sinc indices k ∈ [-N, N] (legacy for-loop body).
    [[nodiscard]] auto node_indices() const;

    /// Map index k → SincSample (asinh map + legacy weight).
    [[nodiscard]] SincSample<T> sample_at(long k) const;

    /// Fold integrand contributions over all nodes via `<ranges>` (Phase 1 implementation).
    [[nodiscard]] T integrate(GvLindhardIntegrand<T> integrand,
                            const GvLindhardSincArgs<T>& args) const;

    /// GV grand-canonical real Lindhard: legacy `get_GVrealLindhard` pipeline entry point.
    [[nodiscard]] T evaluate_gv_real_lindhard(const GvLindhardSincArgs<T>& args) const;

    /// Collect all samples materialized (testing / debugging only).
    [[nodiscard]] std::vector<SincSample<T>> all_samples() const;

private:
    SincQuadratureParams<T> params_;
};

// --- free-function pipeline hooks (declaration only) ---

template<ScalarPhysical T>
[[nodiscard]] T asinh_sinc_abscissa(T kh) noexcept;

template<ScalarPhysical T>
[[nodiscard]] T sinc_weight_factor(T kh) noexcept;

template<ScalarPhysical T>
[[nodiscard]] GvLindhardIntegrand<T> make_gv_lindhard_integrand();

/// Cauchy principal-value kernel [1 - cos(π Δ/h)] / (π Δ/h) with L'Hôpital regularization at grid nodes.
template<ScalarPhysical T>
[[nodiscard]] T cauchy_principal_value_kernel(T diff, T pi_over_h) noexcept
{
    constexpr T resonance_tolerance = T{1.0e-13};
    if (std::abs(diff) <= resonance_tolerance) {
        // lim_{Δ→0} [1 - cos(cΔ)] / (cΔ) = 0 (L'Hôpital); avoids IEEE 754 0/0 → NaN traps.
        return T{0};
    }
    const T argument = pi_over_h * diff;
    return (T{1} - std::cos(argument)) / argument;
}

/// Softplus identity: τ ln(1 + e^{x/τ}) = (x)_+ + τ ln(1 + e^{-|x|/τ}).
/// At x = 0 the residual is τ ln 2, enforcing the Stratonovich midpoint f(0) = 1/2
/// of the reverse-Dedekind (T → 0) Fermi–Dirac cut.
template<ScalarPhysical T>
[[nodiscard]] T softplus_scaled(T x, T tau) noexcept
{
    if (!(tau > T{0})) {
        return (x > T{0}) ? x : T{0};
    }
    const T abs_x = std::abs(x);
    const T singular = (x > T{0}) ? x : T{0};
    // Beyond ~40 τ the residual is smaller than double underflow relative to O(1).
    constexpr T residual_cutoff = T{40};
    if (abs_x > residual_cutoff * tau) {
        return singular;
    }
    return singular + tau * std::log1p(std::exp(-abs_x / tau));
}

/// Smooth residual alone: R_τ(x) = τ ln(1 + e^{-|x|/τ}) (Stratonovich value τ ln 2 at x = 0).
template<ScalarPhysical T>
[[nodiscard]] T softplus_smooth_residual(T x, T tau) noexcept
{
    if (!(tau > T{0})) {
        return T{0};
    }
    const T abs_x = std::abs(x);
    constexpr T residual_cutoff = T{40};
    if (abs_x > residual_cutoff * tau) {
        return T{0};
    }
    return tau * std::log1p(std::exp(-abs_x / tau));
}

/// Continuum Hilbert content of the Stratonovich step kernel g(s) = (γ − s²)_+:
/// Ĩ(ξ) = P∫_{-a}^{a} (γ − s²)/(ξ − s) ds with a = √γ (sinc-compatible sign).
template<ScalarPhysical T>
[[nodiscard]] T stratonovich_step_hilbert_kernel(T xi, T gamma) noexcept
{
    if (!(gamma > T{0})) {
        return T{0};
    }
    const T a = std::sqrt(gamma);
    const T diff_minus = a - xi;
    const T diff_plus = a + xi;
    // Removable edge singularity: (γ − ξ²) ln|(a − ξ)/(a + ξ)| → 0 as ξ → ±a,
    // since γ − ξ² = (a − ξ)(a + ξ) and u ln|u| → 0.
    constexpr T edge_tolerance = T{1.0e-14};
    if (std::abs(diff_minus) <= edge_tolerance || std::abs(diff_plus) <= edge_tolerance) {
        return T{2} * a * xi;
    }
    const T log_ratio = std::log(std::abs(diff_minus) / std::abs(diff_plus));
    // Ĩ(ξ) = −[(γ − ξ²) ln|(a − ξ)/(a + ξ)| − 2 a ξ]
    return -((gamma - xi * xi) * log_ratio - T{2} * a * xi);
}

/// Analytic singular contribution to Re χ̃^L from the reverse-Dedekind step:
/// (π/4q) [H̃[g](ξ_−) − H̃[g](ξ_+)] with H̃ the sinc-compatible Hilbert transform.
template<ScalarPhysical T>
[[nodiscard]] T analytic_stratonovich_re_contribution(T xi_minus,
                                                      T xi_plus,
                                                      T gamma,
                                                      T q) noexcept
{
    const T i_minus = stratonovich_step_hilbert_kernel(xi_minus, gamma);
    const T i_plus = stratonovich_step_hilbert_kernel(xi_plus, gamma);
    return (i_minus - i_plus) / (T{4} * q);
}

/// Stable Fermi–Dirac occupation 1/(e^y + 1) with Stratonovich value 1/2 at y = 0.
template<ScalarPhysical T>
[[nodiscard]] T fermi_dirac_occupation_stable(T y) noexcept
{
    if (y >= T{0}) {
        const T e = std::exp(-y);
        return e / (T{1} + e);
    }
    return T{1} / (T{1} + std::exp(y));
}

}  // namespace mosaiq
