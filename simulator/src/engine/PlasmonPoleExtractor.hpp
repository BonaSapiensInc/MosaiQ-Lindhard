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

#include "core/BrentRootFinder.hpp"
#include "core/Concepts.hpp"
#include "engine/Lindhard.hpp"
#include "engine/PolyLogRPA.hpp"
#include "engine/RPA.hpp"
#include "engine/ZetaRPA.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <limits>
#include <optional>

namespace mosaiq {

/// Real-part dielectric objective Re[ε(q, ω)] for bracketed root finding at fixed q.
template<typename F>
concept DielectricFunction = std::invocable<F, double> && requires(F f, double omega) {
    { f(omega) } -> std::convertible_to<double>;
};

/// Bracket [ω_low, ω_high] for Re[ε(q, ω)] = 0.
template<ScalarPhysical T = double>
struct FrequencyBracket {
    Frequency<T> low{};
    Frequency<T> high{};

    [[nodiscard]] constexpr T width() const noexcept { return high.raw() - low.raw(); }
};

/// Extracted plasmon pole on the Re[ε] = 0 manifold.
template<ScalarPhysical T = double>
struct PlasmonState {
    WaveVector<T> q{};
    Frequency<T> omega_p{};
    T landau_damping{};  ///< Im[ε(q, ω_p)] — dissipation at the pole

    [[nodiscard]] constexpr T q_raw() const noexcept { return q.raw(); }
    [[nodiscard]] constexpr T omega_raw() const noexcept { return omega_p.raw(); }
};

/// Per-species thermodynamic state for two-component ε(q, ω) evaluation.
template<ScalarPhysical T = double>
struct PlasmonSpeciesParams {
    ReducedTemperature<T> tau{};
    ReducedChemicalPotential<T> gamma{};
    T fermi_energy{};  ///< E_F [Ha] — maps electron ω to ion ω via ω_i = ω_e E_F,e / E_F,i
    T mass{constants::electron_mass_amu};
};

/// Inputs for deterministic plasmon-pole extraction at fixed q (StandardRPA / two-component).
template<ScalarPhysical T = double>
struct PlasmonPoleInputs {
    WaveVector<T> q{};
    PlasmonSpeciesParams<T> electron{};
    PlasmonSpeciesParams<T> ion{};
    BarePotentials<T> potentials{};
    T wigner_seitz_radius{};  ///< r_s [a_0] — Bohm–Gross bracket anchor
};

/// Opt-in **scalar** dielectric root inputs (Phase Z2 diagnostic; isolated from manuscript
/// two-component RPA). Supports `StandardRPA` (ε = 1 − v χ^L), PolyLog-RPA (ε^(s)), and
/// Zeta pathways (ε^ζ).
template<ScalarPhysical T = double>
struct PlasmonPoleZetaInputs {
    WaveVector<T> q{};
    ReducedTemperature<T> tau{};
    ReducedChemicalPotential<T> gamma{};  ///< chemical potential γ = μ/E_F
    T bare_potential{};                   ///< scalar v(q)
    T wigner_seitz_radius{};
    CouplingRegime<T> regime{};
    ResponsePathway pathway{ResponsePathway::ZetaRPA};
    bool force_pathway{true};
    BorweinPolicy borwein{};
    ZetaWeightParameters weight_params{};
};

template<ScalarPhysical T = double>
struct PlasmonBracketPolicy {
    T relative_half_width{T{0.20}};       ///< fallback symmetric half-width as fraction of ω_BG
    T absolute_floor{T{1.0e-4}};         ///< minimum scan step in reduced ω units
    T expansion_factor{T{1.618033988749895}};  ///< geometric step growth during scan
    std::size_t expansion_max{32};        ///< symmetric bracket expansions (fallback)
    std::size_t scan_max{512};           ///< deterministic forward scan budget
    T omega_floor{T{1.0e-4}};           ///< lower frequency bound [E_F units]
    T scan_ceiling_factor{T{3.0}};       ///< scan limit = scan_ceiling_factor × ω_BG
};

template<ScalarPhysical T = double>
struct PlasmonPolePolicy {
    BrentRootFinderPolicy<T> root_finder{};
    PlasmonBracketPolicy<T> bracket{};
};

/// Stateless, deterministic plasmon-pole engine (MosaiQ-Lindhard Phase 6 + Phase Z2 scalar Zeta).
struct PlasmonPoleExtractor final {
    PlasmonPoleExtractor() = delete;

    /// Bohm–Gross long-wavelength estimate in electron reduced units:
    /// ω̄² ≈ ω̄_p² + 6 τ q̄² with ω̄_p = ω_p / E_F.
    template<ScalarPhysical T>
    [[nodiscard]] static T bohm_gross_frequency(WaveVector<T> q,
                                                ReducedTemperature<T> tau,
                                                T rs,
                                                T mass = constants::electron_mass_amu) noexcept;

    /// Physics-informed bracket: forward scan from ω_floor with Bohm–Gross ceiling,
    /// then symmetric expansion fallback.
    template<ScalarPhysical T>
    [[nodiscard]] static std::optional<FrequencyBracket<T>> compute_initial_bracket(
        WaveVector<T> q,
        ReducedTemperature<T> tau,
        T rs,
        DielectricFunction auto epsilon_real,
        PlasmonBracketPolicy<T> policy = {}) noexcept;

    /// Locate ω_p such that Re[ε(q, ω_p)] = 0 for **two-component StandardRPA**.
    /// Bit-stable manuscript pathway — do not alter semantics.
    template<ScalarPhysical T>
    [[nodiscard]] static std::optional<PlasmonState<T>> extract(
        PlasmonPoleInputs<T> inputs,
        PlasmonPolePolicy<T> policy = {}) noexcept;

    /// Locate ω_p such that Re[ε_scalar(q, ω_p)] = 0 (StandardRPA or Zeta; opt-in diagnostic).
    template<ScalarPhysical T>
    [[nodiscard]] static std::optional<PlasmonState<T>> extract(
        PlasmonPoleZetaInputs<T> inputs,
        PlasmonPolePolicy<T> policy = {}) noexcept;

    /// Evaluate ε(q, ω) for the two-component RPA stack at electron frequency ω_e.
    template<ScalarPhysical T>
    [[nodiscard]] static std::complex<T> evaluate_epsilon(PlasmonPoleInputs<T> inputs,
                                                          Frequency<T> omega_e) noexcept;

    /// Evaluate scalar dielectric: StandardRPA → 1 − v χ^L; Zeta → ε^ζ via locked W_ζ.
    template<ScalarPhysical T>
    [[nodiscard]] static std::complex<T> evaluate_epsilon_scalar(PlasmonPoleZetaInputs<T> inputs,
                                                                 Frequency<T> omega) noexcept;

    /// Evaluate scalar ε^ζ(q, ω) via locked W_ζ Lindhard dressing (Zeta pathways only).
    template<ScalarPhysical T>
    [[nodiscard]] static std::complex<T> evaluate_epsilon_zeta(PlasmonPoleZetaInputs<T> inputs,
                                                               Frequency<T> omega) noexcept;
};

// --- template definitions (header-only hot path) ---

template<ScalarPhysical T>
std::optional<FrequencyBracket<T>> PlasmonPoleExtractor::compute_initial_bracket(
    WaveVector<T> q,
    ReducedTemperature<T> tau,
    T rs,
    DielectricFunction auto epsilon_real,
    PlasmonBracketPolicy<T> policy) noexcept
{
    const T center = bohm_gross_frequency(q, tau, rs);
    if (!std::isfinite(center) || center <= T{0}) {
        return std::nullopt;
    }

    const T scan_ceiling = center * policy.scan_ceiling_factor;
    T omega = policy.omega_floor;
    T f_prev = epsilon_real(omega);
    if (!std::isfinite(f_prev)) {
        return std::nullopt;
    }

    T step = std::max(policy.absolute_floor, center * policy.relative_half_width / T{4});

    for (std::size_t scan = 0; scan < policy.scan_max; ++scan) {
        const T omega_next = std::min(omega + step, scan_ceiling);
        const T f_next = epsilon_real(omega_next);

        if (std::isfinite(f_next) && f_prev * f_next < T{0}) {
            return FrequencyBracket<T>{
                Frequency<T>{omega},
                Frequency<T>{omega_next},
            };
        }

        if (omega_next >= scan_ceiling) {
            break;
        }

        omega = omega_next;
        f_prev = f_next;
        step = std::min(step * policy.expansion_factor, scan_ceiling - omega);
    }

    // Fallback: symmetric expansion about the Bohm–Gross center.
    T half_width = std::max(policy.absolute_floor, policy.relative_half_width * center);

    for (std::size_t attempt = 0; attempt <= policy.expansion_max; ++attempt) {
        const T low_raw = std::max(policy.omega_floor, center - half_width);
        const T high_raw = std::min(scan_ceiling, center + half_width);

        const T f_low = epsilon_real(low_raw);
        const T f_high = epsilon_real(high_raw);

        if (std::isfinite(f_low) && std::isfinite(f_high) && f_low * f_high < T{0}) {
            return FrequencyBracket<T>{
                Frequency<T>{low_raw},
                Frequency<T>{high_raw},
            };
        }

        half_width *= policy.expansion_factor;
    }

    return std::nullopt;
}

template<ScalarPhysical T>
std::complex<T> PlasmonPoleExtractor::evaluate_epsilon(PlasmonPoleInputs<T> inputs,
                                                       Frequency<T> omega_e) noexcept
{
    const T omega_e_raw = omega_e.raw();
    const T scale = inputs.electron.fermi_energy / inputs.ion.fermi_energy;
    const T omega_i_raw = omega_e_raw * scale;

    const T k_f_e = std::sqrt(T{2} * inputs.electron.mass * inputs.electron.fermi_energy);
    const T k_f_i = std::sqrt(T{2} * inputs.ion.mass * inputs.ion.fermi_energy);
    const T q_i_raw = inputs.q.raw() * (k_f_e / k_f_i);

    const LindhardResult<T> chi_e = evaluate_lindhard(
        inputs.q,
        omega_e,
        inputs.electron.tau,
        inputs.electron.gamma);

    const LindhardResult<T> chi_i = evaluate_lindhard(
        WaveVector<T>{q_i_raw},
        Frequency<T>{omega_i_raw},
        inputs.ion.tau,
        inputs.ion.gamma);

    return evaluate_dielectric(chi_e, chi_i, inputs.potentials);
}

template<ScalarPhysical T>
std::complex<T> PlasmonPoleExtractor::evaluate_epsilon_scalar(PlasmonPoleZetaInputs<T> inputs,
                                                              Frequency<T> omega) noexcept
{
    if (inputs.pathway == ResponsePathway::StandardRPA) {
        const LindhardResult<T> chi_l =
            evaluate_lindhard(inputs.q, omega, inputs.tau, inputs.gamma);
        const std::complex<T> one{T{1}, T{0}};
        return one - inputs.bare_potential * as_complex(chi_l);
    }
    if (inputs.pathway == ResponsePathway::PolyLogRPA) {
        const LindhardResult<T> chi_l =
            evaluate_lindhard(inputs.q, omega, inputs.tau, inputs.gamma);
        PolyLogRpaInputs<double> pl_in{
            .q = WaveVector<double>{static_cast<double>(inputs.q.raw())},
            .omega = Frequency<double>{static_cast<double>(omega.raw())},
            .chi_lindhard =
                LindhardResult<double>{
                    static_cast<double>(chi_l.real()),
                    static_cast<double>(chi_l.imag()),
                },
            .bare_potential = static_cast<double>(inputs.bare_potential),
            .regime =
                CouplingRegime<double>{
                    static_cast<double>(inputs.regime.rs),
                    static_cast<double>(inputs.regime.gamma_plasma),
                    static_cast<double>(inputs.regime.tau),
                },
            .weight_params = inputs.weight_params,
        };
        const auto pl = evaluate_polylog_rpa(pl_in);
        if (!pl) {
            return {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
        }
        return {static_cast<T>(pl->epsilon.real()), static_cast<T>(pl->epsilon.imag())};
    }
    return evaluate_epsilon_zeta(inputs, omega);
}

template<ScalarPhysical T>
std::complex<T> PlasmonPoleExtractor::evaluate_epsilon_zeta(PlasmonPoleZetaInputs<T> inputs,
                                                            Frequency<T> omega) noexcept
{
    if (inputs.pathway != ResponsePathway::ZetaRPA &&
        inputs.pathway != ResponsePathway::ZetaRPA_Experimental) {
        return {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
    }

    const LindhardResult<T> chi_l =
        evaluate_lindhard(inputs.q, omega, inputs.tau, inputs.gamma);

    ZetaRpaInputs<double> zeta_in{
        .q = WaveVector<double>{static_cast<double>(inputs.q.raw())},
        .omega = Frequency<double>{static_cast<double>(omega.raw())},
        .chi_lindhard =
            LindhardResult<double>{
                static_cast<double>(chi_l.real()),
                static_cast<double>(chi_l.imag()),
            },
        .bare_potential = static_cast<double>(inputs.bare_potential),
        .regime =
            CouplingRegime<double>{
                static_cast<double>(inputs.regime.rs),
                static_cast<double>(inputs.regime.gamma_plasma),
                static_cast<double>(inputs.regime.tau),
            },
        .borwein = inputs.borwein,
        .pathway = inputs.pathway,
        .force_pathway = inputs.force_pathway,
        .weight_params = inputs.weight_params,
    };

    const auto zeta = evaluate_zeta_rpa(zeta_in);
    if (!zeta) {
        return {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
    }
    return {static_cast<T>(zeta->epsilon.real()), static_cast<T>(zeta->epsilon.imag())};
}

template<ScalarPhysical T>
std::optional<PlasmonState<T>> PlasmonPoleExtractor::extract(
    PlasmonPoleInputs<T> inputs,
    PlasmonPolePolicy<T> policy) noexcept
{
    const auto epsilon_real = [inputs](double omega_e) noexcept -> double {
        const auto epsilon =
            PlasmonPoleExtractor::evaluate_epsilon(inputs, Frequency<T>{static_cast<T>(omega_e)});
        return epsilon.real();
    };

    const auto bracket = compute_initial_bracket(
        inputs.q,
        inputs.electron.tau,
        inputs.wigner_seitz_radius,
        epsilon_real,
        policy.bracket);

    if (!bracket) {
        return std::nullopt;
    }

    const auto omega_root = brent_root<T>(
        epsilon_real,
        bracket->low.raw(),
        bracket->high.raw(),
        policy.root_finder);

    if (!omega_root) {
        return std::nullopt;
    }

    const std::complex<T> epsilon_at_root =
        evaluate_epsilon(inputs, Frequency<T>{*omega_root});

    return PlasmonState<T>{
        inputs.q,
        Frequency<T>{*omega_root},
        epsilon_at_root.imag(),
    };
}

template<ScalarPhysical T>
std::optional<PlasmonState<T>> PlasmonPoleExtractor::extract(
    PlasmonPoleZetaInputs<T> inputs,
    PlasmonPolePolicy<T> policy) noexcept
{
    const auto epsilon_real = [inputs](double omega) noexcept -> double {
        const auto epsilon = PlasmonPoleExtractor::evaluate_epsilon_scalar(
            inputs, Frequency<T>{static_cast<T>(omega)});
        return epsilon.real();
    };

    const auto bracket = compute_initial_bracket(
        inputs.q,
        inputs.tau,
        inputs.wigner_seitz_radius,
        epsilon_real,
        policy.bracket);

    if (!bracket) {
        return std::nullopt;
    }

    // Scale the bracketed objective so Brent's absolute tolerance / IQI path remains
    // well-conditioned when |Re ε| ≫ 1 (strong-coupling / amplified-v scalar diagnostics).
    // Two-component StandardRPA extract is intentionally untouched.
    const T f_low = static_cast<T>(epsilon_real(static_cast<double>(bracket->low.raw())));
    const T f_high = static_cast<T>(epsilon_real(static_cast<double>(bracket->high.raw())));
    const T scale = std::max(std::abs(f_low), std::abs(f_high));
    if (!(scale > T{0}) || !std::isfinite(scale)) {
        return std::nullopt;
    }

    const auto epsilon_scaled = [epsilon_real, scale](double omega) noexcept -> double {
        return epsilon_real(omega) / static_cast<double>(scale);
    };

    const auto omega_root = brent_root<T>(
        epsilon_scaled,
        bracket->low.raw(),
        bracket->high.raw(),
        policy.root_finder);

    if (!omega_root) {
        return std::nullopt;
    }

    // Unscaled bisection polish: steep strong-coupling crossings can leave
    // |Re ε| ≳ 10^{-8} even when the scaled Brent abscissa has converged.
    T lo = bracket->low.raw();
    T hi = bracket->high.raw();
    T f_lo = f_low;
    T f_hi = f_high;
    T omega = *omega_root;
    for (std::size_t polish = 0; polish < 64; ++polish) {
        const T re = static_cast<T>(epsilon_real(static_cast<double>(omega)));
        if (std::abs(re) <= T{1.0e-8}) {
            break;
        }
        const T mid = (lo + hi) / T{2};
        const T f_mid = static_cast<T>(epsilon_real(static_cast<double>(mid)));
        if (f_lo * f_mid <= T{0}) {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
        omega = (lo + hi) / T{2};
        if (!(std::abs(hi - lo) > T{0})) {
            break;
        }
    }

    const std::complex<T> epsilon_at_root =
        evaluate_epsilon_scalar(inputs, Frequency<T>{omega});
    if (!std::isfinite(epsilon_at_root.real()) || !std::isfinite(epsilon_at_root.imag())) {
        return std::nullopt;
    }
    if (std::abs(epsilon_at_root.real()) > T{1.0e-8}) {
        return std::nullopt;
    }

    return PlasmonState<T>{
        inputs.q,
        Frequency<T>{omega},
        epsilon_at_root.imag(),
    };
}

}  // namespace mosaiq
