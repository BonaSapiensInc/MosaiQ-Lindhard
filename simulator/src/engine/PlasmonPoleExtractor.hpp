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
#include "engine/RPA.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
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

/// Inputs for deterministic plasmon-pole extraction at fixed q.
template<ScalarPhysical T = double>
struct PlasmonPoleInputs {
    WaveVector<T> q{};
    PlasmonSpeciesParams<T> electron{};
    PlasmonSpeciesParams<T> ion{};
    BarePotentials<T> potentials{};
    T wigner_seitz_radius{};  ///< r_s [a_0] — Bohm–Gross bracket anchor
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

/// Stateless, deterministic plasmon-pole engine (MosaiQ-Lindhard Phase 6).
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

    /// Locate ω_p such that Re[ε(q, ω_p)] = 0; returns Im[ε] as Landau damping.
    template<ScalarPhysical T>
    [[nodiscard]] static std::optional<PlasmonState<T>> extract(
        PlasmonPoleInputs<T> inputs,
        PlasmonPolePolicy<T> policy = {}) noexcept;

    /// Evaluate ε(q, ω) for the two-component RPA stack at electron frequency ω_e.
    template<ScalarPhysical T>
    [[nodiscard]] static std::complex<T> evaluate_epsilon(PlasmonPoleInputs<T> inputs,
                                                          Frequency<T> omega_e) noexcept;
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

}  // namespace mosaiq
