/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "engine/Lindhard.hpp"

#include "core/SincQuadrature.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <ranges>

namespace mosaiq {

namespace {

template<ScalarPhysical T>
[[nodiscard]] T branch_xi(Frequency<T> omega, WaveVector<T> q, T sign) noexcept
{
    return ((omega.raw() / q.raw()) + sign * q.raw()) / T{2};
}

template<ScalarPhysical T>
[[nodiscard]] T stable_fermi_log(T arg) noexcept
{
    if (arg >= T{0}) {
        return arg + std::log1p(std::exp(-arg));
    }
    return std::log1p(std::exp(arg));
}

template<ScalarPhysical T>
[[nodiscard]] T fermi_log_kernel(T kh, ReducedChemicalPotential<T> gamma, ReducedTemperature<T> tau) noexcept
{
    const T arg = (gamma.raw() - kh * kh) / tau.raw();
    return stable_fermi_log(arg);
}

/// Legacy high-T path: Σ f(kh) K with f = ln(1+e^{(γ−s²)/τ}), then × πτ/(4q).
template<ScalarPhysical T>
[[nodiscard]] T accumulate_cauchy_principal_value(T xi,
                                                   ReducedChemicalPotential<T> gamma,
                                                   ReducedTemperature<T> tau,
                                                   T h,
                                                   T pi_over_h,
                                                   long half_nodes) noexcept
{
    const auto indices = std::views::iota(-half_nodes, half_nodes + 1);
    return std::transform_reduce(
        std::ranges::begin(indices),
        std::ranges::end(indices),
        T{0},
        std::plus{},
        [&](long k) {
            const T kh = static_cast<T>(k) * h;
            const T f_kh = fermi_log_kernel(kh, gamma, tau);
            return f_kh * cauchy_principal_value_kernel(xi - kh, pi_over_h);
        });
}

/// Smooth residual only: Σ R_τ(γ − kh²) K(ξ − kh), with R_τ the Stratonovich softplus remainder.
template<ScalarPhysical T>
[[nodiscard]] T accumulate_smooth_residual_cpv(T xi,
                                               ReducedChemicalPotential<T> gamma,
                                               ReducedTemperature<T> tau,
                                               T h,
                                               T pi_over_h,
                                               long half_nodes) noexcept
{
    const T gamma_raw = gamma.raw();
    const T tau_raw = tau.raw();
    const auto indices = std::views::iota(-half_nodes, half_nodes + 1);
    return std::transform_reduce(
        std::ranges::begin(indices),
        std::ranges::end(indices),
        T{0},
        std::plus{},
        [&](long k) {
            const T kh = static_cast<T>(k) * h;
            const T residual = softplus_smooth_residual(gamma_raw - kh * kh, tau_raw);
            return residual * cauchy_principal_value_kernel(xi - kh, pi_over_h);
        });
}

template<ScalarPhysical T>
[[nodiscard]] T real_lindhard_kk_legacy(WaveVector<T> q,
                                        Frequency<T> omega,
                                        ReducedTemperature<T> tau,
                                        ReducedChemicalPotential<T> gamma) noexcept
{
    const T q_raw = q.raw();
    const T tau_raw = tau.raw();

    constexpr long half_nodes = static_cast<long>(constants::default_kk_sinc_nodes);
    const T h = std::sqrt(constants::sinc_strip_factor * constants::pi / static_cast<T>(half_nodes));
    const T pi_over_h = constants::pi / h;

    const T xi_n = branch_xi(omega, q, T{-1});
    const T xi_p = branch_xi(omega, q, T{1});

    const T cpv_n = accumulate_cauchy_principal_value(
        xi_n, gamma, tau, h, pi_over_h, half_nodes);
    const T cpv_p = accumulate_cauchy_principal_value(
        xi_p, gamma, tau, h, pi_over_h, half_nodes);

    return (cpv_n - cpv_p) * (constants::pi * tau_raw) / (T{4} * q_raw);
}

/// Reverse-Dedekind singularity excision: analytic Stratonovich step + smooth residual sinc.
template<ScalarPhysical T>
[[nodiscard]] T real_lindhard_kk_excised(WaveVector<T> q,
                                         Frequency<T> omega,
                                         ReducedTemperature<T> tau,
                                         ReducedChemicalPotential<T> gamma) noexcept
{
    const T q_raw = q.raw();
    const T gamma_raw = gamma.raw();

    constexpr long half_nodes = static_cast<long>(constants::default_kk_sinc_nodes);
    const T h = std::sqrt(constants::sinc_strip_factor * constants::pi / static_cast<T>(half_nodes));
    const T pi_over_h = constants::pi / h;

    const T xi_n = branch_xi(omega, q, T{-1});
    const T xi_p = branch_xi(omega, q, T{1});

    const T singular = analytic_stratonovich_re_contribution(xi_n, xi_p, gamma_raw, q_raw);

    if (!(tau.raw() > T{0})) {
        return singular;
    }

    const T residual_n = accumulate_smooth_residual_cpv(
        xi_n, gamma, tau, h, pi_over_h, half_nodes);
    const T residual_p = accumulate_smooth_residual_cpv(
        xi_p, gamma, tau, h, pi_over_h, half_nodes);

    return singular + (residual_n - residual_p) * constants::pi / (T{4} * q_raw);
}

}  // namespace

template<ScalarPhysical T>
T imaginary_lindhard(WaveVector<T> q,
                     Frequency<T> omega,
                     ReducedTemperature<T> tau,
                     ReducedChemicalPotential<T> gamma)
{
    const T q_raw = q.raw();
    const T tau_raw = tau.raw();
    const T gamma_raw = gamma.raw();

    const T nu_p = branch_xi(omega, q, T{1});
    const T nu_n = branch_xi(omega, q, T{-1});

    // τ ln(1 + e^{(γ − ν²)/τ}) via softplus: removes Inf/Inf NaN for τ ≲ 10^{-6}.
    const T soft_n = softplus_scaled(gamma_raw - nu_n * nu_n, tau_raw);
    const T soft_p = softplus_scaled(gamma_raw - nu_p * nu_p, tau_raw);

    return -(constants::pi / q_raw) * (soft_n - soft_p);
}

template<ScalarPhysical T>
T real_lindhard_kk(WaveVector<T> q,
                   Frequency<T> omega,
                   ReducedTemperature<T> tau,
                   ReducedChemicalPotential<T> gamma)
{
    if (tau.raw() < static_cast<T>(constants::singularity_excision_tau_threshold)) {
        return real_lindhard_kk_excised(q, omega, tau, gamma);
    }
    return real_lindhard_kk_legacy(q, omega, tau, gamma);
}

template<ScalarPhysical T>
LindhardResult<T> evaluate_lindhard(WaveVector<T> q,
                                    Frequency<T> omega,
                                    ReducedTemperature<T> tau,
                                    ReducedChemicalPotential<T> gamma)
{
    return LindhardResult<T>{
        real_lindhard_kk(q, omega, tau, gamma),
        imaginary_lindhard(q, omega, tau, gamma),
    };
}

template double imaginary_lindhard<double>(WaveVector<double>,
                                           Frequency<double>,
                                           ReducedTemperature<double>,
                                           ReducedChemicalPotential<double>);
template double real_lindhard_kk<double>(WaveVector<double>,
                                       Frequency<double>,
                                       ReducedTemperature<double>,
                                       ReducedChemicalPotential<double>);
template LindhardResult<double> evaluate_lindhard<double>(WaveVector<double>,
                                                          Frequency<double>,
                                                          ReducedTemperature<double>,
                                                          ReducedChemicalPotential<double>);

}  // namespace mosaiq
