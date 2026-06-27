/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "core/FermiDirac.hpp"

#include "physics/Constants.hpp"

#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <ranges>

namespace mosaiq {

namespace {

template<ScalarPhysical T>
[[nodiscard]] T order_exponent(FermiDiracOrder order) noexcept
{
    switch (order) {
    case FermiDiracOrder::Half:
        return T{0.5};
    case FermiDiracOrder::NegativeHalf:
        return T{-0.5};
    }
    return T{0};
}

template<ScalarPhysical T>
[[nodiscard]] T trapezoidal_contribution(T order, T eta, long N, T h) noexcept
{
    const auto integrand = [&](long n) -> T {
        const T x = static_cast<T>(n) * h;
        return std::pow(x, T{2} * order + T{1}) / (std::exp(x * x - eta) + T{1});
    };

    const T endpoint = integrand(0);
    const auto inner_nodes = std::views::iota(1L, N);

    const T inner_sum = std::transform_reduce(
        std::ranges::begin(inner_nodes),
        std::ranges::end(inner_nodes),
        T{0},
        std::plus{},
        integrand);

    return T{2} * h * (endpoint + inner_sum);
}

template<ScalarPhysical T>
[[nodiscard]] T pole_correction_contribution(T order, T eta, T h, long k_max) noexcept
{
    const T pi = constants::pi;
    const T pi_square = pi * pi;
    const T two_pi = constants::two_pi;
    const T eta_square = eta * eta;

    const bool is_negative_eta = std::signbit(eta);
    const T pole_prefactor = is_negative_eta ? T{4} * pi : T{-4} * pi;

    const auto pole_term = [&](long k) -> T {
        const T dk = static_cast<T>(k);
        const T oddk = T{2} * dk + T{1};
        const T rho_k = std::sqrt(eta_square + pi_square * (oddk * oddk));

        const T theta_k = is_negative_eta ? pi + std::atan(oddk * pi / eta)
                                          : std::atan(oddk * pi / eta);

        const T delta_k = (two_pi / h) * std::sqrt((rho_k + eta) / T{2});
        const T lambda_k = std::exp((two_pi / h) * std::sqrt((rho_k - eta) / T{2}));

        const T numerator =
            lambda_k * std::sin(delta_k + order * theta_k) - std::sin(order * theta_k);
        const T denominator =
            T{1} + lambda_k * lambda_k - T{2} * lambda_k * std::cos(delta_k);

        return std::pow(rho_k, order) * (numerator / denominator);
    };

    const auto pole_indices = std::views::iota(0L, k_max + 1);
    const T pole_sum = std::transform_reduce(
        std::ranges::begin(pole_indices),
        std::ranges::end(pole_indices),
        T{0},
        std::plus{},
        pole_term);

    return pole_prefactor * pole_sum;
}

template<ScalarPhysical T>
[[nodiscard]] long pole_count(T eta, T h) noexcept
{
    const T pi = constants::pi;
    const T d = h * constants::sinc_strip_factor;
    const T d_square = d * d;
    const T eta_square = eta * eta;
    const T two_d_square_plus_eta = T{2} * d_square + eta;
    const T two_d_square_plus_eta_square = two_d_square_plus_eta * two_d_square_plus_eta;

    const T k_max_square =
        (two_d_square_plus_eta_square - eta_square) / (T{4} * pi * pi);

    if (!(k_max_square > T{0})) {
        return 0;
    }

    return static_cast<long>(std::sqrt(k_max_square));
}

}  // namespace

template<ScalarPhysical T>
T fermi_dirac_integral(FermiDiracOrder order_kind, FermiDiracEta<T> eta_arg)
{
    const T order = order_exponent<T>(order_kind);
    const T eta = eta_arg.raw();

    const T x0 = std::signbit(eta) ? T{8} : std::sqrt(eta + T{64});
    constexpr long N = 128;
    const T h = x0 / static_cast<T>(N);

    const T trapezoid = trapezoidal_contribution(order, eta, N, h);
    const long k_max = pole_count(eta, h);
    const T pole = pole_correction_contribution(order, eta, h, k_max);

    return trapezoid + pole;
}

template double fermi_dirac_integral<double>(FermiDiracOrder, FermiDiracEta<double>);

}  // namespace mosaiq
