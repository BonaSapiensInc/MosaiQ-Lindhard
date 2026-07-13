/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "core/Concepts.hpp"
#include "core/SincQuadrature.hpp"
#include "engine/Lindhard.hpp"
#include "physics/Constants.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

using namespace mosaiq;

WaveVector<> reference_q{1.0};
Frequency<> reference_omega{0.5};
ReducedTemperature<> reference_tau{0.05};
ReducedChemicalPotential<> reference_gamma{1.0};

/// Verbatim legacy `get_KKrealLindhard` reference (mutable loop) for pipeline regression.
double reference_legacy_kk_real(double q, double omega, double tau, double gamma)
{
    constexpr long half_nodes = static_cast<long>(constants::default_kk_sinc_nodes);
    const double h = std::sqrt(constants::sinc_strip_factor * constants::pi / static_cast<double>(half_nodes));
    const double pi_over_h = constants::pi / h;

    const double xi_n = ((omega / q) - q) / 2.0;
    const double xi_p = ((omega / q) + q) / 2.0;

    double cpv_n = 0.0;
    double cpv_p = 0.0;
    for (long k = -half_nodes; k <= half_nodes; ++k) {
        const double kh = static_cast<double>(k) * h;
        const double arg4exp = (gamma - kh * kh) / tau;
        const double f_kh = std::log(1.0 + std::exp(arg4exp));

        const double delta_n = xi_n - kh;
        if (std::abs(delta_n) > 1.0e-13) {
            cpv_n += f_kh * ((1.0 - std::cos(pi_over_h * delta_n)) / (pi_over_h * delta_n));
        }

        const double delta_p = xi_p - kh;
        if (std::abs(delta_p) > 1.0e-13) {
            cpv_p += f_kh * ((1.0 - std::cos(pi_over_h * delta_p)) / (pi_over_h * delta_p));
        }
    }

    return (cpv_n - cpv_p) * (constants::pi * tau) / (4.0 * q);
}

void test_imaginary_lindhard_finite()
{
    const double imag = imaginary_lindhard(
        reference_q, reference_omega, reference_tau, reference_gamma);

    assert(std::isfinite(imag));
    assert(!std::isnan(imag));
    assert(!std::isinf(imag));
    assert(imag < 0.0);
    assert(std::abs(imag) < 10.0);

    std::cout << "Im χ^L (GV, q=1, w=0.5, tau=0.05, gamma=1) = " << imag << '\n';
}

void test_evaluate_lindhard_bundle()
{
    const LindhardResult<> chi = evaluate_lindhard(
        reference_q, reference_omega, reference_tau, reference_gamma);

    assert(std::isfinite(chi.real()));
    assert(std::isfinite(chi.imag()));
    assert(chi.imag() < 0.0);
    assert(chi.real() == real_lindhard_kk(
               reference_q, reference_omega, reference_tau, reference_gamma));
    assert(chi.imag() == imaginary_lindhard(
               reference_q, reference_omega, reference_tau, reference_gamma));

    std::cout << "Re χ^L (KK) = " << chi.real() << ", Im χ^L (GV) = " << chi.imag() << '\n';
}

void test_kk_vs_gv_real_isomorphism()
{
    const double real_kk = real_lindhard_kk(
        reference_q, reference_omega, reference_tau, reference_gamma);

    SincQuadrature<> gv_quadrature;
    const double real_gv = gv_quadrature.evaluate_gv_real_lindhard(GvLindhardSincArgs<>{
        reference_q,
        reference_omega,
        reference_tau,
        reference_gamma,
    });

    const double legacy_kk = reference_legacy_kk_real(
        reference_q.raw(),
        reference_omega.raw(),
        reference_tau.raw(),
        reference_gamma.raw());

    assert(std::isfinite(real_kk));
    assert(std::isfinite(real_gv));
    assert(std::isfinite(legacy_kk));

    std::cout << "Re χ^L (KK pipeline) = " << real_kk << '\n';
    std::cout << "Re χ^L (GV direct, Phase 1) = " << real_gv << '\n';
    std::cout << "Re χ^L (legacy KK loop) = " << legacy_kk << '\n';

    // Topological isomorphism: pure `<ranges>` pipeline must reproduce the legacy KK loop.
    const double delta_kk = std::abs(real_kk - legacy_kk);
    assert(delta_kk < 1.0e-6);
    std::cout << "|Δ Re| (KK pipeline vs legacy KK) = " << delta_kk << '\n';

    // Cross-route guardrail at the stable reference point (canonical KK vs grand-canonical GV).
    assert(real_kk < 0.0);
    assert(real_gv < 0.0);
    assert(std::abs(real_kk) < 10.0);
    assert(std::abs(real_gv) < 10.0);

    const double delta_cross = std::abs(real_kk - real_gv);
    std::cout << "|Δ Re| (KK vs GV direct) = " << delta_cross << " (distinct legacy routes)\n";
}

void test_singularity_excision_extreme_cold()
{
    // Below 100 K (and down to 0.01 K-class τ) the reverse-Dedekind excision must stay finite.
    const WaveVector<> q{1.0};
    const Frequency<> omega{0.5};
    const ReducedChemicalPotential<> gamma{1.0};

    for (double tau_raw : {1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-12, 1.0e-20, 0.0}) {
        const LindhardResult<> chi = evaluate_lindhard(
            q, omega, ReducedTemperature<>{tau_raw}, gamma);
        assert(std::isfinite(chi.real()));
        assert(std::isfinite(chi.imag()));
        assert(!std::isnan(chi.real()));
        assert(!std::isnan(chi.imag()));
    }

    const LindhardResult<> chi_t0 =
        evaluate_lindhard(q, omega, ReducedTemperature<>{0.0}, gamma);
    const LindhardResult<> chi_cold =
        evaluate_lindhard(q, omega, ReducedTemperature<>{1.0e-12}, gamma);
    assert(std::abs(chi_cold.real() - chi_t0.real()) < 1.0e-8);
    assert(std::abs(chi_cold.imag() - chi_t0.imag()) < 1.0e-8);

    std::cout << "Singularity excision extreme-cold Re = " << chi_cold.real()
              << ", Im = " << chi_cold.imag() << '\n';
}

}  // namespace

int main()
{
    test_imaginary_lindhard_finite();
    test_evaluate_lindhard_bundle();
    test_kk_vs_gv_real_isomorphism();
    test_singularity_excision_extreme_cold();

    std::cout << "All Lindhard tests passed.\n";
    return 0;
}
