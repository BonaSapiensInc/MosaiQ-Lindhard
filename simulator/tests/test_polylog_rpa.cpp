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
#include "engine/Lindhard.hpp"
#include "engine/PolyLogRPA.hpp"
#include "engine/RPA.hpp"
#include "engine/ZetaRPA.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

namespace {

using namespace mosaiq;

void test_pathway_keeps_polylog()
{
    const CouplingRegime<> weak{.rs = 1.0, .gamma_plasma = 1.0, .tau = 0.1};
    const CouplingRegime<> strong{.rs = 1.0, .gamma_plasma = 50.0, .tau = 0.1};

    assert(select_response_pathway(weak, ResponsePathway::PolyLogRPA) ==
           ResponsePathway::PolyLogRPA);
    assert(select_response_pathway(strong, ResponsePathway::PolyLogRPA) ==
           ResponsePathway::PolyLogRPA);

    std::cout << "PolyLog pathway selection OK.\n";
}

void test_weak_coupling_recovers_rpa()
{
    const LindhardResult<> chi_l = evaluate_lindhard(
        WaveVector<>{1.0},
        Frequency<>{0.5},
        ReducedTemperature<>{0.05},
        ReducedChemicalPotential<>{1.0});

    const double v = 1.2;
    const std::complex<double> chi_c = as_complex(chi_l);
    const std::complex<double> chi_rpa = evaluate_scalar_rpa(chi_c, v);

    // Γ → 0 ⇒ f → 0 ⇒ s → 0 ⇒ ordinary scalar RPA identity gate.
    PolyLogRpaInputs<> inputs{
        .q = WaveVector<>{1.0},
        .omega = Frequency<>{0.5},
        .chi_lindhard = chi_l,
        .bare_potential = v,
        .regime = {.rs = 1.0, .gamma_plasma = 0.0, .tau = 0.05},
    };

    const auto pl = evaluate_polylog_rpa(inputs);
    assert(pl.has_value());
    assert(pl->pathway == ResponsePathway::PolyLogRPA);
    assert(std::abs(pl->s_parameter) <= 1.0e-15);
    assert(std::abs(pl->chi - chi_rpa) < 1.0e-14);

    const std::complex<double> one{1.0, 0.0};
    assert(std::abs(pl->epsilon - (one - v * chi_c)) < 1.0e-14);

    std::cout << "Weak-coupling PolyLog → RPA identity OK.\n";
}

void test_finite_s_finite_result()
{
    const LindhardResult<> chi_l = evaluate_lindhard(
        WaveVector<>{0.5},
        Frequency<>{1.0},
        ReducedTemperature<>{0.05},
        ReducedChemicalPotential<>{1.0});

    // Moderate Γ ⇒ s = f = αΓ/(1+…) > 0; |v χ^L| typically inside series disk.
    PolyLogRpaInputs<> inputs{
        .q = WaveVector<>{0.5},
        .omega = Frequency<>{1.0},
        .chi_lindhard = chi_l,
        .bare_potential = 0.5,
        .regime = {.rs = 1.0, .gamma_plasma = 20.0, .tau = 0.05},
    };

    const auto f = evaluate_coupling_shape_factor(inputs.regime, inputs.weight_params);
    assert(f.has_value());
    assert(*f > 1.0e-14);

    const auto pl = evaluate_polylog_rpa(inputs);
    assert(pl.has_value());
    assert(std::abs(pl->s_parameter - *f) < 1.0e-15);
    assert(std::isfinite(pl->chi.real()) && std::isfinite(pl->chi.imag()));
    assert(std::isfinite(pl->epsilon.real()) && std::isfinite(pl->epsilon.imag()));

    // Finite-s path uses χ = χ^L Li_s(z); must differ from the s→0 RPA gate in general.
    const std::complex<double> chi_rpa = evaluate_scalar_rpa(as_complex(chi_l), inputs.bare_potential);
    assert(std::abs(pl->chi - chi_rpa) > 1.0e-10);

    std::cout << "Finite-s PolyLog evaluation OK (s=" << pl->s_parameter << ").\n";
}

}  // namespace

int main()
{
    test_pathway_keeps_polylog();
    test_weak_coupling_recovers_rpa();
    test_finite_s_finite_result();
    std::cout << "All PolyLog-RPA tests passed.\n";
    return 0;
}
