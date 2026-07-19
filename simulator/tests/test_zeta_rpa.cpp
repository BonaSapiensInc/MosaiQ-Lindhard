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
#include "engine/RPA.hpp"
#include "engine/ZetaRPA.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>

namespace {

using namespace mosaiq;

void test_pathway_selection()
{
    const CouplingRegime<> weak{.rs = 1.0, .gamma_plasma = 1.0, .tau = 0.1};
    const CouplingRegime<> strong{.rs = 1.0, .gamma_plasma = 50.0, .tau = 0.1};

    assert(!is_strong_coupling(weak));
    assert(is_strong_coupling(strong));

    assert(select_response_pathway(weak, ResponsePathway::StandardRPA) ==
           ResponsePathway::StandardRPA);
    assert(select_response_pathway(weak, ResponsePathway::ZetaRPA) ==
           ResponsePathway::StandardRPA);
    assert(select_response_pathway(weak, ResponsePathway::ZetaRPA, /*force=*/true) ==
           ResponsePathway::ZetaRPA);
    assert(select_response_pathway(strong, ResponsePathway::ZetaRPA) == ResponsePathway::ZetaRPA);
    assert(select_response_pathway(strong, ResponsePathway::StandardRPA, false, /*auto=*/true) ==
           ResponsePathway::ZetaRPA);
    assert(select_response_pathway(strong, ResponsePathway::ZetaRPA_Experimental) ==
           ResponsePathway::ZetaRPA_Experimental);
    // Auto-bypass never selects Experimental.
    assert(select_response_pathway(strong, ResponsePathway::StandardRPA, false, true) !=
           ResponsePathway::ZetaRPA_Experimental);

    std::cout << "Pathway selection OK.\n";
}

void test_weak_coupling_identity()
{
    const LindhardResult<> chi_l = evaluate_lindhard(
        WaveVector<>{1.0},
        Frequency<>{0.5},
        ReducedTemperature<>{0.05},
        ReducedChemicalPotential<>{1.0});

    const double v = 1.2;
    const std::complex<double> chi_c = as_complex(chi_l);
    const std::complex<double> chi_rpa = evaluate_scalar_rpa(chi_c, v);

    ZetaRpaInputs<> inputs{
        .q = WaveVector<>{1.0},
        .omega = Frequency<>{0.5},
        .chi_lindhard = chi_l,
        .bare_potential = v,
        .regime = {.rs = 1.0, .gamma_plasma = 0.0, .tau = 0.05},
        .pathway = ResponsePathway::ZetaRPA,
        .force_pathway = true,
    };

    const auto zeta = evaluate_zeta_rpa(inputs);
    assert(zeta.has_value());
    assert(zeta->pathway == ResponsePathway::ZetaRPA);
    assert(std::abs(zeta->zeta_weight - 1.0) < 1.0e-15);
    assert(std::abs(zeta->chi - chi_rpa) < 1.0e-14);

    // Also match decoupled two-component RPA ee channel.
    const BarePotentials<> pot{.v_ee = v, .v_ii = 0.0, .v_ei = 0.0};
    const RpaResult<> rpa = evaluate_rpa_susceptibility(chi_c, std::complex<double>{0.0, 0.0}, pot);
    assert(std::abs(zeta->chi - rpa.chi_ee) < 1.0e-14);

    std::cout << "Weak-coupling identity (W_ζ → 1) OK.\n";
}

void test_experimental_dress_vanishes_at_zero_gamma()
{
    const LindhardResult<> chi_l{.real_part = 0.25, .imag_part = 0.05};
    ZetaRpaInputs<> inputs{
        .q = WaveVector<>{0.8},
        .omega = Frequency<>{0.3},
        .chi_lindhard = chi_l,
        .bare_potential = 0.75,
        .regime = {.rs = 1.0, .gamma_plasma = 0.0, .tau = 0.1},
        .pathway = ResponsePathway::ZetaRPA_Experimental,
        .force_pathway = true,
    };

    const auto weight = evaluate_zeta_weight(ResponsePathway::ZetaRPA_Experimental, inputs.regime);
    assert(weight.has_value());
    assert(std::abs(*weight - 1.0) < 1.0e-15);

    const auto zeta = evaluate_zeta_rpa(inputs);
    assert(zeta.has_value());
    const auto chi_rpa = evaluate_scalar_rpa(as_complex(chi_l), inputs.bare_potential);
    assert(std::abs(zeta->chi - chi_rpa) < 1.0e-14);
    std::cout << "Experimental W_ζ(Γ=0) → 1 OK.\n";
}

void test_no_silent_fallback()
{
    ZetaRpaInputs<> inputs{
        .chi_lindhard = {.real_part = 0.1, .imag_part = 0.0},
        .bare_potential = 1.0,
        .regime = {.rs = 1.0, .gamma_plasma = 1.0, .tau = 0.1},
        .pathway = ResponsePathway::ZetaRPA,
        .force_pathway = false,  // weak regime → selector returns StandardRPA
    };
    assert(!evaluate_zeta_rpa(inputs).has_value());

    inputs.pathway = ResponsePathway::StandardRPA;
    inputs.force_pathway = true;
    assert(!evaluate_zeta_rpa(inputs).has_value());

    std::cout << "No silent RPA fallback OK.\n";
}

void test_nonfinite_rejection()
{
    ZetaRpaInputs<> inputs{
        .chi_lindhard = {.real_part = 0.1, .imag_part = 0.0},
        .bare_potential = std::numeric_limits<double>::quiet_NaN(),
        .regime = {.rs = 1.0, .gamma_plasma = 50.0, .tau = 0.1},
        .pathway = ResponsePathway::ZetaRPA,
    };
    assert(!evaluate_zeta_rpa(inputs).has_value());
    std::cout << "Non-finite rejection OK.\n";
}

}  // namespace

int main()
{
    test_pathway_selection();
    test_weak_coupling_identity();
    test_experimental_dress_vanishes_at_zero_gamma();
    test_no_silent_fallback();
    test_nonfinite_rejection();
    std::cout << "All ZetaRPA tests passed.\n";
    return 0;
}
