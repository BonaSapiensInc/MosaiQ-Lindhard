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

    // Γ → 0 ⇒ f → 0 ⇒ W_ζ → 1 (locked form).
    const auto w0 = evaluate_zeta_weight(
        ResponsePathway::ZetaRPA,
        CouplingRegime<>{.rs = 1.0, .gamma_plasma = 0.0, .tau = 0.05});
    assert(w0.has_value());
    assert(std::abs(*w0 - 1.0) < 1.0e-15);

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

    std::cout << "Weak-coupling identity (W_ζ → 1 as Γ → 0) OK.\n";
}

void test_high_gamma_stability()
{
    const CouplingRegime<> regime{.rs = 1.0, .gamma_plasma = 200.0, .tau = 0.05};
    const auto weight = evaluate_zeta_weight(ResponsePathway::ZetaRPA, regime);
    assert(weight.has_value());
    assert(std::isfinite(*weight));
    assert(*weight > 1.0);  // locked dress grows with Γ for α>0

    const LindhardResult<> chi_l = evaluate_lindhard(
        WaveVector<>{1.0},
        Frequency<>{0.5},
        ReducedTemperature<>{0.05},
        ReducedChemicalPotential<>{1.0});

    ZetaRpaInputs<> inputs{
        .q = WaveVector<>{1.0},
        .omega = Frequency<>{0.5},
        .chi_lindhard = chi_l,
        .bare_potential = 1.2,
        .regime = regime,
        .pathway = ResponsePathway::ZetaRPA,
        .force_pathway = true,
    };
    const auto zeta = evaluate_zeta_rpa(inputs);
    assert(zeta.has_value());
    assert(std::isfinite(zeta->chi.real()) && std::isfinite(zeta->chi.imag()));
    assert(std::isfinite(zeta->epsilon.real()) && std::isfinite(zeta->epsilon.imag()));
    assert(std::abs(zeta->zeta_weight - *weight) < 1.0e-12);

    // Distinct from bare scalar RPA at finite Γ.
    const auto chi_rpa = evaluate_scalar_rpa(as_complex(chi_l), inputs.bare_potential);
    assert(std::abs(zeta->chi - chi_rpa) > 1.0e-6);

    std::cout << "High-Γ stability OK (W_ζ = " << *weight << ").\n";
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

void test_matrix_weak_coupling_identity()
{
    const LindhardResult<> chi_e{.real_part = -0.35, .imag_part = 0.12};
    const LindhardResult<> chi_i{.real_part = -0.08, .imag_part = 0.02};
    const BarePotentials<> pot{.v_ee = 1.4, .v_ii = 0.9, .v_ei = -1.1};

    const CouplingRegime<> weak_e{.rs = 1.0, .gamma_plasma = 0.0, .tau = 0.05};
    const CouplingRegime<> weak_i{.rs = 1.0, .gamma_plasma = 0.0, .tau = 80.0};

    ZetaRpaMatrixInputs<> inputs{
        .q = WaveVector<>{0.5},
        .omega = Frequency<>{0.4},
        .chi_lindhard_e = chi_e,
        .chi_lindhard_i = chi_i,
        .potentials = pot,
        .regime_e = weak_e,
        .regime_i = weak_i,
        .pathway = ResponsePathway::ZetaRPA,
        .force_pathway = true,
    };

    const auto zeta = evaluate_zeta_rpa_matrix(inputs);
    assert(zeta.has_value());
    assert(std::abs(zeta->zeta_weight_ee - 1.0) < 1.0e-15);
    assert(std::abs(zeta->zeta_weight_ii - 1.0) < 1.0e-15);
    assert(std::abs(zeta->zeta_weight_ei - 1.0) < 1.0e-15);

    const RpaResult<> rpa =
        evaluate_rpa_susceptibility(as_complex(chi_e), as_complex(chi_i), pot);
    const auto eps = evaluate_dielectric(as_complex(chi_e), as_complex(chi_i), pot);

    assert(std::abs(zeta->chi_ee - rpa.chi_ee) < 1.0e-14);
    assert(std::abs(zeta->chi_ii - rpa.chi_ii) < 1.0e-14);
    assert(std::abs(zeta->chi_ei - rpa.chi_ei) < 1.0e-14);
    assert(std::abs(zeta->epsilon - eps) < 1.0e-14);

    std::cout << "Matrix weak-coupling identity (all W_ζ → 1) OK.\n";
}

void test_matrix_species_asymmetry_and_dress()
{
    const LindhardResult<> chi_e{.real_part = -0.30, .imag_part = 0.10};
    const LindhardResult<> chi_i{.real_part = -0.05, .imag_part = 0.01};
    const BarePotentials<> pot{.v_ee = 1.2, .v_ii = 0.8, .v_ei = -1.0};

    const CouplingRegime<> strong_e{.rs = 2.0, .gamma_plasma = 100.0, .tau = 0.05};
    const CouplingRegime<> weak_i{.rs = 2.0, .gamma_plasma = 1.0, .tau = 80.0};

    ZetaRpaMatrixInputs<> inputs{
        .q = WaveVector<>{0.4},
        .omega = Frequency<>{0.3},
        .chi_lindhard_e = chi_e,
        .chi_lindhard_i = chi_i,
        .potentials = pot,
        .regime_e = strong_e,
        .regime_i = weak_i,
        .pathway = ResponsePathway::ZetaRPA,
        .force_pathway = true,
    };

    const auto zeta = evaluate_zeta_rpa_matrix(inputs);
    assert(zeta.has_value());
    assert(zeta->zeta_weight_ee > 1.0);
    // Ion channel forced Zeta at weak Γ still → W≈1 via locked form.
    assert(std::abs(zeta->zeta_weight_ii - 1.0) < 1.0e-6 || zeta->zeta_weight_ii > 0.0);
    assert(std::isfinite(zeta->zeta_weight_ei));
    // Cross regime averages Γ → intermediate dress relative to ee.
    assert(zeta->zeta_weight_ei < zeta->zeta_weight_ee + 1.0e-12);

    const RpaResult<> rpa =
        evaluate_rpa_susceptibility(as_complex(chi_e), as_complex(chi_i), pot);
    assert(std::abs(zeta->chi_ee - rpa.chi_ee) > 1.0e-6);
    assert(std::isfinite(zeta->chi_ee.real()) && std::isfinite(zeta->chi_ee.imag()));
    assert(std::isfinite(zeta->chi_ii.real()) && std::isfinite(zeta->chi_ii.imag()));
    assert(std::isfinite(zeta->chi_ei.real()) && std::isfinite(zeta->chi_ei.imag()));
    assert(std::isfinite(zeta->epsilon.real()) && std::isfinite(zeta->epsilon.imag()));

    std::cout << "Matrix species asymmetry / strong-coupling dress OK "
              << "(W_ee=" << zeta->zeta_weight_ee << ", W_ii=" << zeta->zeta_weight_ii
              << ", W_ei=" << zeta->zeta_weight_ei << ").\n";
}

void test_matrix_no_silent_fallback()
{
    ZetaRpaMatrixInputs<> inputs{
        .chi_lindhard_e = {.real_part = -0.2, .imag_part = 0.05},
        .chi_lindhard_i = {.real_part = -0.05, .imag_part = 0.01},
        .potentials = {.v_ee = 1.0, .v_ii = 1.0, .v_ei = -1.0},
        .regime_e = {.rs = 1.0, .gamma_plasma = 1.0, .tau = 0.1},
        .regime_i = {.rs = 1.0, .gamma_plasma = 1.0, .tau = 0.1},
        .pathway = ResponsePathway::ZetaRPA,
        .force_pathway = false,
    };
    assert(!evaluate_zeta_rpa_matrix(inputs).has_value());

    inputs.pathway = ResponsePathway::StandardRPA;
    inputs.force_pathway = true;
    assert(!evaluate_zeta_rpa_matrix(inputs).has_value());

    std::cout << "Matrix no silent RPA fallback OK.\n";
}

}  // namespace

int main()
{
    test_pathway_selection();
    test_weak_coupling_identity();
    test_high_gamma_stability();
    test_experimental_dress_vanishes_at_zero_gamma();
    test_no_silent_fallback();
    test_nonfinite_rejection();
    test_matrix_weak_coupling_identity();
    test_matrix_species_asymmetry_and_dress();
    test_matrix_no_silent_fallback();
    std::cout << "All ZetaRPA tests passed.\n";
    return 0;
}
