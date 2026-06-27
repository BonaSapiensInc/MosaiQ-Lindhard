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

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

namespace {

using namespace mosaiq;

void assert_finite_complex(const char* label, const std::complex<double>& value)
{
    assert(std::isfinite(value.real()));
    assert(std::isfinite(value.imag()));
    assert(!std::isnan(value.real()));
    assert(!std::isnan(value.imag()));
    assert(!std::isinf(value.real()));
    assert(!std::isinf(value.imag()));
    std::cout << label << " = " << value << '\n';
}

void test_decoupling_limit()
{
    const std::complex<double> chi_e{0.42, 0.08};
    const std::complex<double> chi_i{0.31, -0.17};

    const BarePotentials<> potentials{
        .v_ee = 1.20,
        .v_ii = 0.85,
        .v_ei = 0.0,
    };

    const std::complex<double> epsilon = evaluate_dielectric(chi_e, chi_i, potentials);
    const std::complex<double> factored =
        (std::complex<double>{1.0, 0.0} - chi_e * potentials.v_ee) *
        (std::complex<double>{1.0, 0.0} - chi_i * potentials.v_ii);

    assert_finite_complex("epsilon (decoupled)", epsilon);
    assert_finite_complex("factored epsilon", factored);
    assert(std::abs(epsilon - factored) < 1.0e-14);

    const RpaResult<> rpa = evaluate_rpa_susceptibility(chi_e, chi_i, potentials);
    assert_finite_complex("chi_ee^RPA", rpa.chi_ee);
    assert_finite_complex("chi_ii^RPA", rpa.chi_ii);
    assert(rpa.chi_ei == (std::complex<double>{0.0, 0.0}));

    const std::complex<double> chi_ee_expected = chi_e / (std::complex<double>{1.0, 0.0} - chi_e * potentials.v_ee);
    const std::complex<double> chi_ii_expected = chi_i / (std::complex<double>{1.0, 0.0} - chi_i * potentials.v_ii);
    assert(std::abs(rpa.chi_ee - chi_ee_expected) < 1.0e-14);
    assert(std::abs(rpa.chi_ii - chi_ii_expected) < 1.0e-14);

    std::cout << "Decoupling limit verified: epsilon factorizes, chi_ei^RPA = 0.\n";
}

void test_coupled_stability()
{
    const LindhardResult<> chi_e = evaluate_lindhard(
        WaveVector<>{1.0},
        Frequency<>{0.5},
        ReducedTemperature<>{0.05},
        ReducedChemicalPotential<>{1.0});

    const LindhardResult<> chi_i = evaluate_lindhard(
        WaveVector<>{1.0},
        Frequency<>{0.5},
        ReducedTemperature<>{0.02},
        ReducedChemicalPotential<>{0.8});

    const BarePotentials<> potentials{
        .v_ee = 1.0,
        .v_ii = 0.6,
        .v_ei = 0.25,
    };

    const std::complex<double> epsilon = evaluate_dielectric(chi_e, chi_i, potentials);
    const RpaResult<> rpa = evaluate_rpa_susceptibility(chi_e, chi_i, potentials);

    assert_finite_complex("epsilon (coupled)", epsilon);
    assert_finite_complex("chi_ee^RPA (coupled)", rpa.chi_ee);
    assert_finite_complex("chi_ii^RPA (coupled)", rpa.chi_ii);
    assert_finite_complex("chi_ei^RPA (coupled)", rpa.chi_ei);

    const std::complex<double> chi_e_c = as_complex(chi_e);
    const std::complex<double> chi_i_c = as_complex(chi_i);

    assert(std::abs(rpa.chi_ee - (chi_e_c - chi_e_c * potentials.v_ii * chi_i_c) / epsilon) < 1.0e-12);
    assert(std::abs(rpa.chi_ii - (chi_i_c - chi_i_c * potentials.v_ee * chi_e_c) / epsilon) < 1.0e-12);
    assert(std::abs(rpa.chi_ei - (chi_e_c * potentials.v_ei * chi_i_c) / epsilon) < 1.0e-12);

    std::cout << "Coupled RPA response remains in stable complex domain.\n";
}

}  // namespace

int main()
{
    test_decoupling_limit();
    test_coupled_stability();

    std::cout << "All RPA tests passed.\n";
    return 0;
}
