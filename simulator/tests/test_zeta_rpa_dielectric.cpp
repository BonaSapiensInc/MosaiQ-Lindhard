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
#include "engine/PlasmonPoleExtractor.hpp"
#include "engine/RPA.hpp"
#include "engine/StructureFactor.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

namespace {

using namespace mosaiq;

[[nodiscard]] bool build_scalar_state(double rs,
                                      double T_kelvin,
                                      double gamma_plasma,
                                      double q,
                                      PlasmonPoleZetaInputs<>& zeta_out,
                                      PlasmonPoleZetaInputs<>& rpa_out) noexcept
{
    PlasmaContext plasma{};
    if (!build_plasma_context(rs, T_kelvin, plasma)) {
        return false;
    }

    const CouplingRegime<> regime{
        .rs = rs,
        .gamma_plasma = gamma_plasma,
        .tau = plasma.electron.tau.raw(),
    };
    const double v = coulomb_potentials_rational(q).v_ee;

    zeta_out = PlasmonPoleZetaInputs<>{
        .q = WaveVector<>{q},
        .tau = plasma.electron.tau,
        .gamma = plasma.electron.gamma,
        .bare_potential = v,
        .wigner_seitz_radius = rs,
        .regime = regime,
        .pathway = ResponsePathway::ZetaRPA,
        .force_pathway = true,
    };

    rpa_out = zeta_out;
    rpa_out.pathway = ResponsePathway::StandardRPA;
    rpa_out.force_pathway = true;
    return true;
}

[[nodiscard]] PlasmonPolePolicy<> strong_coupling_policy() noexcept
{
    PlasmonPolePolicy<> policy{};
    // W_ζ elevates ω_p above the Bohm–Gross anchor; widen the deterministic scan.
    policy.bracket.scan_ceiling_factor = 8.0;
    return policy;
}

void test_zeta_dielectric_root_gamma_100()
{
    constexpr double rs = 2.0;
    constexpr double T_kelvin = 10000.0;
    constexpr double gamma_plasma = 100.0;
    constexpr double q = 0.35;

    PlasmonPoleZetaInputs<> zeta_in{};
    PlasmonPoleZetaInputs<> rpa_in{};
    assert(build_scalar_state(rs, T_kelvin, gamma_plasma, q, zeta_in, rpa_in));
    assert(is_strong_coupling(zeta_in.regime));

    const auto policy = strong_coupling_policy();
    const auto zeta_state = PlasmonPoleExtractor::extract(zeta_in, policy);
    assert(zeta_state.has_value());

    const auto eps_zeta =
        PlasmonPoleExtractor::evaluate_epsilon_zeta(zeta_in, zeta_state->omega_p);
    assert(std::isfinite(eps_zeta.real()));
    assert(std::isfinite(eps_zeta.imag()));
    assert(std::abs(eps_zeta.real()) < 1.0e-8);

    const auto rpa_state = PlasmonPoleExtractor::extract(rpa_in, policy);
    // Scalar RPA may still root; Zeta must shift away from W=1 at Γ=100.
    if (rpa_state) {
        assert(std::abs(zeta_state->omega_raw() - rpa_state->omega_raw()) > 1.0e-6);
        std::cout << "  Scalar RPA root also found (shifted by Zeta): omega_RPA="
                  << rpa_state->omega_raw() << '\n';
    } else {
        std::cout << "  Scalar RPA failed to bracket/root (Zeta rescued).\n";
    }

    std::cout << "Gamma=100 dielectric root OK: q=" << q
              << "  omega_zeta=" << zeta_state->omega_raw()
              << "  Im(eps^zeta)=" << zeta_state->landau_damping
              << "  Re(eps^zeta)=" << eps_zeta.real() << '\n';
}

void test_amplified_v_rpa_collapse_zeta_rescue()
{
    constexpr double rs = 4.0;
    constexpr double T_kelvin = 1000.0;
    constexpr double gamma_plasma = 200.0;
    constexpr double q = 0.5;
    constexpr double v_scale = 50.0;

    PlasmonPoleZetaInputs<> zeta_in{};
    PlasmonPoleZetaInputs<> rpa_in{};
    assert(build_scalar_state(rs, T_kelvin, gamma_plasma, q, zeta_in, rpa_in));
    zeta_in.bare_potential *= v_scale;
    rpa_in.bare_potential *= v_scale;

    const auto policy = strong_coupling_policy();
    const auto zeta_state = PlasmonPoleExtractor::extract(zeta_in, policy);
    assert(zeta_state.has_value());

    const auto eps_zeta =
        PlasmonPoleExtractor::evaluate_epsilon_zeta(zeta_in, zeta_state->omega_p);
    assert(std::abs(eps_zeta.real()) < 1.0e-8);

    // Probe RPA pathology at the Zeta pole frequency (collapse / extreme response).
    const LindhardResult<> chi_l = evaluate_lindhard(
        zeta_in.q, zeta_state->omega_p, zeta_in.tau, zeta_in.gamma);
    const std::complex<double> chi_c = as_complex(chi_l);
    const std::complex<double> eps_rpa{1.0, 0.0};
    const auto eps_r = eps_rpa - rpa_in.bare_potential * chi_c;
    const auto chi_rpa = evaluate_scalar_rpa(chi_c, rpa_in.bare_potential);

    const bool rpa_pathological =
        !std::isfinite(chi_rpa.real()) || !std::isfinite(chi_rpa.imag()) ||
        !std::isfinite(eps_r.real()) || !std::isfinite(eps_r.imag()) ||
        std::abs(chi_rpa) > 1.0e6 || std::abs(eps_r) < 1.0e-12;

    const auto rpa_state = PlasmonPoleExtractor::extract(rpa_in, policy);
    const bool rpa_failed_root = !rpa_state.has_value();

    assert(rpa_pathological || rpa_failed_root ||
           (rpa_state && std::abs(rpa_state->omega_raw() - zeta_state->omega_raw()) > 1.0e-4));

    std::cout << "Amplified-v rescue OK: Gamma=" << gamma_plasma << "  v_scale=" << v_scale
              << "  omega_zeta=" << zeta_state->omega_raw()
              << "  rpa_pathological=" << (rpa_pathological ? "yes" : "no")
              << "  rpa_root=" << (rpa_failed_root ? "none" : "found") << '\n';
}

void test_standard_two_component_unchanged_api()
{
    // Smoke: manuscript PlasmonPoleInputs overload still compiles and roots.
    PlasmaContext plasma{};
    assert(build_plasma_context(2.0, 10000.0, plasma));

    PlasmonPoleInputs<> inputs{
        .q = WaveVector<>{0.5},
        .electron =
            PlasmonSpeciesParams<>{
                .tau = plasma.electron.tau,
                .gamma = plasma.electron.gamma,
                .fermi_energy = plasma.electron.E_F,
                .mass = plasma.electron.mass,
            },
        .ion =
            PlasmonSpeciesParams<>{
                .tau = plasma.ion.tau,
                .gamma = plasma.ion.gamma,
                .fermi_energy = plasma.ion.E_F,
                .mass = plasma.ion.mass,
            },
        .potentials = coulomb_potentials_rational(0.5),
        .wigner_seitz_radius = 2.0,
    };

    const auto state = PlasmonPoleExtractor::extract(inputs);
    assert(state.has_value());
    const double re =
        PlasmonPoleExtractor::evaluate_epsilon(inputs, state->omega_p).real();
    assert(std::abs(re) < 1.0e-8);
    std::cout << "Two-component StandardRPA extract still green: omega_p=" << state->omega_raw()
              << '\n';
}

}  // namespace

int main()
{
    test_zeta_dielectric_root_gamma_100();
    test_amplified_v_rpa_collapse_zeta_rescue();
    test_standard_two_component_unchanged_api();
    std::cout << "Zeta-RPA dielectric (Phase Z2) tests passed.\n";
    return 0;
}
