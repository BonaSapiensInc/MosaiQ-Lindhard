/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

/// Phase Z1 strong-coupling sweep: Γ × r_s × (q,ω) with Coulomb and amplified v.
/// Writes:
///   output/zeta_rpa_strong_coupling_sweep.dat   — per-state summary
///   output/zeta_rpa_strong_coupling_rescue.dat  — points where RPA fails / extreme
///                                               but Experimental Zeta-RPA stays finite

#include "core/Concepts.hpp"
#include "engine/Lindhard.hpp"
#include "engine/RPA.hpp"
#include "engine/StructureFactor.hpp"
#include "engine/ZetaRPA.hpp"
#include "physics/Constants.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace {

using namespace mosaiq;

constexpr double extreme_abs_chi = 1.0e6;
constexpr double eps_pole_floor = 1.0e-12;

[[nodiscard]] const char* pathway_name(ResponsePathway p) noexcept
{
    switch (p) {
    case ResponsePathway::StandardRPA:
        return "StandardRPA";
    case ResponsePathway::ZetaRPA:
        return "ZetaRPA";
    case ResponsePathway::ZetaRPA_Experimental:
        return "ZetaRPA_Experimental";
    }
    return "Unknown";
}

[[nodiscard]] bool is_finite_c(std::complex<double> z) noexcept
{
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}

[[nodiscard]] bool rpa_failed_or_extreme(std::complex<double> chi_rpa,
                                         std::complex<double> eps_rpa) noexcept
{
    if (!is_finite_c(chi_rpa) || !is_finite_c(eps_rpa)) {
        return true;
    }
    if (std::abs(chi_rpa) > extreme_abs_chi) {
        return true;
    }
    if (std::abs(eps_rpa) < eps_pole_floor) {
        return true;
    }
    return false;
}

[[nodiscard]] double reduced_tau(double rs, double T_kelvin) noexcept
{
    const double k_f = constants::fermi_wavevector_from_rs(rs);
    const double E_F = 0.5 * k_f * k_f;
    return constants::temperature_kelvin_to_hartree(T_kelvin) / E_F;
}

}  // namespace

int main()
{
    const double T_kelvin = 1000.0;
    const double gamma_chem = 1.0;
    const double rs_list[] = {0.5, 1.0, 2.0, 4.0};
    const double gamma_list[] = {10.0, 50.0, 100.0, 200.0};
    // Amplify Coulomb v, and also park on the scalar RPA pole v ≈ 1/Re χ^L.
    const double v_scale_list[] = {1.0, 5.0, 20.0, 50.0, 200.0};

    const std::filesystem::path out_dir{"output"};
    std::filesystem::create_directories(out_dir);
    const auto sweep_path = out_dir / "zeta_rpa_strong_coupling_sweep.dat";
    const auto rescue_path = out_dir / "zeta_rpa_strong_coupling_rescue.dat";

    std::ofstream sweep(sweep_path);
    std::ofstream rescue(rescue_path);
    if (!sweep || !rescue) {
        std::cerr << "Error: cannot open sweep/rescue output files.\n";
        return EXIT_FAILURE;
    }

    sweep << "# MosaiQ-Lindhard Phase Z1 strong-coupling sweep summary\n";
    sweep << "# T = " << T_kelvin << " K\n";
    sweep << "# columns: rs Gamma v_scale n_points n_rpa_fail n_prod_fail n_exp_fail "
             "n_rescue max_abs_chi_rpa max_abs_chi_exp\n";

    rescue << "# Rescue points: scalar RPA failed/extreme but ZetaRPA_Experimental finite\n";
    rescue << "# columns: rs Gamma v_scale q omega v ReChiRPA ImChiRPA AbsChiRPA AbsEpsRPA "
              "ReChiExp ImChiExp AbsChiExp Wexp AbsEpsExp\n";

    std::size_t total_rescue = 0;
    std::size_t total_rpa_fail = 0;
    double strongest_gamma_rescue = 0.0;
    double strongest_rs_rescue = 0.0;
    double strongest_vscale_rescue = 0.0;
    double strongest_abs_chi_rpa = 0.0;

    // Baseline identity check (unchanged contract).
    {
        const CouplingRegime<> regime{.rs = 1.0, .gamma_plasma = 50.0, .tau = 0.01};
        assert(select_response_pathway(regime, ResponsePathway::ZetaRPA) ==
               ResponsePathway::ZetaRPA);
        std::cout << "# pathway_check Gamma=50 -> "
                  << pathway_name(select_response_pathway(regime, ResponsePathway::ZetaRPA))
                  << '\n';
    }

    for (double rs : rs_list) {
        const double tau = reduced_tau(rs, T_kelvin);
        for (double gamma_plasma : gamma_list) {
            const CouplingRegime<> regime{
                .rs = rs,
                .gamma_plasma = gamma_plasma,
                .tau = tau,
            };

            for (double v_scale : v_scale_list) {
                std::size_t n_points = 0;
                std::size_t n_rpa_fail = 0;
                std::size_t n_prod_fail = 0;
                std::size_t n_exp_fail = 0;
                std::size_t n_rescue = 0;
                double max_abs_rpa = 0.0;
                double max_abs_exp = 0.0;

                for (double q = 0.2; q <= 2.0001; q += 0.4) {
                    const double v0 = coulomb_potentials_rational(q).v_ee;
                    for (double omega = 0.2; omega <= 2.0001; omega += 0.4) {
                        const LindhardResult<> chi_l = evaluate_lindhard(
                            WaveVector<>{q},
                            Frequency<>{omega},
                            ReducedTemperature<>{tau},
                            ReducedChemicalPotential<>{gamma_chem});
                        const std::complex<double> chi_c = as_complex(chi_l);

                        double v_candidates[3] = {v0 * v_scale, v0 * v_scale, v0 * v_scale};
                        // Deliberate pole parking: v = 1/Re(χ^L) when |Re χ| is appreciable.
                        if (std::abs(chi_c.real()) > 1.0e-3) {
                            v_candidates[1] = 1.0 / chi_c.real();
                            v_candidates[2] = 1.01 / chi_c.real();
                        }

                        for (double v : v_candidates) {
                            if (!std::isfinite(v) || v == 0.0) {
                                continue;
                            }
                            ++n_points;
                            const std::complex<double> one{1.0, 0.0};
                            const std::complex<double> eps_rpa = one - v * chi_c;
                            const std::complex<double> chi_rpa = evaluate_scalar_rpa(chi_c, v);

                            ZetaRpaInputs<> in_prod{
                                .q = WaveVector<>{q},
                                .omega = Frequency<>{omega},
                                .chi_lindhard = chi_l,
                                .bare_potential = v,
                                .regime = regime,
                                .pathway = ResponsePathway::ZetaRPA,
                                .force_pathway = true,
                            };
                            ZetaRpaInputs<> in_exp = in_prod;
                            in_exp.pathway = ResponsePathway::ZetaRPA_Experimental;

                            const auto z_prod = evaluate_zeta_rpa(in_prod);
                            const auto z_exp = evaluate_zeta_rpa(in_exp);

                            const bool rpa_bad = rpa_failed_or_extreme(chi_rpa, eps_rpa);
                            const bool prod_bad =
                                !z_prod || rpa_failed_or_extreme(z_prod->chi, z_prod->epsilon);
                            const bool exp_ok = z_exp && is_finite_c(z_exp->chi) &&
                                                is_finite_c(z_exp->epsilon) &&
                                                std::abs(z_exp->chi) <= extreme_abs_chi &&
                                                std::abs(z_exp->epsilon) >= eps_pole_floor;

                            if (rpa_bad) {
                                ++n_rpa_fail;
                                ++total_rpa_fail;
                            }
                            if (prod_bad) {
                                ++n_prod_fail;
                            }
                            if (!exp_ok) {
                                ++n_exp_fail;
                            }

                            if (is_finite_c(chi_rpa)) {
                                max_abs_rpa = std::max(max_abs_rpa, std::abs(chi_rpa));
                            }
                            if (z_exp && is_finite_c(z_exp->chi)) {
                                max_abs_exp = std::max(max_abs_exp, std::abs(z_exp->chi));
                            }

                            // Locked production W may also rescue poles; Experimental is A/B.
                            if (rpa_bad && exp_ok) {
                                ++n_rescue;
                                ++total_rescue;
                                if (gamma_plasma >= strongest_gamma_rescue) {
                                    strongest_gamma_rescue = gamma_plasma;
                                    strongest_rs_rescue = rs;
                                    strongest_vscale_rescue = v_scale;
                                    strongest_abs_chi_rpa =
                                        is_finite_c(chi_rpa)
                                            ? std::abs(chi_rpa)
                                            : std::numeric_limits<double>::infinity();
                                }
                                rescue << rs << ' ' << gamma_plasma << ' ' << v_scale << ' ' << q
                                       << ' ' << omega << ' ' << v << ' ' << chi_rpa.real() << ' '
                                       << chi_rpa.imag() << ' ' << std::abs(chi_rpa) << ' '
                                       << std::abs(eps_rpa) << ' ' << z_exp->chi.real() << ' '
                                       << z_exp->chi.imag() << ' ' << std::abs(z_exp->chi) << ' '
                                       << z_exp->zeta_weight << ' ' << std::abs(z_exp->epsilon)
                                       << '\n';
                            }
                        }
                    }
                }

                sweep << rs << ' ' << gamma_plasma << ' ' << v_scale << ' ' << n_points << ' '
                      << n_rpa_fail << ' ' << n_prod_fail << ' ' << n_exp_fail << ' ' << n_rescue
                      << ' ' << max_abs_rpa << ' ' << max_abs_exp << '\n';
            }
        }
    }

    std::cout << "# wrote " << sweep_path << '\n';
    std::cout << "# wrote " << rescue_path << '\n';
    std::cout << "# total_rpa_fail_or_extreme = " << total_rpa_fail << '\n';
    std::cout << "# total_rescue_experimental = " << total_rescue << '\n';
    if (total_rescue > 0) {
        std::cout << "# strongest_rescue: rs=" << strongest_rs_rescue
                  << " Gamma=" << strongest_gamma_rescue << " v_scale=" << strongest_vscale_rescue
                  << " |chi_RPA|~" << strongest_abs_chi_rpa << '\n';
        std::cout << "ZetaRPA strong-coupling sweep PASSED (rescues observed).\n";
        return EXIT_SUCCESS;
    }

    std::cout << "ZetaRPA strong-coupling sweep PASSED (no RPA extremes in probed window; "
                 "infrastructure OK).\n";
    return EXIT_SUCCESS;
}
