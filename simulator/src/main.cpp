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
#include "engine/PlasmonPoleExtractor.hpp"
#include "engine/RPA.hpp"
#include "engine/StructureFactor.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <string_view>

namespace {

using namespace mosaiq;

inline constexpr const char* output_directory = "output";
inline constexpr const char* structure_factor_filename = "output_structure_factor.dat";
inline constexpr const char* lindhard_base_filename = "output_lindhard_base.dat";
inline constexpr const char* dispersion_filename = "output_plasmon_dispersion.dat";
inline constexpr const char* static_sq_gamma_filename = "output_Sq_gamma.dat";

void print_usage(const char* program)
{
    std::cerr << "Usage:\n"
              << "  " << program << " <r_s> <T_kelvin>\n"
              << "  " << program << " --gamma-sweep <r_s> <gamma_1>,<gamma_2>,...\n"
              << "\n"
              << "  r_s       Wigner-Seitz radius in Bohr radii (e.g. 1.0)\n"
              << "  T_kelvin  temperature in kelvin (e.g. 10000.0)\n"
              << "  gamma     plasma coupling parameter Gamma (Eq. plasma-coupling-parameter-Hartree)\n";
}

[[nodiscard]] std::size_t grid_node_count(double min_val, double max_val, double step) noexcept
{
    return static_cast<std::size_t>(std::floor((max_val - min_val) / step)) + 1;
}

[[nodiscard]] PlasmonPoleInputs<> make_plasmon_inputs(double rs,
                                                      double q,
                                                      const PlasmaContext& plasma) noexcept
{
    return PlasmonPoleInputs<>{
        .q = WaveVector<>{q},
        .electron = PlasmonSpeciesParams<>{
            .tau = plasma.electron.tau,
            .gamma = plasma.electron.gamma,
            .fermi_energy = plasma.electron.E_F,
            .mass = plasma.electron.mass,
        },
        .ion = PlasmonSpeciesParams<>{
            .tau = plasma.ion.tau,
            .gamma = plasma.ion.gamma,
            .fermi_energy = plasma.ion.E_F,
            .mass = plasma.ion.mass,
        },
        .potentials = coulomb_potentials_rational(q),
        .wigner_seitz_radius = rs,
    };
}

void write_run_header(std::ostream& out,
                      double rs,
                      double T_kelvin,
                      const PlasmaContext& plasma)
{
    out << std::scientific << std::setprecision(12);
    out << "# r_s = " << rs << "  T_K = " << T_kelvin << '\n';
    out << "# m_e = 1  m_i/m_e = " << ion_mass_ratio << '\n';
    out << "# k_F = " << plasma.electron.k_f << " a_0^-1  n = " << plasma.n << " a_0^-3\n";
    out << "# gamma_e = " << plasma.electron.gamma.raw()
        << "  tau_e = " << plasma.electron.tau.raw()
        << "  gamma_i = " << plasma.ion.gamma.raw()
        << "  tau_i = " << plasma.ion.tau.raw() << '\n';
}

[[nodiscard]] std::size_t sweep_plasmon_dispersion(std::ostream& dispersion_out,
                                                   double rs,
                                                   const PlasmaContext& plasma)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::size_t roots_found = 0;

    for (double q = default_q_min; q <= default_q_max + 0.5 * default_q_step;
         q += default_q_step) {
        const PlasmonPoleInputs<> inputs = make_plasmon_inputs(rs, q, plasma);
        const double bohm_gross = PlasmonPoleExtractor::bohm_gross_frequency(
            WaveVector<>{q},
            plasma.electron.tau,
            rs);

        const auto state = PlasmonPoleExtractor::extract(inputs);
        if (state.has_value()) {
            dispersion_out << q << ' '
                           << state->omega_raw() << ' '
                           << state->landau_damping << ' '
                           << bohm_gross << '\n';
            ++roots_found;
        } else {
            dispersion_out << q << ' ' << nan << ' ' << nan << ' ' << bohm_gross << '\n';
        }
    }

    return roots_found;
}

[[nodiscard]] std::size_t export_bare_lindhard_grid(std::ostream& out,
                                                    const PlasmaContext& plasma)
{
    out << "# MosaiQ-Lindhard CLI — bare Lindhard response manifolds (no RPA screening)\n";
    out << "# omega column: electron reduced units [E_F/e/hbar]\n";
    out << "# chi_i evaluated at ion phase-space (q_i, omega_i); axes use electron (q, omega_e)\n";
    out << "# chi_ei = chi_e * v_ei * chi_i (bare bilinear coupling, Eq. response-ei numerator)\n";
    out << "# columns: q omega_e Re_chi_e Im_chi_e Re_chi_i Im_chi_i Re_chi_ei Im_chi_ei\n";
    out << std::scientific << std::setprecision(12);

    const double q_scale = plasma.electron.k_f / plasma.ion.k_f;
    const double omega_scale = plasma.electron.E_F / plasma.ion.E_F;

    std::size_t rows_written = 0;
    for (double q = default_q_min; q <= dynamic_sq_q_max + 0.5 * default_q_step;
         q += default_q_step) {
        const BarePotentials<> potentials = coulomb_potentials_rational(q);
        const double q_i = q * q_scale;

        for (double omega_e = default_omega_min;
             omega_e <= default_omega_max + 0.5 * default_omega_step;
             omega_e += default_omega_step) {
            const double omega_i = omega_e * omega_scale;

            const LindhardResult<> chi_e = evaluate_lindhard(
                WaveVector<>{q},
                Frequency<>{omega_e},
                plasma.electron.tau,
                plasma.electron.gamma);

            const LindhardResult<> chi_i = evaluate_lindhard(
                WaveVector<>{q_i},
                Frequency<>{omega_i},
                plasma.ion.tau,
                plasma.ion.gamma);

            const std::complex<double> chi_e_c{chi_e.real(), chi_e.imag()};
            const std::complex<double> chi_i_c{chi_i.real(), chi_i.imag()};
            const std::complex<double> chi_ei = chi_e_c * potentials.v_ei * chi_i_c;

            if (!std::isfinite(chi_e.real()) || !std::isfinite(chi_e.imag()) ||
                !std::isfinite(chi_i.real()) || !std::isfinite(chi_i.imag()) ||
                !std::isfinite(chi_ei.real()) || !std::isfinite(chi_ei.imag())) {
                continue;
            }

            out << q << ' ' << omega_e << ' '
                << chi_e.real() << ' ' << chi_e.imag() << ' '
                << chi_i.real() << ' ' << chi_i.imag() << ' '
                << chi_ei.real() << ' ' << chi_ei.imag() << '\n';
            ++rows_written;
        }
    }

    return rows_written;
}

int run_standard_mode(double rs, double T_kelvin)
{
    PlasmaContext plasma{};
    if (!build_plasma_context(rs, T_kelvin, plasma)) {
        std::cerr << "Error: failed to initialize plasma state.\n";
        return EXIT_FAILURE;
    }

    const std::filesystem::path output_dir{output_directory};
    std::error_code ec;
    std::filesystem::create_directories(output_dir, ec);
    if (ec) {
        std::cerr << "Error: cannot create output directory '" << output_dir << "': "
                  << ec.message() << '\n';
        return EXIT_FAILURE;
    }

    const std::filesystem::path structure_factor_path = output_dir / structure_factor_filename;
    const std::filesystem::path lindhard_base_path = output_dir / lindhard_base_filename;
    const std::filesystem::path dispersion_path = output_dir / dispersion_filename;

    std::ofstream structure_factor_output(structure_factor_path);
    if (!structure_factor_output) {
        std::cerr << "Error: cannot open " << structure_factor_path << " for writing.\n";
        return EXIT_FAILURE;
    }

    std::ofstream dispersion_output(dispersion_path);
    if (!dispersion_output) {
        std::cerr << "Error: cannot open " << dispersion_path << " for writing.\n";
        return EXIT_FAILURE;
    }

    std::ofstream lindhard_base_output(lindhard_base_path);
    if (!lindhard_base_output) {
        std::cerr << "Error: cannot open " << lindhard_base_path << " for writing.\n";
        return EXIT_FAILURE;
    }

    structure_factor_output << "# MosaiQ-Lindhard CLI — multi-component RPA structure factors\n";
    write_run_header(structure_factor_output, rs, T_kelvin, plasma);
    write_run_header(lindhard_base_output, rs, T_kelvin, plasma);
    structure_factor_output << "# omega column: electron reduced units [E_F/e/hbar]\n";
    structure_factor_output << "# S_ei uses tau_e (electron thermal reference)\n";
    structure_factor_output << "# columns: q omega Im(chi_ee) S_ee Im(chi_ii) S_ii Im(chi_ei) S_ei\n";

    dispersion_output << "# MosaiQ-Lindhard CLI — plasmon dispersion Re[epsilon]=0 trajectory\n";
    write_run_header(dispersion_output, rs, T_kelvin, plasma);
    dispersion_output << "# omega_p column: electron reduced units [E_F/e/hbar]\n";
    dispersion_output << "# landau_damping = Im[epsilon(q, omega_p)] at the extracted pole\n";
    dispersion_output << "# bohm_gross_estimate: long-wavelength Bohm-Gross anchor (same units)\n";
    dispersion_output << "# missing roots (continuum / no bracket) are written as NaN NaN\n";
    dispersion_output
        << "# columns: q omega_p_electron_units landau_damping bohm_gross_estimate\n";

    const std::size_t total_q =
        grid_node_count(default_q_min, default_q_max, default_q_step);
    const std::size_t dynamic_sq_q_nodes =
        grid_node_count(default_q_min, dynamic_sq_q_max, default_q_step);
    const std::size_t total_omega =
        grid_node_count(default_omega_min, default_omega_max, default_omega_step);
    const std::size_t structure_factor_omega_nodes = grid_node_count(
        default_omega_min, structure_factor_export_omega_max, default_omega_step);

    std::cerr << "Tracing plasmon dispersion over " << total_q << " q points...\n";
    const std::size_t dispersion_roots =
        sweep_plasmon_dispersion(dispersion_output, rs, plasma);
    std::cerr << "  Plasmon roots found: " << dispersion_roots << " / " << total_q << '\n';

    std::cerr << "Exporting bare Lindhard manifolds over " << dynamic_sq_q_nodes
              << " x " << total_omega << " (q, omega) grid...\n";
    const std::size_t lindhard_rows = export_bare_lindhard_grid(lindhard_base_output, plasma);
    if (lindhard_rows == 0) {
        std::cerr << "Error: no finite bare Lindhard grid points were written.\n";
        return EXIT_FAILURE;
    }

    std::cerr << "Scanning dynamic structure factors over " << dynamic_sq_q_nodes
              << " x " << structure_factor_omega_nodes << " (q, omega) grid...\n";

    std::size_t rows_written = 0;
    for (double q = default_q_min; q <= dynamic_sq_q_max + 0.5 * default_q_step;
         q += default_q_step) {
        for (double omega_e = default_omega_min;
             omega_e <= structure_factor_export_omega_max + 0.5 * default_omega_step;
             omega_e += default_omega_step) {
            const double omega_i = omega_e * plasma.electron.E_F / plasma.ion.E_F;
            const RpaResult<> rpa = evaluate_rpa_response(q, omega_e, plasma);

            const double S_ee =
                dynamic_structure_factor(rpa.chi_ee, omega_e, plasma.electron.tau.raw());
            const double S_ii =
                dynamic_structure_factor(rpa.chi_ii, omega_i, plasma.ion.tau.raw());
            const double S_ei =
                dynamic_structure_factor(rpa.chi_ei, omega_e, plasma.electron.tau.raw());

            if (!std::isfinite(S_ee) || !std::isfinite(S_ii) || !std::isfinite(S_ei)) {
                continue;
            }

            structure_factor_output << q << ' ' << omega_e << ' '
                                    << rpa.chi_ee.imag() << ' ' << S_ee << ' '
                                    << rpa.chi_ii.imag() << ' ' << S_ii << ' '
                                    << rpa.chi_ei.imag() << ' ' << S_ei << '\n';
            ++rows_written;
        }
    }

    if (rows_written == 0) {
        std::cerr << "Error: no finite structure-factor grid points were written.\n";
        return EXIT_FAILURE;
    }

    if (dispersion_roots == 0) {
        std::cerr << "Warning: no plasmon roots were extracted on the q grid.\n";
    }

    std::cout << "Wrote " << lindhard_rows << " rows to " << lindhard_base_path << '\n';
    std::cout << "Wrote " << rows_written << " rows to " << structure_factor_path << '\n';
    std::cout << "Wrote " << total_q << " dispersion rows to " << dispersion_path << " ("
              << dispersion_roots << " roots)\n";
    std::cout << "  r_s = " << rs << ", T = " << T_kelvin << " K\n";
    std::cout << "  tau_e = " << plasma.electron.tau.raw()
              << ", tau_i = " << plasma.ion.tau.raw() << '\n';
    return EXIT_SUCCESS;
}

int run_gamma_sweep_mode(double rs, const std::vector<double>& gammas)
{
    const std::filesystem::path output_dir{output_directory};
    std::error_code ec;
    std::filesystem::create_directories(output_dir, ec);
    if (ec) {
        std::cerr << "Error: cannot create output directory '" << output_dir << "': "
                  << ec.message() << '\n';
        return EXIT_FAILURE;
    }

    const std::filesystem::path output_path = output_dir / static_sq_gamma_filename;
    std::ofstream output(output_path);
    if (!output) {
        std::cerr << "Error: cannot open " << output_path << " for writing.\n";
        return EXIT_FAILURE;
    }

    output << "# MosaiQ-Lindhard CLI — static structure factor gamma sweep\n";
    output << "# S(q) = integral S(q, omega) d omega_bar over ["
           << static_omega_min << ", Dynamic Bound] step "
           << static_omega_step << '\n';
    output << std::scientific << std::setprecision(12);

    std::size_t blocks_written = 0;
    for (const double gamma : gammas) {
        const double T_kelvin = temperature_kelvin_from_gamma(rs, gamma);
        PlasmaContext plasma{};
        if (!build_plasma_context(rs, T_kelvin, plasma)) {
            std::cerr << "Error: failed to initialize plasma state for Gamma = " << gamma << ".\n";
            return EXIT_FAILURE;
        }

        output << '\n';
        output << "# Gamma = " << std::fixed << std::setprecision(6) << gamma
               << ", T = " << T_kelvin << " K, r_s = " << rs << '\n';
        output << "# columns: q_bar S_ee(q) S_ii(q) S_ei(q)\n";

        std::cerr << "Integrating static structure factors for Gamma = " << gamma
                  << " (T = " << T_kelvin << " K)...\n";

        std::size_t rows_written = 0;
        for (double q = default_q_min; q <= default_q_max + 0.5 * default_q_step;
             q += default_q_step) {
            std::cout << "\r  [Progress] Calculating q_bar = " << std::fixed << std::setprecision(2)
                      << q << " / " << default_q_max << " ... " << std::flush;

            const StaticStructureFactor static_sq = integrate_static_structure_factor(
                q,
                plasma,
                static_omega_min,
                static_omega_max,
                static_omega_step);

            if (!std::isfinite(static_sq.S_ee) || !std::isfinite(static_sq.S_ii) ||
                !std::isfinite(static_sq.S_ei)) {
                continue;
            }

            output << std::fixed << std::setprecision(6) << q << ' '
                   << std::scientific << std::setprecision(12) << static_sq.S_ee << ' '
                   << static_sq.S_ii << ' ' << static_sq.S_ei << '\n';
            ++rows_written;
        }

        std::cout << '\n';

        if (rows_written == 0) {
            std::cerr << "Error: no finite static structure-factor rows for Gamma = " << gamma
                      << ".\n";
            return EXIT_FAILURE;
        }

        std::cerr << "  Wrote " << rows_written << " q points.\n";
        ++blocks_written;
    }

    std::cout << "Wrote " << blocks_written << " Gamma blocks to " << output_path << '\n';
    std::cout << "  r_s = " << rs << '\n';
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char* argv[])
{
    if (argc == 3) {
        char* end = nullptr;
        const double rs = std::strtod(argv[1], &end);
        if (end == argv[1] || !std::isfinite(rs) || rs <= 0.0) {
            std::cerr << "Error: invalid r_s value.\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        end = nullptr;
        const double T_kelvin = std::strtod(argv[2], &end);
        if (end == argv[2] || !std::isfinite(T_kelvin) || T_kelvin <= 0.0) {
            std::cerr << "Error: invalid temperature value.\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        return run_standard_mode(rs, T_kelvin);
    }

    if (argc == 4 && std::string_view{argv[1]} == "--gamma-sweep") {
        char* end = nullptr;
        const double rs = std::strtod(argv[2], &end);
        if (end == argv[2] || !std::isfinite(rs) || rs <= 0.0) {
            std::cerr << "Error: invalid r_s value.\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        const auto gammas = parse_gamma_list(argv[3]);
        if (!gammas) {
            std::cerr << "Error: invalid Gamma list. Expected comma-separated positive values.\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        return run_gamma_sweep_mode(rs, *gammas);
    }

    print_usage(argv[0]);
    return EXIT_FAILURE;
}
