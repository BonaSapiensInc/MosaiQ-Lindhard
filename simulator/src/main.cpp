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
#include "engine/ZetaRPA.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace mosaiq;

inline constexpr const char* output_directory = "output";
inline constexpr const char* structure_factor_filename = "output_structure_factor.dat";
inline constexpr const char* lindhard_base_filename = "output_lindhard_base.dat";
inline constexpr const char* dispersion_filename = "output_plasmon_dispersion.dat";
inline constexpr const char* static_sq_gamma_filename = "output_Sq_gamma.dat";
inline constexpr const char* zeta_scalar_filename = "output_zeta_rpa_scalar.dat";
inline constexpr const char* zeta_dispersion_filename = "output_zeta_rpa_dispersion.dat";

[[nodiscard]] const char* pathway_label(ResponsePathway pathway) noexcept
{
    switch (pathway) {
    case ResponsePathway::StandardRPA:
        return "StandardRPA";
    case ResponsePathway::ZetaRPA:
        return "ZetaRPA";
    case ResponsePathway::ZetaRPA_Experimental:
        return "ZetaRPA_Experimental";
    }
    return "Unknown";
}

[[nodiscard]] std::optional<ResponsePathway> parse_pathway_token(std::string_view token) noexcept
{
    if (token == "standard-rpa" || token == "rpa") {
        return ResponsePathway::StandardRPA;
    }
    if (token == "zeta-rpa") {
        return ResponsePathway::ZetaRPA;
    }
    if (token == "zeta-rpa-experimental") {
        return ResponsePathway::ZetaRPA_Experimental;
    }
    return std::nullopt;
}

void print_usage(const char* program)
{
    std::cerr << "Usage:\n"
              << "  " << program << " [--pathway NAME] [--gamma GAMMA] [--scalar-diagnostic] "
                 "<r_s> <T_kelvin>\n"
              << "  " << program << " [--pathway NAME] --gamma-sweep <r_s> <gamma_1>,<gamma_2>,...\n"
              << "\n"
              << "  r_s       Wigner-Seitz radius in Bohr radii (e.g. 1.0)\n"
              << "  T_kelvin  temperature in kelvin (e.g. 10000.0)\n"
              << "  --pathway standard-rpa | zeta-rpa | zeta-rpa-experimental\n"
              << "            (default: zeta-rpa — multi-component manuscript pipeline)\n"
              << "  --gamma   plasma coupling Gamma for optional scalar diagnostic export\n"
              << "  --scalar-diagnostic  write scalar Zeta-RPA grids (opt-in; does not replace "
                 "S(q,ω))\n"
              << "  gamma     plasma coupling parameter Gamma (Eq. plasma-coupling-parameter-Hartree)\n";
}

void print_pathway_header(ResponsePathway pathway, bool strong_coupling) noexcept
{
    std::cerr << "ResponsePathway: " << pathway_label(pathway);
    if (pathway == ResponsePathway::StandardRPA) {
        std::cerr << " (legacy undressed two-component RPA)\n";
    } else if (pathway == ResponsePathway::ZetaRPA_Experimental) {
        std::cerr << " (experimental multi-component dress)\n";
    } else if (strong_coupling) {
        std::cerr << " (production multi-component Zeta-RPA; strong-coupling window)\n";
    } else {
        std::cerr << " (production multi-component Zeta-RPA; W_ζ → 1 at weak Γ)\n";
    }
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

int run_standard_mode(double rs, double T_kelvin, ResponsePathway pathway)
{
    PlasmaContext plasma{};
    if (!build_plasma_context(rs, T_kelvin, plasma)) {
        std::cerr << "Error: failed to initialize plasma state.\n";
        return EXIT_FAILURE;
    }

    const double gamma_plasma = plasma_coupling_from_kelvin(rs, T_kelvin);
    const CouplingRegime<> regime{
        .rs = rs,
        .gamma_plasma = gamma_plasma,
        .tau = plasma.electron.tau.raw(),
    };
    print_pathway_header(pathway, is_strong_coupling(regime));

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

    structure_factor_output << "# MosaiQ-Lindhard CLI — multi-component structure factors\n";
    structure_factor_output << "# ResponsePathway = " << pathway_label(pathway) << '\n';
    structure_factor_output << "# Gamma_plasma = " << std::scientific << std::setprecision(12)
                             << gamma_plasma << '\n';
    write_run_header(structure_factor_output, rs, T_kelvin, plasma);
    write_run_header(lindhard_base_output, rs, T_kelvin, plasma);
    structure_factor_output << "# omega column: electron reduced units [E_F/e/hbar]\n";
    structure_factor_output << "# S_ei uses tau_e (electron thermal reference)\n";
    structure_factor_output << "# DOS-restored susceptibilities; core columns unchanged for plots\n";
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
            const ZetaRpaMatrixResult<> rpa =
                evaluate_rpa_response(q, omega_e, plasma, pathway);
            const double dos_e = density_of_states_at_fermi(plasma.n, plasma.electron.E_F);
            const double dos_i = density_of_states_at_fermi(plasma.n, plasma.ion.E_F);
            const double dos_cross = std::sqrt(dos_e * dos_i);
            const double S_ee = dynamic_structure_factor(
                rpa.chi_ee / dos_e, omega_e, plasma.electron.tau.raw());
            const double S_ii = dynamic_structure_factor(
                rpa.chi_ii / dos_i, omega_i, plasma.ion.tau.raw());
            const double S_ei = dynamic_structure_factor(
                rpa.chi_ei / dos_cross, omega_e, plasma.electron.tau.raw());

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
    std::cout << "  ResponsePathway = " << pathway_label(pathway) << '\n';
    std::cout << "  r_s = " << rs << ", T = " << T_kelvin << " K, Gamma = " << gamma_plasma
              << '\n';
    std::cout << "  tau_e = " << plasma.electron.tau.raw()
              << ", tau_i = " << plasma.ion.tau.raw() << '\n';
    return EXIT_SUCCESS;
}

int run_gamma_sweep_mode(double rs,
                         const std::vector<double>& gammas,
                         ResponsePathway pathway)
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
    output << "# ResponsePathway = " << pathway_label(pathway) << '\n';
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
                static_omega_step,
                pathway);

            if (!std::isfinite(static_sq.S_ee) || !std::isfinite(static_sq.S_ii) ||
                !std::isfinite(static_sq.S_ei)) {
                continue;
            }

            output << std::fixed << std::setprecision(6) << q << ' '
                   << std::scientific << std::setprecision(12) << static_sq.S_ee << ' '
                   << static_sq.S_ii << ' ' << static_sq.S_ei << '\n';
            output.flush();
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
    std::cout << "  ResponsePathway = " << pathway_label(pathway) << '\n';
    std::cout << "  r_s = " << rs << '\n';
    return EXIT_SUCCESS;
}

/// Scalar Zeta-RPA diagnostic export (opt-in via --scalar-diagnostic).
[[nodiscard]] int run_zeta_scalar_mode(double rs,
                                       double T_kelvin,
                                       double gamma_plasma,
                                       ResponsePathway requested)
{
    PlasmaContext plasma{};
    if (!build_plasma_context(rs, T_kelvin, plasma)) {
        std::cerr << "Error: failed to initialize plasma state.\n";
        return EXIT_FAILURE;
    }

    const CouplingRegime<> regime{
        .rs = rs,
        .gamma_plasma = gamma_plasma,
        .tau = plasma.electron.tau.raw(),
    };
    const bool strong = is_strong_coupling(regime);
    const bool force = !strong;
    const ResponsePathway selected =
        select_response_pathway(regime, requested, force, /*auto_bypass=*/false);

    print_pathway_header(selected, strong);
    if (selected == ResponsePathway::StandardRPA) {
        std::cerr << "Error: pathway selection resolved to StandardRPA; refusing silent fallback.\n"
                  << "  Hint: raise --gamma above Gamma★=" << constants::zeta_rpa_gamma_star
                  << " or use a Zeta pathway with force in the scalar verifier.\n";
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

    const std::filesystem::path out_path = output_dir / zeta_scalar_filename;
    std::ofstream out(out_path);
    if (!out) {
        std::cerr << "Error: cannot open " << out_path << " for writing.\n";
        return EXIT_FAILURE;
    }

    const std::filesystem::path dispersion_path = output_dir / zeta_dispersion_filename;
    std::ofstream dispersion_out(dispersion_path);
    if (!dispersion_out) {
        std::cerr << "Error: cannot open " << dispersion_path << " for writing.\n";
        return EXIT_FAILURE;
    }

    out << "# MosaiQ-Lindhard CLI — scalar Zeta-RPA diagnostic (manuscript pipelines untouched)\n";
    out << "# ResponsePathway = " << pathway_label(selected) << '\n';
    write_run_header(out, rs, T_kelvin, plasma);
    out << "# Gamma_plasma = " << std::scientific << std::setprecision(12) << gamma_plasma << '\n';
    out << "# columns: q omega ReChiL ImChiL ReChiRPA ImChiRPA ReChiZeta ImChiZeta Wzeta "
           "ReEps ImEps\n";

    PlasmonPolePolicy<> pole_policy{};
    pole_policy.bracket.scan_ceiling_factor = 8.0;

    dispersion_out
        << "# MosaiQ-Lindhard CLI — scalar Zeta-RPA vs RPA plasmon comparison (Phase Z2)\n";
    dispersion_out << "# ResponsePathway = " << pathway_label(selected) << '\n';
    write_run_header(dispersion_out, rs, T_kelvin, plasma);
    dispersion_out << "# Gamma_plasma = " << std::scientific << std::setprecision(12)
                    << gamma_plasma << '\n';
    dispersion_out << "# omega columns: electron reduced units [E_F/hbar]\n";
    dispersion_out << "# ImEps evaluated at each pathway's extracted pole (NaN if no root)\n";
    dispersion_out << "# bohm_gross_estimate: long-wavelength Bohm-Gross anchor (same units)\n";
    dispersion_out << "# columns: q omega_p_RPA omega_p_zeta ImEps_RPA ImEps_zeta "
                      "bohm_gross_estimate\n";

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::size_t rows = 0;
    std::size_t finite_ok = 0;
    std::size_t zeta_roots = 0;
    std::size_t rpa_roots = 0;
    const std::size_t total_q = grid_node_count(default_q_min, default_q_max, default_q_step);

    std::cerr << "Tracing scalar RPA / Zeta plasmon comparison over " << total_q
              << " q points...\n";

    for (double q = default_q_min; q <= default_q_max + 0.5 * default_q_step; q += default_q_step) {
        const double v = coulomb_potentials_rational(q).v_ee;
        const double bohm_gross = PlasmonPoleExtractor::bohm_gross_frequency(
            WaveVector<>{q}, plasma.electron.tau, rs);

        PlasmonPoleZetaInputs<> zeta_pole{
            .q = WaveVector<>{q},
            .tau = plasma.electron.tau,
            .gamma = plasma.electron.gamma,
            .bare_potential = v,
            .wigner_seitz_radius = rs,
            .regime = regime,
            .pathway = selected,
            .force_pathway = force,
        };
        PlasmonPoleZetaInputs<> rpa_pole = zeta_pole;
        rpa_pole.pathway = ResponsePathway::StandardRPA;
        rpa_pole.force_pathway = true;

        const auto zeta_state = PlasmonPoleExtractor::extract(zeta_pole, pole_policy);
        const auto rpa_state = PlasmonPoleExtractor::extract(rpa_pole, pole_policy);

        const double omega_rpa = rpa_state ? rpa_state->omega_raw() : nan;
        const double omega_zeta = zeta_state ? zeta_state->omega_raw() : nan;
        const double im_rpa = rpa_state ? rpa_state->landau_damping : nan;
        const double im_zeta = zeta_state ? zeta_state->landau_damping : nan;

        if (rpa_state) {
            ++rpa_roots;
        }
        if (zeta_state) {
            ++zeta_roots;
        }

        dispersion_out << q << ' ' << omega_rpa << ' ' << omega_zeta << ' ' << im_rpa << ' '
                       << im_zeta << ' ' << bohm_gross << '\n';
    }

    std::cerr << "  Scalar RPA roots: " << rpa_roots << " / " << total_q << '\n';
    std::cerr << "  Zeta roots:       " << zeta_roots << " / " << total_q << '\n';

    for (double q = default_q_min; q <= 4.0 + 0.5 * default_q_step; q += default_q_step) {
        const double v = coulomb_potentials_rational(q).v_ee;
        for (double omega = default_omega_min; omega <= default_omega_max + 0.5 * default_omega_step;
             omega += default_omega_step) {
            const LindhardResult<> chi_l = evaluate_lindhard(
                WaveVector<>{q},
                Frequency<>{omega},
                plasma.electron.tau,
                plasma.electron.gamma);
            const std::complex<double> chi_c = as_complex(chi_l);
            const std::complex<double> chi_rpa = evaluate_scalar_rpa(chi_c, v);

            ZetaRpaInputs<> inputs{
                .q = WaveVector<>{q},
                .omega = Frequency<>{omega},
                .chi_lindhard = chi_l,
                .bare_potential = v,
                .regime = regime,
                .pathway = selected,
                .force_pathway = force,
            };
            const auto zeta = evaluate_zeta_rpa(inputs);
            if (!zeta) {
                continue;
            }

            out << q << ' ' << omega << ' ' << chi_c.real() << ' ' << chi_c.imag() << ' '
                << chi_rpa.real() << ' ' << chi_rpa.imag() << ' ' << zeta->chi.real() << ' '
                << zeta->chi.imag() << ' ' << zeta->zeta_weight << ' ' << zeta->epsilon.real()
                << ' ' << zeta->epsilon.imag() << '\n';
            ++rows;
            if (std::isfinite(chi_rpa.real()) && std::isfinite(chi_rpa.imag()) &&
                std::isfinite(zeta->chi.real()) && std::isfinite(zeta->chi.imag())) {
                ++finite_ok;
            }
        }
    }

    if (rows == 0) {
        std::cerr << "Error: no scalar Zeta-RPA rows written.\n";
        return EXIT_FAILURE;
    }

    if (zeta_roots == 0) {
        std::cerr << "Warning: no scalar Zeta plasmon roots were extracted on the q grid.\n";
    }

    std::cout << "Wrote " << rows << " rows to " << out_path << " (" << finite_ok
              << " fully finite)\n";
    std::cout << "Wrote " << total_q << " dispersion rows to " << dispersion_path << " ("
              << rpa_roots << " RPA / " << zeta_roots << " Zeta roots)\n";
    std::cout << "  ResponsePathway = " << pathway_label(selected) << '\n';
    std::cout << "  r_s = " << rs << "  T = " << T_kelvin << " K  Gamma = " << gamma_plasma
              << '\n';
    return EXIT_SUCCESS;
}

struct CliOptions {
    ResponsePathway pathway{ResponsePathway::ZetaRPA};
    bool pathway_explicit{false};
    bool scalar_diagnostic{false};
    std::optional<double> gamma_plasma{};
    std::vector<std::string_view> positionals{};
};

[[nodiscard]] std::optional<CliOptions> parse_cli(int argc, char* argv[])
{
    CliOptions opts;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--pathway") {
            if (i + 1 >= argc) {
                return std::nullopt;
            }
            const auto pathway = parse_pathway_token(argv[++i]);
            if (!pathway) {
                return std::nullopt;
            }
            opts.pathway = *pathway;
            opts.pathway_explicit = true;
            continue;
        }
        if (arg == "--gamma") {
            if (i + 1 >= argc) {
                return std::nullopt;
            }
            char* end = nullptr;
            const double g = std::strtod(argv[++i], &end);
            if (end == argv[i] || !std::isfinite(g) || g <= 0.0) {
                return std::nullopt;
            }
            opts.gamma_plasma = g;
            continue;
        }
        if (arg == "--scalar-diagnostic") {
            opts.scalar_diagnostic = true;
            continue;
        }
        if (arg == "--gamma-sweep") {
            // Handled as a dedicated mode; stash token for caller.
            opts.positionals.push_back(arg);
            for (++i; i < argc; ++i) {
                opts.positionals.push_back(argv[i]);
            }
            return opts;
        }
        if (!arg.empty() && arg.front() == '-') {
            return std::nullopt;
        }
        opts.positionals.push_back(arg);
    }
    return opts;
}

}  // namespace

int main(int argc, char* argv[])
{
    const auto parsed = parse_cli(argc, argv);
    if (!parsed) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (!parsed->positionals.empty() && parsed->positionals[0] == "--gamma-sweep") {
        if (parsed->positionals.size() != 3) {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
        char* end = nullptr;
        const double rs = std::strtod(std::string(parsed->positionals[1]).c_str(), &end);
        if (end == std::string(parsed->positionals[1]).c_str() || !std::isfinite(rs) || rs <= 0.0) {
            std::cerr << "Error: invalid r_s value.\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
        const auto gammas = parse_gamma_list(std::string(parsed->positionals[2]));
        if (!gammas) {
            std::cerr << "Error: invalid Gamma list. Expected comma-separated positive values.\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
        std::cerr << "ResponsePathway: " << pathway_label(parsed->pathway)
                  << " (gamma-sweep static S(q))\n";
        return run_gamma_sweep_mode(rs, *gammas, parsed->pathway);
    }

    if (parsed->positionals.size() != 2) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    char* end = nullptr;
    const std::string rs_token{parsed->positionals[0]};
    const double rs = std::strtod(rs_token.c_str(), &end);
    if (end == rs_token.c_str() || !std::isfinite(rs) || rs <= 0.0) {
        std::cerr << "Error: invalid r_s value.\n";
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    end = nullptr;
    const std::string t_token{parsed->positionals[1]};
    const double T_kelvin = std::strtod(t_token.c_str(), &end);
    if (end == t_token.c_str() || !std::isfinite(T_kelvin) || T_kelvin <= 0.0) {
        std::cerr << "Error: invalid temperature value.\n";
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (parsed->scalar_diagnostic) {
        if (parsed->pathway == ResponsePathway::StandardRPA) {
            std::cerr << "Error: --scalar-diagnostic requires a Zeta pathway "
                         "(--pathway zeta-rpa[...]).\n";
            return EXIT_FAILURE;
        }
        const double gamma = parsed->gamma_plasma.value_or(constants::zeta_rpa_gamma_star);
        return run_zeta_scalar_mode(rs, T_kelvin, gamma, parsed->pathway);
    }

    return run_standard_mode(rs, T_kelvin, parsed->pathway);
}
