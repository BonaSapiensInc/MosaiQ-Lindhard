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
#include "physics/Constants.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

using namespace mosaiq;

inline constexpr double q_bar = 1.0;
inline constexpr double omega_bar = 0.5;
inline constexpr double gamma_value = 1.0;

inline constexpr const char* output_directory = "output";
inline constexpr const char* output_filename = "output_lindhard_rs_t_sweep.dat";

[[nodiscard]] double fermi_temperature_kelvin(double rs) noexcept
{
    const double k_f = constants::fermi_wavevector_from_rs(rs);
    const double e_f = k_f * k_f / (2.0 * constants::electron_mass_amu);
    return e_f / constants::boltzmann_hartree_per_kelvin;
}

[[nodiscard]] double reduced_temperature_from_kelvin(double t_kelvin, double rs) noexcept
{
    return t_kelvin / fermi_temperature_kelvin(rs);
}

[[nodiscard]] std::filesystem::path resolve_output_path()
{
    const std::filesystem::path relative =
        std::filesystem::path{output_directory} / output_filename;
    if (std::filesystem::exists(relative.parent_path()) ||
        std::filesystem::create_directories(relative.parent_path())) {
        return relative;
    }
    return std::filesystem::path{".."} / relative;
}

[[nodiscard]] std::vector<double> logspace(double t_min, double t_max, std::size_t count)
{
    std::vector<double> values;
    values.reserve(count);
    if (count == 0) {
        return values;
    }
    if (count == 1) {
        values.push_back(t_min);
        return values;
    }
    const double log_min = std::log(t_min);
    const double log_max = std::log(t_max);
    for (std::size_t i = 0; i < count; ++i) {
        const double frac = static_cast<double>(i) / static_cast<double>(count - 1);
        values.push_back(std::exp(log_min + frac * (log_max - log_min)));
    }
    return values;
}

[[nodiscard]] std::vector<double> linspace(double x_min, double x_max, std::size_t count)
{
    std::vector<double> values;
    values.reserve(count);
    if (count == 0) {
        return values;
    }
    if (count == 1) {
        values.push_back(x_min);
        return values;
    }
    for (std::size_t i = 0; i < count; ++i) {
        const double frac = static_cast<double>(i) / static_cast<double>(count - 1);
        values.push_back(x_min + frac * (x_max - x_min));
    }
    return values;
}

}  // namespace

int main()
{
    // Metallic density window r_s ∈ [1, 5]; T from room temperature to extreme cold.
    const auto rs_grid = linspace(1.0, 5.0, 17);
    const auto t_grid = logspace(0.01, 300.0, 49);

    const std::filesystem::path output_path = resolve_output_path();
    std::ofstream output(output_path);
    if (!output) {
        std::cerr << "Failed to open " << output_path << '\n';
        return 1;
    }

    output << std::scientific << std::setprecision(16);
    output << "# MosaiQ-Lindhard (r_s, T) singularity-excision validation sweep\n";
    output << "# Fixed probe: q_bar = " << q_bar << "  omega_bar = " << omega_bar
           << "  gamma = " << gamma_value << " (mu = epsilon_F)\n";
    output << "# T=0 reference: reverse-Dedekind / Stratonovich causality-pinned Lindhard\n";
    output << "# columns: r_s T_K tau re_engine im_engine re_t0 im_t0 "
              "rel_err_re rel_err_im finite_flag\n";

    const WaveVector<> q{q_bar};
    const Frequency<> omega{omega_bar};
    const ReducedChemicalPotential<> gamma_wrapped{gamma_value};
    const ReducedTemperature<> tau_zero{0.0};
    const LindhardResult<> chi_t0 =
        evaluate_lindhard(q, omega, tau_zero, gamma_wrapped);

    std::size_t nan_count = 0;
    double max_rel_re = 0.0;
    double max_rel_im = 0.0;

    for (double rs : rs_grid) {
        for (double t_kelvin : t_grid) {
            const double tau = reduced_temperature_from_kelvin(t_kelvin, rs);
            const LindhardResult<> chi = evaluate_lindhard(
                q, omega, ReducedTemperature<>{tau}, gamma_wrapped);

            const bool finite =
                std::isfinite(chi.real()) && std::isfinite(chi.imag());
            if (!finite) {
                ++nan_count;
            }

            const double denom_re = std::max(std::abs(chi_t0.real()), 1.0e-30);
            const double denom_im = std::max(std::abs(chi_t0.imag()), 1.0e-30);
            const double rel_re = finite ? std::abs(chi.real() - chi_t0.real()) / denom_re
                                         : std::numeric_limits<double>::quiet_NaN();
            const double rel_im = finite ? std::abs(chi.imag() - chi_t0.imag()) / denom_im
                                         : std::numeric_limits<double>::quiet_NaN();

            if (finite) {
                max_rel_re = std::max(max_rel_re, rel_re);
                max_rel_im = std::max(max_rel_im, rel_im);
            }

            output << rs << ' ' << t_kelvin << ' ' << tau << ' ' << chi.real() << ' '
                   << chi.imag() << ' ' << chi_t0.real() << ' ' << chi_t0.imag() << ' '
                   << rel_re << ' ' << rel_im << ' ' << (finite ? 1 : 0) << '\n';
        }
    }

    std::cout << "Wrote " << output_path << '\n';
    std::cout << "grid: n_rs = " << rs_grid.size() << "  n_T = " << t_grid.size() << '\n';
    std::cout << "NaN/Inf cells = " << nan_count << '\n';
    std::cout << "max rel_err_re (finite) = " << max_rel_re << '\n';
    std::cout << "max rel_err_im (finite) = " << max_rel_im << '\n';
    std::cout << "T=0 reference Re = " << chi_t0.real() << "  Im = " << chi_t0.imag() << '\n';

    return nan_count == 0 ? 0 : 2;
}
