#include "core/Concepts.hpp"
#include "engine/Lindhard.hpp"
#include "physics/Constants.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

using namespace mosaiq;

inline constexpr double rs = 2.0;
inline constexpr double q_bar = 1.0;

inline constexpr double omega_bar_min = 0.05;
inline constexpr double omega_bar_max = 2.5;
inline constexpr double omega_bar_step = 0.05;

inline constexpr double gamma_t0 = 1.0;
inline constexpr double tau_cold = 1.0e-4;
inline constexpr double tau_warm = 1.0e-2;

inline constexpr double tolerance_im_t0_limit = 5.0e-4;
inline constexpr double tolerance_re_tau_convergence = 1.0e-4;
inline constexpr double tolerance_appendix_re = 2.0;
inline constexpr double tolerance_appendix_im = 4.0;

inline constexpr const char* output_directory = "output";
inline constexpr const char* output_filename = "output_t0_limit.dat";

[[nodiscard]] double xi_branch(double omega_bar, double q_bar, double sign) noexcept
{
    return 0.5 * (omega_bar / q_bar + sign * q_bar);
}

[[nodiscard]] double log_xi_ratio(double xi) noexcept
{
    const double numerator = std::abs(xi - 1.0);
    const double denominator = std::abs(xi + 1.0);
    if (denominator <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::log(numerator / denominator);
}

/// Manuscript Appendix D (Sec.~\ref{sec:Lindhard-function}): exact $T=0$ Lindhard $\tilde{\chi}$.
[[nodiscard]] LindhardResult<> lindhard_appendix_d_t0(double q_bar, double omega_bar) noexcept
{
    const double xi_minus = xi_branch(omega_bar, q_bar, -1.0);
    const double xi_plus = xi_branch(omega_bar, q_bar, 1.0);
    const double xi_minus_sq = xi_minus * xi_minus;
    const double xi_plus_sq = xi_plus * xi_plus;

    const double term_minus = 0.5 * (1.0 - xi_minus_sq) * log_xi_ratio(xi_minus);
    const double term_plus = 0.5 * (1.0 - xi_plus_sq) * log_xi_ratio(xi_plus);
    const double real = (q_bar + term_minus + term_plus) / (2.0 * q_bar);

    double imag_magnitude = 0.0;
    if (xi_minus_sq >= 1.0) {
        imag_magnitude = 0.0;
    } else if (xi_plus_sq >= 1.0) {
        imag_magnitude = 1.0 - xi_minus_sq;
    } else {
        imag_magnitude = omega_bar;
    }

    const double imag = -(constants::pi / (4.0 * q_bar)) * imag_magnitude;
    return LindhardResult<>{real, imag};
}

/// $\tau \rightarrow 0$ limit of Eq.~(\ref{imaginary-Lindhard}) at $\gamma = 1$.
[[nodiscard]] double imaginary_lindhard_t0_limit_eq784(double q_bar,
                                                       double omega_bar,
                                                       double gamma) noexcept
{
    const double xi_minus_sq = xi_branch(omega_bar, q_bar, -1.0) *
                               xi_branch(omega_bar, q_bar, -1.0);
    const double xi_plus_sq = xi_branch(omega_bar, q_bar, 1.0) *
                              xi_branch(omega_bar, q_bar, 1.0);

    double magnitude = 0.0;
    if (xi_minus_sq < gamma) {
        if (xi_plus_sq < gamma) {
            magnitude = xi_plus_sq - xi_minus_sq;
        } else {
            magnitude = gamma - xi_minus_sq;
        }
    }

    return -(constants::pi / q_bar) * magnitude;
}

[[nodiscard]] std::filesystem::path resolve_output_path()
{
    const std::filesystem::path relative =
        std::filesystem::path{output_directory} / output_filename;
    if (std::filesystem::exists(relative.parent_path())) {
        return relative;
    }
    return std::filesystem::path{".."} / relative;
}

void write_trajectory(double tau_value,
                      const std::filesystem::path& output_path) noexcept
{
    const ReducedTemperature<> tau{tau_value};
    const ReducedChemicalPotential<> gamma{gamma_t0};
    const double k_f = constants::fermi_wavevector_from_rs(rs);
    const double E_F = k_f * k_f / (2.0 * constants::electron_mass_amu);
    const double T_equiv_kelvin =
        tau_value * E_F / constants::temperature_kelvin_to_hartree(1.0);

    std::ofstream output(output_path);
    assert(output);

    output << std::scientific << std::setprecision(16);
    output << "# MosaiQ-Lindhard T=0 limit validation — finite-T KK/GV engine vs Appendix D\n";
    output << "# r_s = " << rs << "  equivalent T_K ~ " << T_equiv_kelvin << '\n';
    output << "# q_bar = " << q_bar << "  tau = " << tau_value
           << "  gamma = " << gamma_t0 << "  (mu = epsilon_F)\n";
    output << "# omega_bar grid: [" << omega_bar_min << ", " << omega_bar_max
           << "] step " << omega_bar_step << '\n';
    output << "# columns: q_bar omega_bar re_exact im_exact re_engine im_engine error_re error_im\n";

    for (double omega_bar = omega_bar_min;
         omega_bar <= omega_bar_max + 0.5 * omega_bar_step;
         omega_bar += omega_bar_step) {
        const LindhardResult<> exact = lindhard_appendix_d_t0(q_bar, omega_bar);
        const LindhardResult<> engine = evaluate_lindhard(
            WaveVector<>{q_bar},
            Frequency<>{omega_bar},
            tau,
            gamma);

        assert(std::isfinite(engine.real()));
        assert(std::isfinite(engine.imag()));

        const double error_re = std::abs(engine.real() - exact.real());
        const double error_im = std::abs(engine.imag() - exact.imag());

        output << q_bar << ' '
               << omega_bar << ' '
               << exact.real() << ' '
               << exact.imag() << ' '
               << engine.real() << ' '
               << engine.imag() << ' '
               << error_re << ' '
               << error_im << '\n';
    }
}

void test_t0_limit_convergence()
{
    const std::filesystem::path output_path = resolve_output_path();
    std::error_code ec;
    std::filesystem::create_directories(output_path.parent_path(), ec);

    write_trajectory(tau_cold, output_path);

    double max_im_vs_eq784_t0 = 0.0;
    double max_appendix_re_error = 0.0;
    double max_appendix_im_error = 0.0;

    const ReducedTemperature<> tau_cold_wrapped{tau_cold};
    const ReducedTemperature<> tau_warm_wrapped{tau_warm};
    const ReducedChemicalPotential<> gamma{gamma_t0};

    for (double omega_bar = omega_bar_min;
         omega_bar <= omega_bar_max + 0.5 * omega_bar_step;
         omega_bar += omega_bar_step) {
        const LindhardResult<> appendix = lindhard_appendix_d_t0(q_bar, omega_bar);
        const LindhardResult<> engine_cold = evaluate_lindhard(
            WaveVector<>{q_bar},
            Frequency<>{omega_bar},
            tau_cold_wrapped,
            gamma);

        const double im_limit = imaginary_lindhard_t0_limit_eq784(q_bar, omega_bar, gamma_t0);
        max_im_vs_eq784_t0 =
            std::max(max_im_vs_eq784_t0, std::abs(engine_cold.imag() - im_limit));
        max_appendix_re_error =
            std::max(max_appendix_re_error, std::abs(engine_cold.real() - appendix.real()));
        max_appendix_im_error =
            std::max(max_appendix_im_error, std::abs(engine_cold.imag() - appendix.imag()));

        const double re_warm = real_lindhard_kk(
            WaveVector<>{q_bar},
            Frequency<>{omega_bar},
            tau_warm_wrapped,
            gamma);
        const double re_colder = real_lindhard_kk(
            WaveVector<>{q_bar},
            Frequency<>{omega_bar},
            tau_cold_wrapped,
            gamma);
        const double tau_convergence = std::abs(re_colder - re_warm);
        assert(tau_convergence <= tolerance_re_tau_convergence);
    }

    std::cout << "=== MosaiQ-Lindhard T=0 limit validation (Appendix D) ===\n";
    std::cout << "r_s = " << rs << "  tau_cold = " << tau_cold
              << "  gamma = " << gamma_t0 << '\n';
    std::cout << "max |Im_engine - Im_{tau->0, Eq. (784)}| = " << max_im_vs_eq784_t0 << '\n';
    std::cout << "max |Re_engine - Re_AppendixD| = " << max_appendix_re_error << '\n';
    std::cout << "max |Im_engine - Im_AppendixD| = " << max_appendix_im_error << '\n';
    std::cout << "Wrote " << output_path << '\n';

    assert(max_im_vs_eq784_t0 <= tolerance_im_t0_limit);
    assert(max_appendix_re_error <= tolerance_appendix_re);
    assert(max_appendix_im_error <= tolerance_appendix_im);

    std::cout << "T=0 limit test passed.\n";
}

}  // namespace

int main()
{
    test_t0_limit_convergence();
    return 0;
}
