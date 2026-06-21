#include "core/Concepts.hpp"
#include "engine/Lindhard.hpp"
#include "engine/PlasmonPoleExtractor.hpp"
#include "engine/RPA.hpp"
#include "physics/ChemicalPotential.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>

namespace {

using namespace sophus;

/// Proton mass in electron-mass units (legacy two-component hydrogenic plasma).
inline constexpr double ion_mass_ratio = 1836.15267343;

inline constexpr double q_min = 0.1;
inline constexpr double q_max = 4.0;
inline constexpr double q_step = 0.05;

inline constexpr double omega_min = 0.01;
inline constexpr double omega_max = 3.0;
inline constexpr double omega_step = 0.02;

inline constexpr double coulomb_softening = 0.01;

inline constexpr const char* output_directory = "output";
inline constexpr const char* structure_factor_filename = "output_structure_factor.dat";
inline constexpr const char* dispersion_filename = "output_plasmon_dispersion.dat";

void print_usage(const char* program)
{
    std::cerr << "Usage: " << program << " <r_s> <T_kelvin>\n"
              << "  r_s       Wigner-Seitz radius in Bohr radii (e.g. 1.0)\n"
              << "  T_kelvin  temperature in kelvin (e.g. 10000.0)\n";
}

[[nodiscard]] double number_density_from_rs(double rs) noexcept
{
    return 3.0 / (4.0 * constants::pi * rs * rs * rs);
}

[[nodiscard]] double fermi_energy_from_kf(double k_f, double mass) noexcept
{
    return k_f * k_f / (2.0 * mass);
}

[[nodiscard]] double reduced_temperature_from_kelvin(double T_kelvin, double E_F) noexcept
{
    const double temp_hartree = constants::temperature_kelvin_to_hartree(T_kelvin);
    return temp_hartree / E_F;
}

[[nodiscard]] BarePotentials<> coulomb_potentials_rational(double q_rational) noexcept
{
    const double q_sq = q_rational * q_rational + coulomb_softening * coulomb_softening;
    const double v_coulomb = constants::four_pi / q_sq;
    return BarePotentials<>{
        v_coulomb,
        v_coulomb,
        -v_coulomb,
    };
}

/// Bose–Einstein factor 1 − exp(−ω̄/τ) with small-ω̄ thermal limit to avoid 0/0.
[[nodiscard]] double bose_einstein_factor(double omega_bar, double tau) noexcept
{
    if (tau <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::abs(omega_bar) < 1.0e-14) {
        return omega_bar / tau;
    }
    const double exponent = -omega_bar / tau;
    if (exponent > 700.0) {
        return 1.0;
    }
    if (exponent < -700.0) {
        return -std::exp(exponent);
    }
    return 1.0 - std::exp(exponent);
}

/// Fluctuation–dissipation theorem: S_st ∝ −Im χ_st^RPA / (1 − exp(−ω̄/τ)).
[[nodiscard]] double dynamic_structure_factor(
    const std::complex<double>& chi_rpa,
    double omega_bar,
    double tau) noexcept
{
    const double denominator = bose_einstein_factor(omega_bar, tau);
    if (!std::isfinite(denominator) || std::abs(denominator) < 1.0e-300) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double imag = chi_rpa.imag();
    if (!std::isfinite(imag)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return -imag / denominator;
}

struct ComponentState {
    double mass{};
    double k_f{};
    double E_F{};
    ReducedTemperature<> tau{};
    ReducedChemicalPotential<> gamma{};
};

[[nodiscard]] bool initialize_component(double rs,
                                        double T_kelvin,
                                        double mass,
                                        ComponentState& out)
{
    const double density = number_density_from_rs(rs);
    out.mass = mass;
    out.k_f = constants::fermi_wavevector_from_rs(rs);
    out.E_F = fermi_energy_from_kf(out.k_f, mass);
    out.tau = ReducedTemperature<>{reduced_temperature_from_kelvin(T_kelvin, out.E_F)};

    const auto gamma = find_chemical_potential(ChemicalPotentialInputs<>{
        out.tau,
        Density<>{density},
        WaveVector<>{out.k_f},
    });
    if (!gamma) {
        return false;
    }
    out.gamma = *gamma;
    return true;
}

[[nodiscard]] std::size_t grid_node_count(double min_val, double max_val, double step) noexcept
{
    return static_cast<std::size_t>(std::floor((max_val - min_val) / step)) + 1;
}

[[nodiscard]] PlasmonPoleInputs<> make_plasmon_inputs(double rs,
                                                      double q,
                                                      const ComponentState& electron,
                                                      const ComponentState& ion) noexcept
{
    return PlasmonPoleInputs<>{
        .q = WaveVector<>{q},
        .electron = PlasmonSpeciesParams<>{
            .tau = electron.tau,
            .gamma = electron.gamma,
            .fermi_energy = electron.E_F,
            .mass = electron.mass,
        },
        .ion = PlasmonSpeciesParams<>{
            .tau = ion.tau,
            .gamma = ion.gamma,
            .fermi_energy = ion.E_F,
            .mass = ion.mass,
        },
        .potentials = coulomb_potentials_rational(q),
        .wigner_seitz_radius = rs,
    };
}

void write_run_header(std::ostream& out,
                      double rs,
                      double T_kelvin,
                      const ComponentState& electron,
                      const ComponentState& ion)
{
    out << std::scientific << std::setprecision(12);
    out << "# r_s = " << rs << "  T_K = " << T_kelvin << '\n';
    out << "# m_e = 1  m_i/m_e = " << ion_mass_ratio << '\n';
    out << "# k_F = " << electron.k_f << " a_0^-1  n = " << number_density_from_rs(rs)
        << " a_0^-3\n";
    out << "# gamma_e = " << electron.gamma.raw() << "  tau_e = " << electron.tau.raw()
        << "  gamma_i = " << ion.gamma.raw() << "  tau_i = " << ion.tau.raw() << '\n';
}

[[nodiscard]] std::size_t sweep_plasmon_dispersion(std::ostream& dispersion_out,
                                                   double rs,
                                                   const ComponentState& electron,
                                                   const ComponentState& ion)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::size_t roots_found = 0;

    for (double q = q_min; q <= q_max + 0.5 * q_step; q += q_step) {
        const PlasmonPoleInputs<> inputs = make_plasmon_inputs(rs, q, electron, ion);
        const double bohm_gross = PlasmonPoleExtractor::bohm_gross_frequency(
            WaveVector<>{q},
            electron.tau,
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

}  // namespace

int main(int argc, char* argv[])
{
    if (argc != 3) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

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

    ComponentState electron{};
    ComponentState ion{};
    if (!initialize_component(rs, T_kelvin, constants::electron_mass_amu, electron)) {
        std::cerr << "Error: failed to invert chemical potential for electrons.\n";
        return EXIT_FAILURE;
    }
    if (!initialize_component(rs, T_kelvin, ion_mass_ratio, ion)) {
        std::cerr << "Error: failed to invert chemical potential for ions.\n";
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

    structure_factor_output << "# SOPHUS-Lindhard CLI — multi-component RPA structure factors\n";
    write_run_header(structure_factor_output, rs, T_kelvin, electron, ion);
    structure_factor_output << "# omega column: electron reduced units [E_F/e/hbar]\n";
    structure_factor_output << "# S_ei uses tau_e (electron thermal reference)\n";
    structure_factor_output << "# columns: q omega Im(chi_ee) S_ee Im(chi_ii) S_ii Im(chi_ei) S_ei\n";

    dispersion_output << "# SOPHUS-Lindhard CLI — plasmon dispersion Re[epsilon]=0 trajectory\n";
    write_run_header(dispersion_output, rs, T_kelvin, electron, ion);
    dispersion_output << "# omega_p column: electron reduced units [E_F/e/hbar]\n";
    dispersion_output << "# landau_damping = Im[epsilon(q, omega_p)] at the extracted pole\n";
    dispersion_output << "# bohm_gross_estimate: long-wavelength Bohm-Gross anchor (same units)\n";
    dispersion_output << "# missing roots (continuum / no bracket) are written as NaN NaN\n";
    dispersion_output
        << "# columns: q omega_p_electron_units landau_damping bohm_gross_estimate\n";

    const std::size_t total_q = grid_node_count(q_min, q_max, q_step);
    const std::size_t total_omega = grid_node_count(omega_min, omega_max, omega_step);

    std::cerr << "Tracing plasmon dispersion over " << total_q << " q points...\n";
    const std::size_t dispersion_roots =
        sweep_plasmon_dispersion(dispersion_output, rs, electron, ion);
    std::cerr << "  Plasmon roots found: " << dispersion_roots << " / " << total_q << '\n';

    std::cerr << "Scanning " << total_q << " x " << total_omega << " (q, omega) grid...\n";

    std::size_t rows_written = 0;
    for (double q = q_min; q <= q_max + 0.5 * q_step; q += q_step) {
        const BarePotentials<> potentials = coulomb_potentials_rational(q);

        for (double omega_e = omega_min; omega_e <= omega_max + 0.5 * omega_step;
             omega_e += omega_step) {
            const double omega_i = omega_e * electron.E_F / ion.E_F;

            const LindhardResult<> chi_e = evaluate_lindhard(
                WaveVector<>{q},
                Frequency<>{omega_e},
                electron.tau,
                electron.gamma);

            const LindhardResult<> chi_i = evaluate_lindhard(
                WaveVector<>{q},
                Frequency<>{omega_i},
                ion.tau,
                ion.gamma);

            const RpaResult<> rpa = evaluate_rpa_susceptibility(chi_e, chi_i, potentials);

            const double S_ee =
                dynamic_structure_factor(rpa.chi_ee, omega_e, electron.tau.raw());
            const double S_ii =
                dynamic_structure_factor(rpa.chi_ii, omega_i, ion.tau.raw());
            const double S_ei =
                dynamic_structure_factor(rpa.chi_ei, omega_e, electron.tau.raw());

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

    std::cout << "Wrote " << rows_written << " rows to " << structure_factor_path << '\n';
    std::cout << "Wrote " << total_q << " dispersion rows to " << dispersion_path << " ("
              << dispersion_roots << " roots)\n";
    std::cout << "  r_s = " << rs << ", T = " << T_kelvin << " K\n";
    std::cout << "  tau_e = " << electron.tau.raw() << ", tau_i = " << ion.tau.raw() << '\n';
    return EXIT_SUCCESS;
}
