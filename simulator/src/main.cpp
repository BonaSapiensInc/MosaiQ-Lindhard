#include "core/Concepts.hpp"
#include "engine/Lindhard.hpp"
#include "engine/RPA.hpp"
#include "physics/ChemicalPotential.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

namespace {

using namespace sophus;

/// Proton mass in electron-mass units (legacy two-component hydrogenic plasma).
inline constexpr double ion_mass_ratio = 1836.15267343;

inline constexpr double q_min = 0.1;
inline constexpr double q_max = 3.0;
inline constexpr double q_step = 0.5;

inline constexpr double omega_min = 0.01;
inline constexpr double omega_max = 2.0;
inline constexpr double omega_step = 0.25;

inline constexpr double coulomb_softening = 0.01;

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
        v_coulomb,   // v_ee  (|e|²)
        v_coulomb,   // v_ii  (|Z_i|²)
        -v_coulomb,  // v_ei  (e · Z_i < 0 for opposite charges)
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

/// Fluctuation–dissipation theorem: S_ee ∝ −Im χ_ee^RPA / (1 − exp(−ω̄/τ)).
[[nodiscard]] double dynamic_structure_factor_ee(
    const std::complex<double>& chi_ee_rpa,
    double omega_bar_e,
    double tau_e) noexcept
{
    const double denominator = bose_einstein_factor(omega_bar_e, tau_e);
    if (!std::isfinite(denominator) || std::abs(denominator) < 1.0e-300) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return -chi_ee_rpa.imag() / denominator;
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

    std::ofstream output("output_structure_factor.dat");
    if (!output) {
        std::cerr << "Error: cannot open output_structure_factor.dat for writing.\n";
        return EXIT_FAILURE;
    }

    output << std::scientific << std::setprecision(12);
    output << "# SOPHUS-Lindhard CLI — RPA electron–ion structure factor\n";
    output << "# r_s = " << rs << "  T_K = " << T_kelvin << '\n';
    output << "# m_e = 1  m_i/m_e = " << ion_mass_ratio << '\n';
    output << "# k_F = " << electron.k_f << " a_0^-1  n = " << number_density_from_rs(rs)
           << " a_0^-3\n";
    output << "# gamma_e = " << electron.gamma.raw() << "  tau_e = " << electron.tau.raw()
           << "  gamma_i = " << ion.gamma.raw() << "  tau_i = " << ion.tau.raw() << '\n';
    output << "# columns: q  omega  Re(chi_ee^RPA)  Im(chi_ee^RPA)  S_ee\n";

    std::size_t rows_written = 0;
    const std::size_t total_q =
        static_cast<std::size_t>(std::floor((q_max - q_min) / q_step)) + 1;
    const std::size_t total_omega =
        static_cast<std::size_t>(std::floor((omega_max - omega_min) / omega_step)) + 1;
    std::cerr << "Scanning " << total_q << " x " << total_omega << " (q, omega) grid...\n";

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

            const double S_ee = dynamic_structure_factor_ee(rpa.chi_ee, omega_e, electron.tau.raw());
            if (!std::isfinite(S_ee)) {
                continue;
            }

            output << q << ' ' << omega_e << ' ' << rpa.chi_ee.real() << ' ' << rpa.chi_ee.imag()
                   << ' ' << S_ee << '\n';
            ++rows_written;
        }
    }

    if (rows_written == 0) {
        std::cerr << "Error: no finite grid points were written.\n";
        return EXIT_FAILURE;
    }

    std::cout << "Wrote " << rows_written << " rows to output_structure_factor.dat\n";
    std::cout << "  r_s = " << rs << ", T = " << T_kelvin << " K\n";
    std::cout << "  tau_e = " << electron.tau.raw() << ", tau_i = " << ion.tau.raw() << '\n';
    return EXIT_SUCCESS;
}
