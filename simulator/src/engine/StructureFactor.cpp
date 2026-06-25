#include "engine/StructureFactor.hpp"

#include "core/TrapezoidalIntegration.hpp"
#include "engine/Lindhard.hpp"
#include "physics/ChemicalPotential.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <limits>
#include <sstream>

namespace mosaiq {

double number_density_from_rs(double rs) noexcept
{
    return 3.0 / (4.0 * constants::pi * rs * rs * rs);
}

double fermi_energy_from_kf(double k_f, double mass) noexcept
{
    return k_f * k_f / (2.0 * mass);
}

double reduced_temperature_from_kelvin(double T_kelvin, double E_F) noexcept
{
    return constants::temperature_kelvin_to_hartree(T_kelvin) / E_F;
}

double temperature_kelvin_from_gamma(double rs, double gamma) noexcept
{
    return constants::debye_prefactor_numerator / (rs * gamma);
}

BarePotentials<> coulomb_potentials_rational(double q_rational) noexcept
{
    const double q_sq = q_rational * q_rational + coulomb_softening * coulomb_softening;
    const double v_coulomb = constants::four_pi / q_sq;
    return BarePotentials<>{
        v_coulomb,
        v_coulomb,
        -v_coulomb,
    };
}

double bose_einstein_factor(double omega_bar, double tau) noexcept
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
        // ω̄ ≫ τ: 1 − e^{−ω̄/τ} → 1 (Bose denominator in fluctuation–dissipation).
        return 1.0;
    }
    return 1.0 - std::exp(exponent);
}

double dynamic_structure_factor(const std::complex<double>& chi_rpa,
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

namespace {

[[nodiscard]] bool initialize_species(double rs,
                                      double T_kelvin,
                                      double mass,
                                      PlasmaSpecies& out) noexcept
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

bool build_plasma_context(double rs, double T_kelvin, PlasmaContext& out) noexcept
{
    out.rs = rs;
    out.n = number_density_from_rs(rs);
    if (!initialize_species(rs, T_kelvin, constants::electron_mass_amu, out.electron)) {
        return false;
    }
    if (!initialize_species(rs, T_kelvin, ion_mass_ratio, out.ion)) {
        return false;
    }
    return true;
}

std::optional<std::vector<double>> parse_gamma_list(const std::string& text)
{
    std::vector<double> values;
    std::stringstream stream(text);
    std::string token;

    while (std::getline(stream, token, ',')) {
        if (token.empty()) {
            continue;
        }

        char* end = nullptr;
        const double value = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || !std::isfinite(value) || value <= 0.0) {
            return std::nullopt;
        }
        values.push_back(value);
    }

    if (values.empty()) {
        return std::nullopt;
    }
    return values;
}

DynamicStructureFactorSample evaluate_dynamic_sample(double q,
                                                     double omega_e,
                                                     const PlasmaContext& plasma) noexcept
{
    const RpaResult<> rpa = evaluate_rpa_response(q, omega_e, plasma);
    const double omega_i = omega_e * plasma.electron.E_F / plasma.ion.E_F;

    return DynamicStructureFactorSample{
        omega_e,
        dynamic_structure_factor(rpa.chi_ee, omega_e, plasma.electron.tau.raw()),
        dynamic_structure_factor(rpa.chi_ii, omega_i, plasma.ion.tau.raw()),
        dynamic_structure_factor(rpa.chi_ei, omega_e, plasma.electron.tau.raw()),
    };
}

RpaResult<> evaluate_rpa_response(double q,
                                  double omega_e,
                                  const PlasmaContext& plasma) noexcept
{
    const double omega_i = omega_e * plasma.electron.E_F / plasma.ion.E_F;
    const double q_i = q * (plasma.electron.k_f / plasma.ion.k_f);
    const BarePotentials<> potentials = coulomb_potentials_rational(q);

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

    return evaluate_rpa_susceptibility(chi_e, chi_i, potentials);
}

StaticStructureFactor integrate_static_structure_factor(double q,
                                                        const PlasmaContext& plasma,
                                                        double omega_min,
                                                        double omega_max,
                                                        double omega_step) noexcept
{
    std::vector<double> omega_grid;
    std::vector<double> s_ee;
    std::vector<double> s_ii;
    std::vector<double> s_ei;

    // Calculate the dynamic upper bound based on physical kinematics: omega_max = q^2 + 2q * sqrt(gamma + 36 * tau)
    // We ensure it does not fall below the baseline static omega_max for small q.
    const double dynamic_limit =
        (q * q) + 2.0 * q * std::sqrt(std::max(0.0, plasma.electron.gamma.raw()) +
                                      36.0 * plasma.electron.tau.raw());
    const double actual_omega_max = std::max(omega_max, dynamic_limit);

    for (double omega_e = omega_min; omega_e <= actual_omega_max + 0.5 * omega_step;
         omega_e += omega_step) {
        const DynamicStructureFactorSample sample =
            evaluate_dynamic_sample(q, omega_e, plasma);

        if (!std::isfinite(sample.S_ee) || !std::isfinite(sample.S_ii) ||
            !std::isfinite(sample.S_ei)) {
            continue;
        }

        omega_grid.push_back(omega_e);
        s_ee.push_back(sample.S_ee);
        s_ii.push_back(sample.S_ii * (plasma.electron.E_F / plasma.ion.E_F));
        s_ei.push_back(sample.S_ei);
    }

    if (omega_grid.size() < 2) {
        return StaticStructureFactor{
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
        };
    }

    return StaticStructureFactor{
        trapezoidal_integrate(omega_grid, s_ee),
        trapezoidal_integrate(omega_grid, s_ii),
        trapezoidal_integrate(omega_grid, s_ei),
    };
}

}  // namespace mosaiq
