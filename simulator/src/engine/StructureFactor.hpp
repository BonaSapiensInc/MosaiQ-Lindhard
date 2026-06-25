#pragma once

#include "core/Concepts.hpp"
#include "engine/RPA.hpp"

#include <complex>
#include <optional>
#include <string>
#include <vector>

namespace mosaiq {

/// Proton mass in electron-mass units (legacy two-component hydrogenic plasma).
inline constexpr double ion_mass_ratio = 1836.15267343;

inline constexpr double default_q_min = 0.1;
inline constexpr double default_q_max = 15.0;
inline constexpr double default_q_step = 0.05;

inline constexpr double default_omega_min = 0.01;
inline constexpr double default_omega_max = 3.0;
inline constexpr double default_omega_step = 0.02;

inline constexpr double static_omega_min = 0.01;
inline constexpr double static_omega_max = 30.0;
inline constexpr double static_omega_step = 0.002;

inline constexpr double coulomb_softening = 0.01;

struct PlasmaSpecies {
    double mass{};
    double k_f{};
    double E_F{};
    ReducedTemperature<> tau{};
    ReducedChemicalPotential<> gamma{};
};

struct PlasmaContext {
    double rs{};
    double n{};
    PlasmaSpecies electron{};
    PlasmaSpecies ion{};
};

struct DynamicStructureFactorSample {
    double omega_e{};
    double S_ee{};
    double S_ii{};
    double S_ei{};
};

struct StaticStructureFactor {
    double S_ee{};
    double S_ii{};
    double S_ei{};
};

[[nodiscard]] double number_density_from_rs(double rs) noexcept;

[[nodiscard]] double fermi_energy_from_kf(double k_f, double mass) noexcept;

[[nodiscard]] double reduced_temperature_from_kelvin(double T_kelvin, double E_F) noexcept;

[[nodiscard]] double temperature_kelvin_from_gamma(double rs, double gamma) noexcept;

[[nodiscard]] BarePotentials<> coulomb_potentials_rational(double q_rational) noexcept;

[[nodiscard]] double bose_einstein_factor(double omega_bar, double tau) noexcept;

[[nodiscard]] double dynamic_structure_factor(const std::complex<double>& chi_rpa,
                                              double omega_bar,
                                              double tau) noexcept;

[[nodiscard]] bool build_plasma_context(double rs, double T_kelvin, PlasmaContext& out) noexcept;

[[nodiscard]] std::optional<std::vector<double>> parse_gamma_list(const std::string& text);

[[nodiscard]] DynamicStructureFactorSample evaluate_dynamic_sample(
    double q,
    double omega_e,
    const PlasmaContext& plasma) noexcept;

[[nodiscard]] StaticStructureFactor integrate_static_structure_factor(
    double q,
    const PlasmaContext& plasma,
    double omega_min,
    double omega_max,
    double omega_step) noexcept;

}  // namespace mosaiq
