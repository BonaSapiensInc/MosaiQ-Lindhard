/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#pragma once

#include "core/Concepts.hpp"
#include "engine/RPA.hpp"
#include "engine/ZetaRPA.hpp"

#include <complex>
#include <optional>
#include <string>
#include <vector>

namespace mosaiq {

/// Proton mass in electron-mass units (legacy two-component hydrogenic plasma).
inline constexpr double ion_mass_ratio = 1836.15267343;

inline constexpr double default_q_min = 0.1;
inline constexpr double default_q_max = 50.0;
inline constexpr double default_q_step = 0.05;

/// Upper $\bar{q}$ limit for the dynamic $S(\bar{q},\bar{\omega})$ mesh (Figure 3).
/// Electron panels crop to $\bar{q}\le 4$ in plotting; $S_{ii}$ uses the full macroscopic extent.
inline constexpr double dynamic_sq_q_max = 50.0;

inline constexpr double default_omega_min = 0.01;
inline constexpr double default_omega_max = 3.0;
inline constexpr double default_omega_step = 0.02;

/// Dynamic $S_{ee}$ export ceiling for Figure~5 plasmon background ($\bar{\omega} \le 6.5$).
inline constexpr double structure_factor_export_omega_max = 6.5;

inline constexpr double static_omega_min = 0.01;
inline constexpr double static_omega_max = 30.0;
inline constexpr double static_omega_step = 0.002;

inline constexpr double coulomb_softening = 0.01;

/// Upper ceiling for the dynamic $S(\bar{q},\bar{\omega})$ mesh at $\bar{q}=50$ ($\bar{\omega}\approx\bar{q}^2$).
inline constexpr double dynamic_sq_omega_ceiling = 2600.0;

/// Target number of $\bar{\omega}$ nodes per $\bar{q}$ on the wide-range dynamic mesh.
inline constexpr double dynamic_sq_omega_mesh_target = 400.0;

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

/// Plasma coupling Γ = 311775 / (r_s T[K]) — inverse of temperature_kelvin_from_gamma.
[[nodiscard]] double plasma_coupling_from_kelvin(double rs, double T_kelvin) noexcept;

/// Recover absolute temperature [K] from electron reduced thermodynamics.
[[nodiscard]] double temperature_kelvin_from_plasma(const PlasmaContext& plasma) noexcept;

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
    const PlasmaContext& plasma,
    ResponsePathway pathway = ResponsePathway::ZetaRPA) noexcept;

/// Lindhard (reduced) → DOS-restored response at (q̄_e, ω̄_e).
/// Default production pathway is multi-component Zeta-RPA (Phase Z4 versioned switch).
[[nodiscard]] ZetaRpaMatrixResult<> evaluate_rpa_response(
    double q,
    double omega_e,
    const PlasmaContext& plasma,
    ResponsePathway pathway = ResponsePathway::ZetaRPA) noexcept;

[[nodiscard]] StaticStructureFactor integrate_static_structure_factor(
    double q,
    const PlasmaContext& plasma,
    double omega_min,
    double omega_max,
    double omega_step,
    ResponsePathway pathway = ResponsePathway::ZetaRPA) noexcept;

/// Kinematic $\bar{\omega}$ upper bound for the wide-range dynamic structure-factor mesh.
[[nodiscard]] double dynamic_structure_factor_omega_max(double q,
                                                          const PlasmaContext& plasma) noexcept;

/// Adaptive $\bar{\omega}$ step keeping the mesh near \c dynamic_sq_omega_mesh_target nodes.
[[nodiscard]] double dynamic_structure_factor_omega_step(double omega_max) noexcept;

}  // namespace mosaiq
