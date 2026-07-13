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

#include <numbers>

namespace mosaiq::constants {

/// Mathematical constants
inline constexpr double pi = std::numbers::pi;
inline constexpr double two_pi = 2.0 * pi;
inline constexpr double four_pi = 4.0 * pi;

/// Legacy: `Component.cpp` — `_pi(4.0 * atan(1.0))`
inline constexpr double legacy_pi = pi;

/// Legacy: `TemperatureFermionComponent`, `Plasma.cpp` — `_constantBoltzmann(3.1668114e-6)`
/// Boltzmann constant in Hartree / kelvin.
inline constexpr double boltzmann_hartree_per_kelvin = 3.1668114e-6;

/// Legacy: `Potential.cpp` — `_Ehe(-0.90767486e-1)` statvolts expressed as E_H / e.
inline constexpr double statvolt_hartree_per_e = -0.90767486e-1;

/// Atomic units (Hartree system) — canonical reference values for Phase 2+ thermodynamics.
inline constexpr double electron_mass_amu = 1.0;  // m_e as unit mass
inline constexpr double hartree_energy = 1.0;     // 1 Ha
inline constexpr double bohr_radius = 1.0;        // 1 a_0

/// Legacy Debye prefactor fragment: `sqrt(3 Z^2 * 311775 / (T rs^3))` in TemperatureFermionComponent.
inline constexpr double debye_prefactor_numerator = 311775.0;

/// Sinc-quadrature strip parameter d = h * 7 (legacy: `sqrt(pi * 7 / N)` for step h).
inline constexpr double sinc_strip_factor = 7.0;

/// Default GV grand-canonical sinc node count (legacy: `get_GVrealLindhard`, N = 256).
inline constexpr std::size_t default_gv_sinc_nodes = 256;

/// Default KK / plasma-dispersion sinc node count (legacy: `get_KKrealLindhard`, N = 128).
inline constexpr std::size_t default_kk_sinc_nodes = 128;

/// Below this reduced temperature, Kramers–Krönig Re χ^L uses reverse-Dedekind
/// singularity excision (analytic Stratonovich step + smooth residual sinc).
/// Above it, the legacy bit-reproducible softplus/log path is retained.
inline constexpr double singularity_excision_tau_threshold = 2.0e-2;

/// (9 pi / 4)^(1/3) — legacy FermionComponent fermi wavevector prefactor.
inline constexpr double nine_pi_over_four_cbrt = 1.7044777697860926;

/// Legacy Fermi wavevector convention: (9 pi / 4)^(1/3) / r_s  [a_0^-1].
[[nodiscard]] constexpr double fermi_wavevector_from_rs(double rs) noexcept
{
    return nine_pi_over_four_cbrt / rs;
}

/// Convert absolute temperature [K] to Hartree energy (legacy: `_tempHartree = k_B * T`).
[[nodiscard]] constexpr double temperature_kelvin_to_hartree(double T_kelvin) noexcept
{
    return boltzmann_hartree_per_kelvin * T_kelvin;
}

}  // namespace mosaiq::constants
