#pragma once

#include "core/Concepts.hpp"

#include <complex>

namespace mosaiq {

/// Bare two-component interaction potentials v_st(q) at fixed (q, ω).
template<ScalarPhysical T = double>
struct BarePotentials {
    T v_ee{};
    T v_ii{};
    T v_ei{};
};

/// RPA susceptibilities χ^RPA_st — manuscript Eqs. (response-ee), (response-ii), (response-ei).
template<ScalarPhysical T = double>
struct RpaResult {
    std::complex<T> chi_ee{};
    std::complex<T> chi_ii{};
    std::complex<T> chi_ei{};

    [[nodiscard]] constexpr const std::complex<T>& ee() const noexcept { return chi_ee; }
    [[nodiscard]] constexpr const std::complex<T>& ii() const noexcept { return chi_ii; }
    [[nodiscard]] constexpr const std::complex<T>& ei() const noexcept { return chi_ei; }
};

/// Map Lindhard Re/Im pair to complex χ^L(q, ω).
template<ScalarPhysical T>
[[nodiscard]] constexpr std::complex<T> as_complex(LindhardResult<T> chi) noexcept
{
    return {chi.real(), chi.imag()};
}

/// D(ε_F) = 3n / (2 E_F) — manuscript Eq. (dos-at-Fermi-energy); ℏ = 1.
[[nodiscard]] constexpr double density_of_states_at_fermi(double n, double E_F) noexcept
{
    return 1.5 * n / E_F;
}

/// Legacy complexLindhard convention: χ_nat = D(E_F) × χ̃_reduced.
template<ScalarPhysical T>
[[nodiscard]] constexpr std::complex<T> susceptibility_in_natural_units(
    LindhardResult<T> chi_reduced,
    T density_of_states) noexcept
{
    return as_complex(chi_reduced) * density_of_states;
}

/// Two-component dielectric function ε(q, ω) — manuscript Eq. (dielectric-function) / legacy `get_dielectric`.
template<ScalarPhysical T = double>
[[nodiscard]] std::complex<T> evaluate_dielectric(std::complex<T> chi_e,
                                                std::complex<T> chi_i,
                                                BarePotentials<T> potentials);

template<ScalarPhysical T = double>
[[nodiscard]] std::complex<T> evaluate_dielectric(LindhardResult<T> chi_e,
                                                LindhardResult<T> chi_i,
                                                BarePotentials<T> potentials);

/// RPA susceptibility matrix {χ^RPA_ee, χ^RPA_ii, χ^RPA_ei} — manuscript Eqs. (response-ee)–(response-ei).
template<ScalarPhysical T = double>
[[nodiscard]] RpaResult<T> evaluate_rpa_susceptibility(std::complex<T> chi_e,
                                                       std::complex<T> chi_i,
                                                       BarePotentials<T> potentials);

template<ScalarPhysical T = double>
[[nodiscard]] RpaResult<T> evaluate_rpa_susceptibility(LindhardResult<T> chi_e,
                                                       LindhardResult<T> chi_i,
                                                       BarePotentials<T> potentials);

/// Assemble RPA from reduced Lindhard outputs after per-species DOS restoration.
template<ScalarPhysical T = double>
[[nodiscard]] RpaResult<T> evaluate_rpa_susceptibility(LindhardResult<T> chi_e,
                                                       LindhardResult<T> chi_i,
                                                       BarePotentials<T> potentials,
                                                       T electron_density_of_states,
                                                       T ion_density_of_states);

}  // namespace mosaiq
