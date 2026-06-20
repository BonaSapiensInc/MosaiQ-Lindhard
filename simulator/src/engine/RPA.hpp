#pragma once

#include "core/Concepts.hpp"

#include <complex>

namespace sophus {

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

}  // namespace sophus
