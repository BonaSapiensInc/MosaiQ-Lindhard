/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "engine/RPA.hpp"

#include "core/Concepts.hpp"

#include <complex>

namespace mosaiq {

namespace {

template<ScalarPhysical T>
[[nodiscard]] std::complex<T> unity() noexcept
{
    return {T{1}, T{0}};
}

template<ScalarPhysical T>
[[nodiscard]] std::complex<T> dielectric_from_susceptibilities(std::complex<T> chi_e,
                                                               std::complex<T> chi_i,
                                                               BarePotentials<T> v) noexcept
{
    // ε = (1 − χ_e v_ee)(1 − χ_i v_ii) − χ_e χ_i v_ei²
    //   = 1 − χ_e v_ee − χ_i v_ii + χ_e χ_i (v_ee v_ii − v_ei²)   [Eq. (dielectric-function)]
    const std::complex<T> one = unity<T>();
    const std::complex<T> block_e = one - chi_e * v.v_ee;
    const std::complex<T> block_i = one - chi_i * v.v_ii;
    return block_e * block_i - chi_e * chi_i * v.v_ei * v.v_ei;
}

template<ScalarPhysical T>
[[nodiscard]] RpaResult<T> rpa_from_susceptibilities(std::complex<T> chi_e,
                                                       std::complex<T> chi_i,
                                                       BarePotentials<T> v)
{
    const std::complex<T> epsilon = dielectric_from_susceptibilities(chi_e, chi_i, v);

    RpaResult<T> result;
    result.chi_ee = (chi_e - chi_e * v.v_ii * chi_i) / epsilon;
    result.chi_ii = (chi_i - chi_i * v.v_ee * chi_e) / epsilon;
    result.chi_ei = (chi_e * v.v_ei * chi_i) / epsilon;
    return result;
}

}  // namespace

ResponsePathway select_response_pathway(CouplingRegime<double> regime,
                                        ResponsePathway requested,
                                        bool force_pathway,
                                        bool auto_bypass) noexcept
{
    const bool strong = is_strong_coupling(regime);

    switch (requested) {
    case ResponsePathway::StandardRPA:
        if (auto_bypass && strong) {
            // Auto-bypass selects production ZetaRPA only — never Experimental.
            return ResponsePathway::ZetaRPA;
        }
        return ResponsePathway::StandardRPA;

    case ResponsePathway::PolyLogRPA:
        // Scalar Pathway A: keep the requested pathway so the evaluator can apply
        // s = f (including the weak-coupling RPA identity gate). Never auto-promote.
        return ResponsePathway::PolyLogRPA;

    case ResponsePathway::ZetaRPA:
    case ResponsePathway::ZetaRPA_Experimental:
        if (strong || force_pathway) {
            return requested;
        }
        return ResponsePathway::StandardRPA;
    }

    return ResponsePathway::StandardRPA;
}

template<ScalarPhysical T>
std::complex<T> evaluate_dielectric(std::complex<T> chi_e,
                                    std::complex<T> chi_i,
                                    BarePotentials<T> potentials)
{
    return dielectric_from_susceptibilities(chi_e, chi_i, potentials);
}

template<ScalarPhysical T>
std::complex<T> evaluate_dielectric(LindhardResult<T> chi_e,
                                    LindhardResult<T> chi_i,
                                    BarePotentials<T> potentials)
{
    return evaluate_dielectric(as_complex(chi_e), as_complex(chi_i), potentials);
}

template<ScalarPhysical T>
RpaResult<T> evaluate_rpa_susceptibility(std::complex<T> chi_e,
                                         std::complex<T> chi_i,
                                         BarePotentials<T> potentials)
{
    return rpa_from_susceptibilities(chi_e, chi_i, potentials);
}

template<ScalarPhysical T>
RpaResult<T> evaluate_rpa_susceptibility(LindhardResult<T> chi_e,
                                         LindhardResult<T> chi_i,
                                         BarePotentials<T> potentials)
{
    return evaluate_rpa_susceptibility(as_complex(chi_e), as_complex(chi_i), potentials);
}

template<ScalarPhysical T>
RpaResult<T> evaluate_rpa_susceptibility(LindhardResult<T> chi_e,
                                         LindhardResult<T> chi_i,
                                         BarePotentials<T> potentials,
                                         T electron_density_of_states,
                                         T ion_density_of_states)
{
    return evaluate_rpa_susceptibility(
        susceptibility_in_natural_units(chi_e, electron_density_of_states),
        susceptibility_in_natural_units(chi_i, ion_density_of_states),
        potentials);
}

template std::complex<double> evaluate_dielectric<double>(std::complex<double>,
                                                        std::complex<double>,
                                                        BarePotentials<double>);
template std::complex<double> evaluate_dielectric<double>(LindhardResult<double>,
                                                        LindhardResult<double>,
                                                        BarePotentials<double>);
template RpaResult<double> evaluate_rpa_susceptibility<double>(std::complex<double>,
                                                               std::complex<double>,
                                                               BarePotentials<double>);
template RpaResult<double> evaluate_rpa_susceptibility<double>(LindhardResult<double>,
                                                               LindhardResult<double>,
                                                               BarePotentials<double>);
template RpaResult<double> evaluate_rpa_susceptibility<double>(LindhardResult<double>,
                                                               LindhardResult<double>,
                                                               BarePotentials<double>,
                                                               double,
                                                               double);

}  // namespace mosaiq
