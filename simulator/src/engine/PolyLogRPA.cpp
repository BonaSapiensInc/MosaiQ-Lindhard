/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "engine/PolyLogRPA.hpp"

#include "physics/Constants.hpp"

#include <cmath>

namespace mosaiq {

std::optional<PolyLogRpaResult<double>> evaluate_polylog_rpa(
    const PolyLogRpaInputs<double>& inputs)
{
    const double v = inputs.bare_potential;
    const std::complex<double> chi_l = as_complex(inputs.chi_lindhard);

    if (!std::isfinite(v) || !std::isfinite(chi_l.real()) || !std::isfinite(chi_l.imag())) {
        return std::nullopt;
    }

    const auto f = evaluate_coupling_shape_factor(inputs.regime, inputs.weight_params);
    if (!f) {
        return std::nullopt;
    }

    const double s = *f;
    const std::complex<double> z = v * chi_l;

    // Weak-coupling identity gate: s = f → 0 recovers ordinary scalar RPA.
    // Bypass Li_s near s = 0 (analytic Li_0(z) = z/(1−z) is numerically fragile when
    // composed as χ^L Li_0, and the production identity demanded by the CLI/tests is
    // the full RPA susceptibility χ^L/(1−z)).
    if (s <= constants::zeta_weight_f_floor) {
        const std::complex<double> chi = evaluate_scalar_rpa(chi_l, v);
        const std::complex<double> one{1.0, 0.0};
        const std::complex<double> epsilon = one - z;
        if (!std::isfinite(chi.real()) || !std::isfinite(chi.imag()) ||
            !std::isfinite(epsilon.real()) || !std::isfinite(epsilon.imag())) {
            return std::nullopt;
        }
        return PolyLogRpaResult<double>{
            .chi = chi,
            .epsilon = epsilon,
            .s_parameter = s,
            .pathway = ResponsePathway::PolyLogRPA,
        };
    }

    const auto li = evaluate_complex_polylog(s, z, inputs.polylog_policy);
    if (!li) {
        return std::nullopt;
    }
    if (li->real() == 0.0 && li->imag() == 0.0) {
        return std::nullopt;
    }

    // Manuscript Pathway A: χ^(s) = χ^L Li_s(z),  ε^(s) = χ^L / χ^(s) = 1/Li_s(z).
    const std::complex<double> chi = chi_l * (*li);
    const std::complex<double> epsilon = 1.0 / (*li);

    if (!std::isfinite(chi.real()) || !std::isfinite(chi.imag()) ||
        !std::isfinite(epsilon.real()) || !std::isfinite(epsilon.imag())) {
        return std::nullopt;
    }

    return PolyLogRpaResult<double>{
        .chi = chi,
        .epsilon = epsilon,
        .s_parameter = s,
        .pathway = ResponsePathway::PolyLogRPA,
    };
}

}  // namespace mosaiq
