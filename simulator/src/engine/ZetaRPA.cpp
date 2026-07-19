/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "engine/ZetaRPA.hpp"

#include "physics/Constants.hpp"

#include <cmath>

namespace mosaiq {

std::optional<double> evaluate_zeta_weight(ResponsePathway pathway,
                                           CouplingRegime<double> regime,
                                           BorweinPolicy borwein)
{
    switch (pathway) {
    case ResponsePathway::StandardRPA:
    case ResponsePathway::ZetaRPA:
        // Production identity until manuscript locks W_ζ.
        return 1.0;

    case ResponsePathway::ZetaRPA_Experimental: {
        if (!std::isfinite(regime.gamma_plasma) || regime.gamma_plasma < 0.0) {
            return std::nullopt;
        }
        const auto zeta2 = riemann_zeta_borwein(2.0, borwein);
        const auto zeta3 = riemann_zeta_borwein(3.0, borwein);
        if (!zeta2 || !zeta3 || *zeta2 == 0.0) {
            return std::nullopt;
        }
        const double gamma = regime.gamma_plasma;
        const double dress = constants::zeta_rpa_experimental_kappa * (gamma / (1.0 + gamma)) *
                             (*zeta3 / *zeta2);
        const double weight = 1.0 + dress;
        if (!std::isfinite(weight)) {
            return std::nullopt;
        }
        return weight;
    }
    }

    return std::nullopt;
}

std::optional<ZetaRpaResult<double>> evaluate_zeta_rpa(const ZetaRpaInputs<double>& inputs)
{
    if (inputs.pathway == ResponsePathway::StandardRPA) {
        return std::nullopt;
    }

    const ResponsePathway selected = select_response_pathway(
        inputs.regime, inputs.pathway, inputs.force_pathway, /*auto_bypass=*/false);

    if (selected == ResponsePathway::StandardRPA) {
        // Explicit policy: do not silently evaluate as RPA inside Zeta-RPA.
        return std::nullopt;
    }

    const double v = inputs.bare_potential;
    const std::complex<double> chi_l = as_complex(inputs.chi_lindhard);

    if (!std::isfinite(v) || !std::isfinite(chi_l.real()) || !std::isfinite(chi_l.imag())) {
        return std::nullopt;
    }

    const auto weight = evaluate_zeta_weight(selected, inputs.regime, inputs.borwein);
    if (!weight) {
        return std::nullopt;
    }

    const std::complex<double> one{1.0, 0.0};
    const std::complex<double> epsilon = one - v * (*weight) * chi_l;
    if (!std::isfinite(epsilon.real()) || !std::isfinite(epsilon.imag())) {
        return std::nullopt;
    }
    if (epsilon.real() == 0.0 && epsilon.imag() == 0.0) {
        return std::nullopt;
    }

    const std::complex<double> chi = chi_l / epsilon;
    if (!std::isfinite(chi.real()) || !std::isfinite(chi.imag())) {
        return std::nullopt;
    }

    return ZetaRpaResult<double>{
        .chi = chi,
        .epsilon = epsilon,
        .zeta_weight = *weight,
        .pathway = selected,
    };
}

}  // namespace mosaiq
