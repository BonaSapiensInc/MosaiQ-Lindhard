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

namespace {

/// f = α Γ^β / (1 + γ r_s^{-δ} τ)
[[nodiscard]] std::optional<double> coupling_shape_factor(CouplingRegime<double> regime,
                                                         ZetaWeightParameters params) noexcept
{
    if (!std::isfinite(regime.gamma_plasma) || regime.gamma_plasma < 0.0 ||
        !std::isfinite(regime.rs) || regime.rs <= 0.0 || !std::isfinite(regime.tau) ||
        regime.tau < 0.0) {
        return std::nullopt;
    }
    if (!std::isfinite(params.alpha) || !std::isfinite(params.beta) ||
        !std::isfinite(params.gamma) || !std::isfinite(params.delta) || params.alpha < 0.0) {
        return std::nullopt;
    }

    const double gamma_pow = std::pow(regime.gamma_plasma, params.beta);
    const double rs_pow = std::pow(regime.rs, -params.delta);
    const double denom = 1.0 + params.gamma * rs_pow * regime.tau;
    if (!std::isfinite(gamma_pow) || !std::isfinite(rs_pow) || !std::isfinite(denom) ||
        denom <= 0.0) {
        return std::nullopt;
    }

    const double f = params.alpha * gamma_pow / denom;
    if (!std::isfinite(f) || f < 0.0) {
        return std::nullopt;
    }
    return f;
}

/// Laurent-regularized ζ(1+f)/ζ(1) ≡ f·ζ(1+f) for f>0; identity 1 at f=0.
[[nodiscard]] std::optional<double> locked_zeta_weight(double f, BorweinPolicy borwein)
{
    if (f <= constants::zeta_weight_f_floor) {
        return 1.0;
    }

    const double s = 1.0 + f;
    const auto zeta = riemann_zeta_borwein(s, borwein);
    if (!zeta) {
        return std::nullopt;
    }

    const double weight = f * (*zeta);
    if (!std::isfinite(weight) || weight <= 0.0) {
        return std::nullopt;
    }
    return weight;
}

}  // namespace

std::optional<double> evaluate_zeta_weight(ResponsePathway pathway,
                                           CouplingRegime<double> regime,
                                           BorweinPolicy borwein,
                                           ZetaWeightParameters params)
{
    switch (pathway) {
    case ResponsePathway::StandardRPA:
        return 1.0;

    case ResponsePathway::ZetaRPA: {
        const auto f = coupling_shape_factor(regime, params);
        if (!f) {
            return std::nullopt;
        }
        return locked_zeta_weight(*f, borwein);
    }

    case ResponsePathway::ZetaRPA_Experimental: {
        // Provisional A/B probe (not the locked production form).
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

    const auto weight =
        evaluate_zeta_weight(selected, inputs.regime, inputs.borwein, inputs.weight_params);
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
