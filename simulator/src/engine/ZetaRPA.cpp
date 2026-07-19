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

std::optional<double> evaluate_coupling_shape_factor(CouplingRegime<double> regime,
                                                      ZetaWeightParameters params)
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

std::optional<double> evaluate_zeta_weight(ResponsePathway pathway,
                                           CouplingRegime<double> regime,
                                           BorweinPolicy borwein,
                                           ZetaWeightParameters params)
{
    switch (pathway) {
    case ResponsePathway::StandardRPA:
        return 1.0;

    case ResponsePathway::PolyLogRPA:
        // Pathway A does not use W_ζ; refuse rather than invent a dress.
        return std::nullopt;

    case ResponsePathway::ZetaRPA: {
        const auto f = evaluate_coupling_shape_factor(regime, params);
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

std::optional<ZetaRpaMatrixResult<double>> evaluate_zeta_rpa_matrix(
    const ZetaRpaMatrixInputs<double>& inputs)
{
    if (inputs.pathway == ResponsePathway::StandardRPA) {
        return std::nullopt;
    }

    const std::complex<double> chi_e = as_complex(inputs.chi_lindhard_e);
    const std::complex<double> chi_i = as_complex(inputs.chi_lindhard_i);
    const BarePotentials<double>& v = inputs.potentials;

    if (!std::isfinite(v.v_ee) || !std::isfinite(v.v_ii) || !std::isfinite(v.v_ei) ||
        !std::isfinite(chi_e.real()) || !std::isfinite(chi_e.imag()) ||
        !std::isfinite(chi_i.real()) || !std::isfinite(chi_i.imag())) {
        return std::nullopt;
    }

    // Cross channel: arithmetic mean of species regimes (species-asymmetry contract).
    const CouplingRegime<double> regime_ei{
        .rs = 0.5 * (inputs.regime_e.rs + inputs.regime_i.rs),
        .gamma_plasma = 0.5 * (inputs.regime_e.gamma_plasma + inputs.regime_i.gamma_plasma),
        .tau = 0.5 * (inputs.regime_e.tau + inputs.regime_i.tau),
    };

    if (!std::isfinite(regime_ei.rs) || regime_ei.rs <= 0.0 ||
        !std::isfinite(regime_ei.gamma_plasma) || regime_ei.gamma_plasma < 0.0 ||
        !std::isfinite(regime_ei.tau) || regime_ei.tau < 0.0) {
        return std::nullopt;
    }

    const ResponsePathway selected_e = select_response_pathway(
        inputs.regime_e, inputs.pathway, inputs.force_pathway, /*auto_bypass=*/false);
    const ResponsePathway selected_i = select_response_pathway(
        inputs.regime_i, inputs.pathway, inputs.force_pathway, /*auto_bypass=*/false);
    const ResponsePathway selected_ei = select_response_pathway(
        regime_ei, inputs.pathway, inputs.force_pathway, /*auto_bypass=*/false);

    // Refuse silent full-matrix RPA fallback inside the Zeta pathway.
    if (selected_e == ResponsePathway::StandardRPA &&
        selected_i == ResponsePathway::StandardRPA &&
        selected_ei == ResponsePathway::StandardRPA) {
        return std::nullopt;
    }

    const auto w_ee =
        evaluate_zeta_weight(selected_e, inputs.regime_e, inputs.borwein, inputs.weight_params);
    const auto w_ii =
        evaluate_zeta_weight(selected_i, inputs.regime_i, inputs.borwein, inputs.weight_params);
    const auto w_ei =
        evaluate_zeta_weight(selected_ei, regime_ei, inputs.borwein, inputs.weight_params);

    if (!w_ee || !w_ii || !w_ei) {
        return std::nullopt;
    }

    // Dress bare potentials; algebra matches RPA.cpp with v → W v.
    const BarePotentials<double> v_eff{
        .v_ee = v.v_ee * (*w_ee),
        .v_ii = v.v_ii * (*w_ii),
        .v_ei = v.v_ei * (*w_ei),
    };

    const std::complex<double> one{1.0, 0.0};
    const std::complex<double> block_e = one - chi_e * v_eff.v_ee;
    const std::complex<double> block_i = one - chi_i * v_eff.v_ii;
    const std::complex<double> epsilon =
        block_e * block_i - chi_e * chi_i * v_eff.v_ei * v_eff.v_ei;

    if (!std::isfinite(epsilon.real()) || !std::isfinite(epsilon.imag()) ||
        (epsilon.real() == 0.0 && epsilon.imag() == 0.0)) {
        return std::nullopt;
    }

    const std::complex<double> chi_ee = (chi_e - chi_e * v_eff.v_ii * chi_i) / epsilon;
    const std::complex<double> chi_ii = (chi_i - chi_i * v_eff.v_ee * chi_e) / epsilon;
    const std::complex<double> chi_ei = (chi_e * v_eff.v_ei * chi_i) / epsilon;

    if (!std::isfinite(chi_ee.real()) || !std::isfinite(chi_ee.imag()) ||
        !std::isfinite(chi_ii.real()) || !std::isfinite(chi_ii.imag()) ||
        !std::isfinite(chi_ei.real()) || !std::isfinite(chi_ei.imag())) {
        return std::nullopt;
    }

    return ZetaRpaMatrixResult<double>{
        .chi_ee = chi_ee,
        .chi_ii = chi_ii,
        .chi_ei = chi_ei,
        .epsilon = epsilon,
        .zeta_weight_ee = *w_ee,
        .zeta_weight_ii = *w_ii,
        .zeta_weight_ei = *w_ei,
        .pathway = inputs.pathway,
    };
}

}  // namespace mosaiq
