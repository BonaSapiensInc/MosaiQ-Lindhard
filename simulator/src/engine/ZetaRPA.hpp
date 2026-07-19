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
#include "core/RiemannZetaBorwein.hpp"
#include "engine/RPA.hpp"
#include "physics/Constants.hpp"

#include <complex>
#include <optional>

namespace mosaiq {

/// Theory coefficients for locked W_ζ (manuscript Appendix C).
struct ZetaWeightParameters {
    double alpha{constants::zeta_weight_alpha};
    double beta{constants::zeta_weight_beta};
    double gamma{constants::zeta_weight_gamma};
    double delta{constants::zeta_weight_delta};
};

/// Scalar Zeta-RPA result (Phase Z1 — single component only).
template<ScalarPhysical T = double>
struct ZetaRpaResult {
    std::complex<T> chi{};
    std::complex<T> epsilon{};
    T zeta_weight{};
    ResponsePathway pathway{ResponsePathway::ZetaRPA};
};

/// Inputs for scalar Zeta-RPA. χ^L must already be in DOS-restored (natural) units.
template<ScalarPhysical T = double>
struct ZetaRpaInputs {
    WaveVector<T> q{};
    Frequency<T> omega{};
    LindhardResult<T> chi_lindhard{};
    T bare_potential{};
    CouplingRegime<T> regime{};
    BorweinPolicy borwein{};
    ResponsePathway pathway{ResponsePathway::ZetaRPA};
    bool force_pathway{false};
    ZetaWeightParameters weight_params{};
};

/// Locked production weight:
///   f = α Γ^β / (1 + γ r_s^{-δ} τ),
///   W = ζ(1+f)/ζ(1)  ≅  f·ζ(1+f)  (Laurent regularization; W=1 at f=0).
/// Experimental pathway retains a provisional A/B dress.
[[nodiscard]] std::optional<double> evaluate_zeta_weight(
    ResponsePathway pathway,
    CouplingRegime<double> regime,
    BorweinPolicy borwein = {},
    ZetaWeightParameters params = {});

/// Pure scalar Zeta-RPA evaluation. Returns nullopt on non-finite inputs or failed ζ.
/// Does not renegotiate causality: χ^L is an immutable upstream Lindhard input.
[[nodiscard]] std::optional<ZetaRpaResult<double>> evaluate_zeta_rpa(
    const ZetaRpaInputs<double>& inputs);

}  // namespace mosaiq
