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

/// Multi-component Zeta-RPA result (Phase Z3 — {ee, ii, ei} channel layout).
template<ScalarPhysical T = double>
struct ZetaRpaMatrixResult {
    std::complex<T> chi_ee{};
    std::complex<T> chi_ii{};
    std::complex<T> chi_ei{};
    std::complex<T> epsilon{};  ///< two-component dielectric determinant ε^ζ
    T zeta_weight_ee{};
    T zeta_weight_ii{};
    T zeta_weight_ei{};
    ResponsePathway pathway{ResponsePathway::ZetaRPA};
};

/// Inputs for two-component Zeta-RPA. χ^L_e/i must already be DOS-restored.
template<ScalarPhysical T = double>
struct ZetaRpaMatrixInputs {
    WaveVector<T> q{};
    Frequency<T> omega{};
    LindhardResult<T> chi_lindhard_e{};
    LindhardResult<T> chi_lindhard_i{};
    BarePotentials<T> potentials{};
    CouplingRegime<T> regime_e{};
    CouplingRegime<T> regime_i{};
    BorweinPolicy borwein{};
    ResponsePathway pathway{ResponsePathway::ZetaRPA};
    bool force_pathway{false};
    ZetaWeightParameters weight_params{};
};

/// Thermodynamic distance for locked W_ζ / PolyLog s-parameter:
///   f = α Γ^β / (1 + γ r_s^{-δ} τ).
[[nodiscard]] std::optional<double> evaluate_coupling_shape_factor(
    CouplingRegime<double> regime,
    ZetaWeightParameters params = {});

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

/// Two-component Zeta-RPA matrix inversion with per-channel W_ζ dressing.
/// Algebraically continuous with `evaluate_rpa_susceptibility` when all W_ζ → 1.
/// Does not rewrite manuscript StandardRPA pipelines.
[[nodiscard]] std::optional<ZetaRpaMatrixResult<double>> evaluate_zeta_rpa_matrix(
    const ZetaRpaMatrixInputs<double>& inputs);

}  // namespace mosaiq
