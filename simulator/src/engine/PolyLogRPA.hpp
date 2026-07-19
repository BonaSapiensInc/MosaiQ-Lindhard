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

#include "core/ComplexPolylogarithm.hpp"
#include "core/Concepts.hpp"
#include "engine/RPA.hpp"
#include "engine/ZetaRPA.hpp"

#include <complex>
#include <optional>

namespace mosaiq {

/// Scalar PolyLog-RPA result (Pathway A — formal diagrammatic definition).
/// χ^(s) = χ^L Li_s(v χ^L),  ε^(s) = 1 / Li_s(v χ^L)  (manuscript Eq. polylog-rpa).
template<ScalarPhysical T = double>
struct PolyLogRpaResult {
    std::complex<T> chi{};
    std::complex<T> epsilon{};
    T s_parameter{};  ///< Topological penalty index s (= thermodynamic distance f)
    ResponsePathway pathway{ResponsePathway::PolyLogRPA};
};

/// Inputs for scalar PolyLog-RPA. χ^L must already be in DOS-restored (natural) units.
template<ScalarPhysical T = double>
struct PolyLogRpaInputs {
    WaveVector<T> q{};
    Frequency<T> omega{};
    LindhardResult<T> chi_lindhard{};
    T bare_potential{};
    CouplingRegime<T> regime{};
    PolylogPolicy polylog_policy{};
    ZetaWeightParameters weight_params{};
};

/// Pure scalar PolyLog-RPA evaluation (Pathway A).
/// Sets s = f from the same coupling shape factor as Zeta-RPA.
/// For s → 0 (weak coupling) returns ordinary scalar RPA to avoid Li_s precision loss.
/// Returns nullopt on non-finite inputs or failed polylog evaluation.
/// Does not extend to the coupled two-component RPA tensor.
[[nodiscard]] std::optional<PolyLogRpaResult<double>> evaluate_polylog_rpa(
    const PolyLogRpaInputs<double>& inputs);

}  // namespace mosaiq
