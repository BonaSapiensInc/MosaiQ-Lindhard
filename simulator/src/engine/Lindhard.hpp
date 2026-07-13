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

namespace mosaiq {

/// Grand-canonical imaginary part — Giuliani & Vignale Eq. (4.45) / manuscript Eq. (38).
template<ScalarPhysical T = double>
[[nodiscard]] T imaginary_lindhard(WaveVector<T> q,
                                   Frequency<T> omega,
                                   ReducedTemperature<T> tau,
                                   ReducedChemicalPotential<T> gamma);

/// Canonical real part via Kramers–Krönig / Hilbert sinc-quadrature (legacy `get_KKrealLindhard`).
/// For τ < singularity_excision_tau_threshold, uses reverse-Dedekind singularity excision:
/// analytic Stratonovich step Hilbert transform + smooth softplus residual only in the sinc sum.
template<ScalarPhysical T = double>
[[nodiscard]] T real_lindhard_kk(WaveVector<T> q,
                                 Frequency<T> omega,
                                 ReducedTemperature<T> tau,
                                 ReducedChemicalPotential<T> gamma);

/// Unified χ^L real and imaginary parts in rational units (Re/Im before DOS scaling).
template<ScalarPhysical T = double>
[[nodiscard]] LindhardResult<T> evaluate_lindhard(WaveVector<T> q,
                                                  Frequency<T> omega,
                                                  ReducedTemperature<T> tau,
                                                  ReducedChemicalPotential<T> gamma);

}  // namespace mosaiq
