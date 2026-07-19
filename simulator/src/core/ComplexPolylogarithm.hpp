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

#include <complex>
#include <cstddef>
#include <optional>

namespace mosaiq {

/// Policy for deterministic complex polylogarithm Li_s(z) (Pathway A / PolyLog-RPA).
/// Fast path: power series inside |z| < series_radius.
/// Heavy path: Mellin integral on the real line after t = e^u, summed by a
/// truncated sinc / trapezoidal rule on u ∈ [-N h, N h].
struct PolylogPolicy {
    std::size_t max_series_terms{1000};
    double tolerance{1.0e-12};
    /// Use the Taylor series only for |z| strictly below this radius.
    double series_radius{0.75};
    /// Half-range truncation for the infinite-line sinc sum (k = -N … N).
    /// Defaults target ≈ double precision on the double-exponentially decaying integrand.
    std::size_t sinc_N{128};
    /// Uniform step h in the u-variable (t = e^u).
    double sinc_h{0.1};
    /// Relative floor for |e^{e^u} − z|; below this the real-axis pole is too close.
    double pole_guard{1.0e-14};
};

/// Evaluate Li_s(z) for real s > 0 and complex z.
/// Returns nullopt if s ≤ 0, inputs are non-finite, the series fails to converge,
/// z = 1 with s ≤ 1 (ζ pole / harmonic divergence), the integrand encounters a
/// near-pole on the positive real axis (branch cut), or the heavy-path quadrature
/// fails finiteness checks. The z → 1 identity Li_s(1) = ζ(s) requires s > 1.
[[nodiscard]] std::optional<std::complex<double>>
evaluate_complex_polylog(double s,
                         std::complex<double> z,
                         PolylogPolicy policy = {});

/// Convenience overload for ScalarPhysical s with complex z in the same float type.
template<ScalarPhysical T>
[[nodiscard]] std::optional<std::complex<T>>
evaluate_complex_polylog(T s, std::complex<T> z, PolylogPolicy policy = {})
{
    const auto out = evaluate_complex_polylog(static_cast<double>(s),
                                              std::complex<double>{static_cast<double>(z.real()),
                                                                   static_cast<double>(z.imag())},
                                              policy);
    if (!out) {
        return std::nullopt;
    }
    return std::complex<T>{static_cast<T>(out->real()), static_cast<T>(out->imag())};
}

}  // namespace mosaiq
