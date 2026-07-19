/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "core/RiemannZetaBorwein.hpp"

#include <cmath>
#include <vector>

namespace mosaiq {

namespace {

/// Borwein coefficients d_k for truncation order n (P. Borwein, CMAA 2000).
/// Uses the stable ratio a_i / a_{i-1} to avoid factorial overflow.
///   a_0 = 1/n,
///   a_i / a_{i-1} = 4(n+i-1) / ( (n-i+1) (2i) (2i-1) ),
///   d_k = n Σ_{i=0}^{min(k,n)} a_i.
[[nodiscard]] std::vector<double> borwein_coefficients(std::size_t n)
{
    std::vector<double> d(n + 1, 0.0);
    double a = 1.0 / static_cast<double>(n);
    double partial = 0.0;

    for (std::size_t k = 0; k <= n; ++k) {
        if (k >= 1) {
            // a_k / a_{k-1} = 4 (n+k-1)(n-k+1) / ( (2k)(2k-1) )
            const double num = 4.0 * static_cast<double>(n + k - 1) *
                               static_cast<double>(n - k + 1);
            const double den = static_cast<double>(2 * k) *
                               static_cast<double>(2 * k - 1);
            a *= num / den;
        }
        partial += a;
        d[k] = static_cast<double>(n) * partial;
    }
    return d;
}

}  // namespace

std::optional<double> riemann_zeta_borwein(double s, BorweinPolicy policy)
{
    if (!std::isfinite(s) || s <= 1.0 || policy.truncation_order == 0) {
        return std::nullopt;
    }

    const std::size_t n = policy.truncation_order;
    const std::vector<double> d = borwein_coefficients(n);
    const double d_n = d[n];
    if (!std::isfinite(d_n) || d_n == 0.0) {
        return std::nullopt;
    }

    const double two_pow = std::pow(2.0, 1.0 - s);
    const double denom = d_n * (1.0 - two_pow);
    if (!std::isfinite(denom) || denom == 0.0) {
        return std::nullopt;
    }

    double series = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const double sign = ((k % 2) == 0) ? 1.0 : -1.0;
        // Sign convention: (d_n − d_k) yields ζ(s) > 0 for s > 1 (equivalent to
        // flipping the classic (d_k − d_n) form under the same alternating prefactor).
        series += sign * (d_n - d[k]) / std::pow(static_cast<double>(k + 1), s);
    }

    const double zeta = series / denom;
    if (!std::isfinite(zeta)) {
        return std::nullopt;
    }
    return zeta;
}

}  // namespace mosaiq
