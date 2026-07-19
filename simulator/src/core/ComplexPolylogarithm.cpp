/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "core/ComplexPolylogarithm.hpp"

#include "core/RiemannZetaBorwein.hpp"

#include <cmath>

namespace mosaiq {

namespace {

[[nodiscard]] bool is_finite_complex(std::complex<double> z) noexcept
{
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}

/// Direct series Li_s(z) = Σ_{n=1}^N z^n / n^s  (absolute convergence for |z| < 1).
[[nodiscard]] std::optional<std::complex<double>>
polylog_series(double s, std::complex<double> z, const PolylogPolicy& policy)
{
    if (z == std::complex<double>{0.0, 0.0}) {
        return std::complex<double>{0.0, 0.0};
    }

    std::complex<double> sum{0.0, 0.0};
    std::complex<double> z_pow{1.0, 0.0};  // will multiply by z each step → z^n

    for (std::size_t n = 1; n <= policy.max_series_terms; ++n) {
        z_pow *= z;
        const double n_s = std::pow(static_cast<double>(n), s);
        if (!(n_s > 0.0) || !std::isfinite(n_s)) {
            return std::nullopt;
        }
        const std::complex<double> term = z_pow / n_s;
        if (!is_finite_complex(term)) {
            return std::nullopt;
        }
        sum += term;
        if (std::abs(term) < policy.tolerance * (1.0 + std::abs(sum))) {
            return sum;
        }
    }
    return std::nullopt;  // failed to meet tolerance within max_series_terms
}

}  // namespace

std::optional<std::complex<double>> evaluate_complex_polylog(double s,
                                                             std::complex<double> z,
                                                             PolylogPolicy policy)
{
    if (!std::isfinite(s) || s <= 1.0 || !is_finite_complex(z)) {
        return std::nullopt;
    }
    if (!(policy.series_radius > 0.0) || !(policy.tolerance > 0.0) ||
        policy.max_series_terms == 0) {
        return std::nullopt;
    }

    const double abs_z = std::abs(z);

    // Catastrophe boundary on the real axis: Li_s(1) = ζ(s) for Re(s) > 1.
    if (std::abs(z - 1.0) < 1.0e-14) {
        const auto zeta = riemann_zeta_borwein(s);
        if (!zeta) {
            return std::nullopt;
        }
        return std::complex<double>{*zeta, 0.0};
    }

    // Fast path — interior of the disk of absolute convergence (shrunk by policy).
    if (abs_z < policy.series_radius) {
        return polylog_series(s, z, policy);
    }

    // -------------------------------------------------------------------------
    // TODO(Pathway-A / PolyLog-RPA heavy path):
    // Analytic continuation for |z| >= series_radius via the Mellin-type integral
    //
    //   Li_s(z) = z / Γ(s) ∫_0^∞ t^{s-1} / (e^t − z) dt   (appropriate branch of z)
    //
    // to be evaluated with the existing Sinc-Quadrature / deterministic quadrature
    // stack. Until that contour/integral path is pinned, refuse to invent a value.
    // -------------------------------------------------------------------------
    (void)abs_z;
    return std::nullopt;
}

}  // namespace mosaiq
