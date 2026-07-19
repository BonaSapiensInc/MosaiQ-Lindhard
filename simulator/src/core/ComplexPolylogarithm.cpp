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
#include <complex>

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

/// Integrand node for Li_s after t = e^u:
///   f(u) = e^{u s} / (e^{e^u} − z)
/// with overflow-safe evaluation in the double-exponential tails.
[[nodiscard]] std::optional<std::complex<double>>
mellin_integrand_at(double s, std::complex<double> z, double u, double pole_guard)
{
    // exp(u) overflows double near u ≳ 709; the integrand is then negligible.
    if (u > 6.5) {
        return std::complex<double>{0.0, 0.0};
    }

    const double exp_u = std::exp(u);
    if (!std::isfinite(exp_u)) {
        return std::complex<double>{0.0, 0.0};
    }

    // Double-exponential regime: e^{e^u} ≫ 1.
    // Rewrite  e^{u s}/(e^{e^u}−z) = e^{u s − e^u} / (1 − z e^{−e^u}).
    constexpr double exp_u_rewrite = 40.0;
    if (exp_u >= exp_u_rewrite) {
        const double expo = u * s - exp_u;
        if (expo < -700.0) {
            return std::complex<double>{0.0, 0.0};
        }
        const std::complex<double> num{std::exp(expo), 0.0};
        const double em = std::exp(-exp_u);
        const std::complex<double> den = 1.0 - z * em;
        if (std::abs(den) < pole_guard * (1.0 + std::abs(z))) {
            return std::nullopt;
        }
        const std::complex<double> term = num / den;
        if (!is_finite_complex(term)) {
            return std::nullopt;
        }
        return term;
    }

    const double exp_exp_u = std::exp(exp_u);
    if (!std::isfinite(exp_exp_u)) {
        // Should be unreachable for exp_u < 40, but refuse rather than invent.
        return std::nullopt;
    }

    const std::complex<double> den = exp_exp_u - z;
    if (std::abs(den) < pole_guard * (1.0 + std::abs(exp_exp_u) + std::abs(z))) {
        return std::nullopt;
    }

    const double us = u * s;
    if (us < -700.0) {
        return std::complex<double>{0.0, 0.0};
    }
    if (us > 700.0) {
        return std::nullopt;
    }

    const std::complex<double> num{std::exp(us), 0.0};
    const std::complex<double> term = num / den;
    if (!is_finite_complex(term)) {
        return std::nullopt;
    }
    return term;
}

/// Heavy path: Li_s(z) = z/Γ(s) ∫_{-∞}^{∞} e^{u s}/(e^{e^u}−z) du
/// approximated by a truncated sinc / trapezoidal sum on u = k h, k ∈ [-N, N].
[[nodiscard]] std::optional<std::complex<double>>
polylog_sinc_mellin(double s, std::complex<double> z, const PolylogPolicy& policy)
{
    // Principal branch cut of Li_s lies on the real axis for z ∈ [1, ∞).
    // Exact z = 1 is handled upstream via ζ(s); refuse the rest of the cut.
    if (std::abs(z.imag()) < 1.0e-15 && z.real() >= 1.0 - 1.0e-14) {
        return std::nullopt;
    }

    if (policy.sinc_N == 0 || !(policy.sinc_h > 0.0) || !(policy.pole_guard > 0.0)) {
        return std::nullopt;
    }

    const double gamma_s = std::tgamma(s);
    if (!std::isfinite(gamma_s) || gamma_s == 0.0) {
        return std::nullopt;
    }

    std::complex<double> acc{0.0, 0.0};
    const long N = static_cast<long>(policy.sinc_N);
    for (long k = -N; k <= N; ++k) {
        const double u = static_cast<double>(k) * policy.sinc_h;
        const auto term = mellin_integrand_at(s, z, u, policy.pole_guard);
        if (!term) {
            return std::nullopt;
        }
        acc += *term;
    }

    const std::complex<double> integral = policy.sinc_h * acc;
    const std::complex<double> value = z * integral / gamma_s;
    if (!is_finite_complex(value)) {
        return std::nullopt;
    }
    return value;
}

}  // namespace

std::optional<std::complex<double>> evaluate_complex_polylog(double s,
                                                             std::complex<double> z,
                                                             PolylogPolicy policy)
{
    // Pathway A uses s = f with 0 < f ≲ O(1); series/Mellin both need s > 0.
    if (!std::isfinite(s) || s <= 0.0 || !is_finite_complex(z)) {
        return std::nullopt;
    }
    if (!(policy.series_radius > 0.0) || !(policy.tolerance > 0.0) ||
        policy.max_series_terms == 0) {
        return std::nullopt;
    }

    const double abs_z = std::abs(z);

    // Catastrophe boundary on the real axis: Li_s(1) = ζ(s) only for Re(s) > 1.
    if (std::abs(z - 1.0) < 1.0e-14) {
        if (s <= 1.0) {
            return std::nullopt;
        }
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

    // Heavy path — Mellin integral via truncated sinc sum on the u-line.
    return polylog_sinc_mellin(s, z, policy);
}

}  // namespace mosaiq
