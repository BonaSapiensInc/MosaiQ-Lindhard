#pragma once

#include "core/Concepts.hpp"

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <optional>

namespace mosaiq {

/// Bracketed scalar objective f(x) used by Brent's method (no derivatives).
template<typename F>
concept BracketedRootFunction = std::invocable<F, double> && requires(F f, double x) {
    { f(x) } -> std::convertible_to<double>;
};

template<ScalarPhysical T = double>
struct BrentRootFinderPolicy {
    std::size_t iter_max{128};
    T tolerance{1.0e-12};
};

/// Brent's method: bracketed root of f(x) = 0 on [a, b] with f(a) f(b) < 0.
/// Pure, deterministic, and derivative-free.
template<ScalarPhysical T, BracketedRootFunction F>
[[nodiscard]] std::optional<T> brent_root(F f,
                                          T a,
                                          T b,
                                          BrentRootFinderPolicy<T> policy = {}) noexcept
{
    T fa = f(a);
    T fb = f(b);

    if (fa == T{0}) {
        return a;
    }
    if (fb == T{0}) {
        return b;
    }
    if (fa * fb > T{0}) {
        return std::nullopt;
    }

    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    T c = a;
    T fc = fa;
    T d = b - a;
    T e = d;

    for (std::size_t iter = 0; iter < policy.iter_max; ++iter) {
        if (fb == T{0} || std::abs(b - a) <= policy.tolerance) {
            return b;
        }

        T s{};
        const bool distinct_ab = std::abs(fa - fb) > std::numeric_limits<T>::epsilon();
        const bool distinct_cb = std::abs(fb - fc) > std::numeric_limits<T>::epsilon();

        if (distinct_ab && distinct_cb) {
            // Inverse quadratic interpolation.
            s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb));
        } else if (distinct_ab) {
            // Secant step.
            s = b - fb * (b - a) / (fb - fa);
        } else {
            // Bisection fallback.
            s = (a + b) / T{2};
        }

        const T mid = (a + b) / T{2};
        const bool must_bisect =
            (s < mid && s < b) || (s > mid && s > b) ||
            (e < policy.tolerance && std::abs(b - s) >= policy.tolerance / T{2}) ||
            (e >= policy.tolerance && std::abs(b - s) >= std::abs(e) / T{2});

        if (must_bisect) {
            d = mid - b;
            e = d;
        } else {
            e = d;
            d = s - b;
        }

        c = b;
        fc = fb;

        if (std::abs(d) > policy.tolerance) {
            b += d;
        } else {
            b += (b > a ? policy.tolerance : -policy.tolerance);
        }

        fb = f(b);

        if (fa * fb > T{0}) {
            a = c;
            fa = fc;
            d = b - a;
            e = d;
        }
    }

    if (std::abs(fb) <= policy.tolerance * T{10}) {
        return b;
    }

    return std::nullopt;
}

}  // namespace mosaiq
