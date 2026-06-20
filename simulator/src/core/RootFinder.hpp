#pragma once

#include "core/Concepts.hpp"

#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <optional>

namespace sophus {

/// Value f(x) and derivative f'(x) returned together for Newton steps.
template<ScalarPhysical T>
struct NewtonEval {
    T value{};
    T derivative{};
};

template<ScalarPhysical T = double>
struct RootFinderPolicy {
    std::size_t iter_max{100};
    T epsilon{1.0e-14};
};

/// Solve f(x) = target using Newton–Raphson: x ← x + (target − f(x)) / f'(x).
/// Returns the final iterate (legacy semantics: usable even if ε is not reached within iter_max).
template<ScalarPhysical T, typename Eval>
    requires requires(Eval eval, T x) {
        { eval(x) } -> std::convertible_to<NewtonEval<T>>;
    }
[[nodiscard]] std::optional<T> newton_raphson(Eval eval,
                                              T initial_guess,
                                              T target,
                                              RootFinderPolicy<T> policy = {})
{
    T x = initial_guess;

    for (std::size_t iter = 0; iter < policy.iter_max; ++iter) {
        const NewtonEval<T> probe = eval(x);
        if (probe.derivative == T{0}) {
            return std::nullopt;
        }

        const T delta = (target - probe.value) / probe.derivative;
        x += delta;

        if (std::abs(delta) <= policy.epsilon) {
            break;
        }
    }

    return x;
}

}  // namespace sophus
