#pragma once

#include "core/Concepts.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <functional>
#include <ranges>
#include <type_traits>
#include <vector>

namespace mosaiq {

/// Parameters for the modified sinc rule (Davis–Rabinowitz Eq. 3.4.6.21; legacy N=256, h = sqrt(7 pi / N)).
template<ScalarPhysical T = double>
struct SincQuadratureParams {
    std::size_t half_nodes{constants::default_gv_sinc_nodes};
    T strip_factor{constants::sinc_strip_factor};

    [[nodiscard]] constexpr T step() const noexcept
    {
        return static_cast<T>(std::sqrt(strip_factor * constants::pi /
                                        static_cast<T>(half_nodes)));
    }

    [[nodiscard]] constexpr std::size_t node_count() const noexcept
    {
        return 2 * half_nodes + 1;
    }
};

/// Arguments for the GV Eq. (4.44) grand-canonical real Lindhard sinc integrand.
template<ScalarPhysical T = double>
struct GvLindhardSincArgs {
    WaveVector<T> q;
    Frequency<T> omega;
    ReducedTemperature<T> tau;
    ReducedChemicalPotential<T> gamma;
};

/// One sinc abscissa sample: index k, mapped x = asinh(e^{kh}), and quadrature weight factor.
template<ScalarPhysical T = double>
struct SincSample {
    long k{};
    T kh{};
    T x{};
    T weight{};  // 1 / sqrt(1 + exp(-2 kh)) in legacy get_GVrealLindhard
};

/// Callable integrand at a single sinc node (pure; no mutation).
template<ScalarPhysical T = double>
using GvLindhardIntegrand = std::function<T(const SincSample<T>&, const GvLindhardSincArgs<T>&)>;

/// Declarations only — implementations deferred to SincQuadrature.cpp (Phase 1).
template<ScalarPhysical T = double>
class SincQuadrature {
public:
    explicit SincQuadrature(SincQuadratureParams<T> params = {});

    [[nodiscard]] SincQuadratureParams<T> params() const noexcept;

    /// Lazy view of sinc indices k ∈ [-N, N] (legacy for-loop body).
    [[nodiscard]] auto node_indices() const;

    /// Map index k → SincSample (asinh map + legacy weight).
    [[nodiscard]] SincSample<T> sample_at(long k) const;

    /// Fold integrand contributions over all nodes via `<ranges>` (Phase 1 implementation).
    [[nodiscard]] T integrate(GvLindhardIntegrand<T> integrand,
                            const GvLindhardSincArgs<T>& args) const;

    /// GV grand-canonical real Lindhard: legacy `get_GVrealLindhard` pipeline entry point.
    [[nodiscard]] T evaluate_gv_real_lindhard(const GvLindhardSincArgs<T>& args) const;

    /// Collect all samples materialized (testing / debugging only).
    [[nodiscard]] std::vector<SincSample<T>> all_samples() const;

private:
    SincQuadratureParams<T> params_;
};

// --- free-function pipeline hooks (declaration only) ---

template<ScalarPhysical T>
[[nodiscard]] T asinh_sinc_abscissa(T kh) noexcept;

template<ScalarPhysical T>
[[nodiscard]] T sinc_weight_factor(T kh) noexcept;

template<ScalarPhysical T>
[[nodiscard]] GvLindhardIntegrand<T> make_gv_lindhard_integrand();

/// Cauchy principal-value kernel [1 - cos(π Δ/h)] / (π Δ/h) with L'Hôpital regularization at grid nodes.
template<ScalarPhysical T>
[[nodiscard]] T cauchy_principal_value_kernel(T diff, T pi_over_h) noexcept
{
    constexpr T resonance_tolerance = T{1.0e-13};
    if (std::abs(diff) <= resonance_tolerance) {
        // lim_{Δ→0} [1 - cos(cΔ)] / (cΔ) = 0 (L'Hôpital); avoids IEEE 754 0/0 → NaN traps.
        return T{0};
    }
    const T argument = pi_over_h * diff;
    return (T{1} - std::cos(argument)) / argument;
}

}  // namespace mosaiq
