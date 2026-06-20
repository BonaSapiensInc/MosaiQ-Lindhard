#pragma once

#include <concepts>
#include <cstddef>
#include <type_traits>

namespace sophus {

/// Legacy code passed bare `double` for q, ω, τ — compile-time guard against that class of bug.
template<typename T>
concept ScalarPhysical = std::floating_point<T>;

namespace detail {

template<typename Tag, ScalarPhysical T>
struct StrongScalar {
    T value{};

    constexpr StrongScalar() = default;

    explicit constexpr StrongScalar(T v) noexcept : value{v} {}

    [[nodiscard]] constexpr T raw() const noexcept { return value; }

    friend constexpr auto operator<=>(const StrongScalar&, const StrongScalar&) = default;
};

}  // namespace detail

struct WaveVectorTag {};
struct FrequencyTag {};
struct TemperatureTag {};
struct ReducedTemperatureTag {};
struct ReducedChemicalPotentialTag {};

/// Wave vector q (legacy: rational units w.r.t. k_F, or Hartree atomic units depending on layer).
template<ScalarPhysical T = double>
using WaveVector = detail::StrongScalar<WaveVectorTag, T>;

/// Frequency ω (legacy: rational units w.r.t. E_F).
template<ScalarPhysical T = double>
using Frequency = detail::StrongScalar<FrequencyTag, T>;

/// Absolute temperature T [K] (legacy: `_temperature` in TemperatureFermionComponent).
template<ScalarPhysical T = double>
using Temperature = detail::StrongScalar<TemperatureTag, T>;

/// Reduced temperature τ = k_B T / E_F (legacy: `_reducedTemperature`).
template<ScalarPhysical T = double>
using ReducedTemperature = detail::StrongScalar<ReducedTemperatureTag, T>;

/// Reduced chemical potential γ = μ / E_F (legacy: `_reducedChemicalPotential`).
template<ScalarPhysical T = double>
using ReducedChemicalPotential = detail::StrongScalar<ReducedChemicalPotentialTag, T>;

/// Lindhard susceptibility returned in natural (DOS-scaled) units.
template<ScalarPhysical T = double>
struct LindhardResult {
    T real_part{};
    T imag_part{};

    [[nodiscard]] constexpr T real() const noexcept { return real_part; }
    [[nodiscard]] constexpr T imag() const noexcept { return imag_part; }
};

// --- compile-time dimensional firewall (WaveVector ≠ Frequency ≠ Temperature) ---

template<typename TagA, typename TagB, ScalarPhysical T>
constexpr bool same_dimension_v =
    std::is_same_v<TagA, TagB>;

template<typename TagA, typename TagB, ScalarPhysical T>
    requires(!same_dimension_v<TagA, TagB, T>)
constexpr detail::StrongScalar<TagA, T> operator+(
    detail::StrongScalar<TagA, T>,
    detail::StrongScalar<TagB, T>) = delete;

template<typename TagA, typename TagB, ScalarPhysical T>
    requires(!same_dimension_v<TagA, TagB, T>)
constexpr detail::StrongScalar<TagA, T> operator-(
    detail::StrongScalar<TagA, T>,
    detail::StrongScalar<TagB, T>) = delete;

}  // namespace sophus
