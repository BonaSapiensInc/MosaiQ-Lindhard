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
