#pragma once

#include "core/Concepts.hpp"

namespace sophus {

/// F_{1/2}(η) and F_{-1/2}(η) orders used in chemical-potential inversion.
enum class FermiDiracOrder {
    Half,          ///< order +1/2  (legacy: orderNumerator =  1, orderDenominator = 2)
    NegativeHalf,  ///< order −1/2  (legacy: orderNumerator = -1, orderDenominator = 2)
};

/// Dimensionless argument η = γ / τ entering the legacy `theFermiDiracIntegral`.
template<ScalarPhysical T = double>
struct FermiDiracEtaTag {};
template<ScalarPhysical T = double>
using FermiDiracEta = detail::StrongScalar<FermiDiracEtaTag<T>, T>;

template<ScalarPhysical T = double>
[[nodiscard]] FermiDiracEta<T> make_fermi_dirac_eta(ReducedChemicalPotential<T> gamma,
                                                    ReducedTemperature<T> tau) noexcept
{
    return FermiDiracEta<T>{gamma.raw() / tau.raw()};
}

/// Mohankumar / sinc-strip Fermi–Dirac integral (legacy `theFermiDiracIntegral`).
template<ScalarPhysical T = double>
[[nodiscard]] T fermi_dirac_integral(FermiDiracOrder order, FermiDiracEta<T> eta);

}  // namespace sophus
