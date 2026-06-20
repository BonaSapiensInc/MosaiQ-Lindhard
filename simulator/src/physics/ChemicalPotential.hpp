#pragma once

#include "core/Concepts.hpp"
#include "core/FermiDirac.hpp"
#include "core/RootFinder.hpp"
#include "physics/Constants.hpp"

#include <cmath>
#include <limits>
#include <optional>

namespace sophus {

struct DensityTag {};

/// Number density n (legacy: `_density` in atomic units).
template<ScalarPhysical T = double>
using Density = detail::StrongScalar<DensityTag, T>;

template<ScalarPhysical T = double>
struct ChemicalPotentialInputs {
    ReducedTemperature<T> tau;
    Density<T> density;
    WaveVector<T> fermi_wavevector;  ///< k_F
};

/// Legacy density prefactor: τ^{3/2} k_F^3 / (2π²).
template<ScalarPhysical T = double>
[[nodiscard]] T density_prefactor(ReducedTemperature<T> tau, WaveVector<T> k_f) noexcept
{
    const T k_cube = k_f.raw() * k_f.raw() * k_f.raw();
    return std::pow(tau.raw(), T{1.5}) * k_cube / (T{2} * constants::pi * constants::pi);
}

/// Invert F_{1/2}(η) = n / prefactor via Newton–Raphson (legacy `set_ChemicalPotential`).
template<ScalarPhysical T = double>
[[nodiscard]] std::optional<ReducedChemicalPotential<T>>
find_chemical_potential(ChemicalPotentialInputs<T> inputs,
                        RootFinderPolicy<T> policy = {})
{
    const T tau = inputs.tau.raw();
    const T prefactor = density_prefactor(inputs.tau, inputs.fermi_wavevector);
    const T target = inputs.density.raw() / prefactor;

    // Legacy: γ = 1 ⇒ η = 1/τ; classical (large τ, small target) needs η ≪ 0.
    T eta_seed = T{1} / tau;
    const T F_at_seed =
        fermi_dirac_integral(FermiDiracOrder::Half, FermiDiracEta<T>{eta_seed});
    if (F_at_seed > target) {
        eta_seed = std::log(std::max(target, std::numeric_limits<T>::min()));
    } else {
        for (std::size_t expand = 0;
             expand < 64 &&
             fermi_dirac_integral(FermiDiracOrder::Half, FermiDiracEta<T>{eta_seed}) < target;
             ++expand) {
            eta_seed *= T{2};
        }
    }

    const auto eval = [&](T eta) -> NewtonEval<T> {
        const FermiDiracEta<T> eta_wrapped{eta};
        return NewtonEval<T>{
            fermi_dirac_integral(FermiDiracOrder::Half, eta_wrapped),
            fermi_dirac_integral(FermiDiracOrder::NegativeHalf, eta_wrapped),
        };
    };

    const std::optional<T> eta_root = newton_raphson(eval, eta_seed, target, policy);
    if (!eta_root) {
        return std::nullopt;
    }

    return ReducedChemicalPotential<T>{(*eta_root) * tau};
}

}  // namespace sophus
