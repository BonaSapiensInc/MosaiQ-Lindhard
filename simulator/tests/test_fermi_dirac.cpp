#include "core/Concepts.hpp"
#include "core/FermiDirac.hpp"
#include "physics/ChemicalPotential.hpp"
#include "physics/Constants.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

using namespace mosaiq;

void test_fermi_dirac_half_at_zero_eta()
{
    // Legacy Mohankumar/sinc-strip quadrature at η = 0 (reference value ≈ 0.678).
    const double F = fermi_dirac_integral(FermiDiracOrder::Half, FermiDiracEta<>{0.0});
    assert(std::isfinite(F));
    assert(F > 0.5 && F < 1.0);
    std::cout << "F_{1/2}(0) = " << F << '\n';
}

void test_fermi_dirac_negative_half_at_zero_eta()
{
    const double Fd = fermi_dirac_integral(FermiDiracOrder::NegativeHalf, FermiDiracEta<>{0.0});
    assert(std::isfinite(Fd));
    assert(Fd > 0.0);
    std::cout << "F_{-1/2}(0) = " << Fd << '\n';
}

void test_chemical_potential_zero_temperature_limit()
{
    // Legacy: at T→0, μ → E_F ⇒ γ → 1 for the unperturbed electron gas (rs = 3).
    constexpr double rs = 3.0;
    const double k_f = constants::fermi_wavevector_from_rs(rs);
    const double density_value = 3.0 / (4.0 * constants::pi * rs * rs * rs);

    const ReducedTemperature<> tau{1.0e-4};

    const auto gamma = find_chemical_potential(ChemicalPotentialInputs<>{
        tau,
        Density<>{density_value},
        WaveVector<>{k_f},
    });

    assert(gamma.has_value());
    assert(std::isfinite(gamma->raw()));
    // Low-τ limit: γ remains O(1) and approaches unity from above for rs = 3.
    assert(gamma->raw() > 1.0);
    assert(gamma->raw() < 1.5);
    std::cout << "γ (τ=1e-4, rs=3) = " << gamma->raw() << '\n';
}

void test_chemical_potential_finite_temperature()
{
    const double k_f = constants::fermi_wavevector_from_rs(3.0);
    const ReducedTemperature<> tau{0.05};
    const double prefactor = density_prefactor(tau, WaveVector<>{k_f});

    const auto gamma = find_chemical_potential(ChemicalPotentialInputs<>{
        tau,
        Density<>{prefactor * fermi_dirac_integral(FermiDiracOrder::Half, FermiDiracEta<>{1.0})},
        WaveVector<>{k_f},
    });

    assert(gamma.has_value());
    assert(std::isfinite(gamma->raw()));
    assert(gamma->raw() > 0.0);
    assert(gamma->raw() < 2.0);
    std::cout << "γ (τ=0.05, inverted) = " << gamma->raw() << '\n';
}

void test_newton_roundtrip_consistency()
{
    const double k_f = constants::fermi_wavevector_from_rs(2.5);
    const ReducedTemperature<> tau{0.08};
    const FermiDiracEta<> eta{1.2};
    const double target =
        fermi_dirac_integral(FermiDiracOrder::Half, eta);

    const auto gamma = find_chemical_potential(ChemicalPotentialInputs<>{
        tau,
        Density<>{density_prefactor(tau, WaveVector<>{k_f}) * target},
        WaveVector<>{k_f},
    });

    assert(gamma.has_value());
    const FermiDiracEta<> recovered_eta{gamma->raw() / tau.raw()};
    const double reproduced = fermi_dirac_integral(FermiDiracOrder::Half, recovered_eta);
    assert(std::abs(reproduced - target) < 1.0e-10);
    std::cout << "Newton round-trip η = " << recovered_eta.raw() << '\n';
}

}  // namespace

int main()
{
    test_fermi_dirac_half_at_zero_eta();
    test_fermi_dirac_negative_half_at_zero_eta();
    test_chemical_potential_zero_temperature_limit();
    test_chemical_potential_finite_temperature();
    test_newton_roundtrip_consistency();

    std::cout << "All Fermi–Dirac / chemical-potential tests passed.\n";
    return 0;
}
