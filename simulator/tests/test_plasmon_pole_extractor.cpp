#include "core/Concepts.hpp"
#include "engine/PlasmonPoleExtractor.hpp"
#include "physics/ChemicalPotential.hpp"
#include "physics/Constants.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

using namespace sophus;

inline constexpr double ion_mass_ratio = 1836.15267343;
inline constexpr double coulomb_softening = 0.01;

[[nodiscard]] double number_density_from_rs(double rs) noexcept
{
    return 3.0 / (4.0 * constants::pi * rs * rs * rs);
}

[[nodiscard]] double fermi_energy_from_kf(double k_f, double mass) noexcept
{
    return k_f * k_f / (2.0 * mass);
}

[[nodiscard]] double reduced_temperature_from_kelvin(double T_kelvin, double E_F) noexcept
{
    return constants::temperature_kelvin_to_hartree(T_kelvin) / E_F;
}

[[nodiscard]] BarePotentials<> coulomb_potentials_rational(double q_rational) noexcept
{
    const double q_sq = q_rational * q_rational + coulomb_softening * coulomb_softening;
    const double v_coulomb = constants::four_pi / q_sq;
    return BarePotentials<>{
        v_coulomb,
        v_coulomb,
        -v_coulomb,
    };
}

[[nodiscard]] bool build_inputs(double rs,
                                double T_kelvin,
                                double q,
                                PlasmonPoleInputs<>& out) noexcept
{
    const double density = number_density_from_rs(rs);
    const double k_f = constants::fermi_wavevector_from_rs(rs);

    const double E_F_e = fermi_energy_from_kf(k_f, constants::electron_mass_amu);
    const double E_F_i = fermi_energy_from_kf(k_f, ion_mass_ratio);

    const auto gamma_e = find_chemical_potential(ChemicalPotentialInputs<>{
        ReducedTemperature<>{reduced_temperature_from_kelvin(T_kelvin, E_F_e)},
        Density<>{density},
        WaveVector<>{k_f},
    });
    const auto gamma_i = find_chemical_potential(ChemicalPotentialInputs<>{
        ReducedTemperature<>{reduced_temperature_from_kelvin(T_kelvin, E_F_i)},
        Density<>{density},
        WaveVector<>{k_f},
    });

    if (!gamma_e || !gamma_i) {
        return false;
    }

    out = PlasmonPoleInputs<>{
        .q = WaveVector<>{q},
        .electron = PlasmonSpeciesParams<>{
            .tau = ReducedTemperature<>{reduced_temperature_from_kelvin(T_kelvin, E_F_e)},
            .gamma = *gamma_e,
            .fermi_energy = E_F_e,
            .mass = constants::electron_mass_amu,
        },
        .ion = PlasmonSpeciesParams<>{
            .tau = ReducedTemperature<>{reduced_temperature_from_kelvin(T_kelvin, E_F_i)},
            .gamma = *gamma_i,
            .fermi_energy = E_F_i,
            .mass = ion_mass_ratio,
        },
        .potentials = coulomb_potentials_rational(q),
        .wigner_seitz_radius = rs,
    };
    return true;
}

void test_bohm_gross_positive()
{
    const double rs = 2.0;
    const double tau = 0.05;
    const double q = 0.5;

    const double omega_bg =
        PlasmonPoleExtractor::bohm_gross_frequency(WaveVector<>{q}, ReducedTemperature<>{tau}, rs);

    assert(std::isfinite(omega_bg));
    assert(omega_bg > 0.0);
    std::cout << "Bohm–Gross estimate at q=" << q << ": omega_BG=" << omega_bg << '\n';
}

void test_brent_on_known_root()
{
    const auto linear = [](double x) noexcept { return x - 1.25; };
    const auto root = brent_root(linear, 0.5, 2.0);
    assert(root.has_value());
    assert(std::abs(*root - 1.25) < 1.0e-10);
    std::cout << "Brent verified on linear function: root=" << *root << '\n';
}

void test_plasmon_extraction()
{
    const double rs = 2.0;
    const double T_kelvin = 10000.0;
    const double q = 0.5;

    PlasmonPoleInputs<> inputs{};
    assert(build_inputs(rs, T_kelvin, q, inputs));

    const auto state = PlasmonPoleExtractor::extract(inputs);
    assert(state.has_value());

    const double omega_p = state->omega_raw();
    const double re_at_pole =
        PlasmonPoleExtractor::evaluate_epsilon(inputs, state->omega_p).real();

    assert(std::isfinite(omega_p));
    assert(std::isfinite(state->landau_damping));
    assert(std::abs(re_at_pole) < 1.0e-8);

    std::cout << "Plasmon pole at q=" << q << ": omega_p=" << omega_p
              << "  Im(eps)=" << state->landau_damping
              << "  Re(eps)=" << re_at_pole << '\n';
}

}  // namespace

int main()
{
    test_bohm_gross_positive();
    test_brent_on_known_root();
    test_plasmon_extraction();
    std::cout << "PlasmonPoleExtractor tests passed.\n";
    return 0;
}
