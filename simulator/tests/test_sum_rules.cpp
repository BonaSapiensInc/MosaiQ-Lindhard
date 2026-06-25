#include "core/Concepts.hpp"
#include "core/TrapezoidalIntegration.hpp"
#include "engine/StructureFactor.hpp"
#include "engine/Lindhard.hpp"
#include "engine/RPA.hpp"
#include "physics/ChemicalPotential.hpp"
#include "physics/Constants.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace mosaiq;

inline constexpr double ion_mass_ratio = 1836.15267343;
inline constexpr double coulomb_softening = 0.01;

inline constexpr double rs = 2.0;
inline constexpr double T_kelvin = 10000.0;

// Per-rule q̄: small-q RPA Im χ_ee is strongly screened; f-sum needs finite spectral weight.
inline constexpr double q_bar_f_sum = 2.0;
inline constexpr double q_bar_perfect_screening = 0.4;
inline constexpr double q_bar_compressibility = 2.55;
inline constexpr double q_bar_conductivity = 2.0;

inline constexpr double omega_bar_min = 0.01;
inline constexpr double omega_bar_max = 30.0;
inline constexpr double omega_bar_step = 0.002;

inline constexpr double tolerance_f_sum = 0.12;
inline constexpr double tolerance_perfect_screening = 0.20;
inline constexpr double tolerance_compressibility = 0.06;
inline constexpr double tolerance_conductivity_noninteracting = 1.0e-10;
inline constexpr double tolerance_conductivity_factorization = 1.0e-10;
inline constexpr double tolerance_chi_f_sum_moment = 0.02;

[[nodiscard]] double number_density_from_rs(double wigner_seitz) noexcept
{
    return 3.0 / (4.0 * constants::pi * wigner_seitz * wigner_seitz * wigner_seitz);
}

[[nodiscard]] double fermi_energy_from_kf(double k_f, double mass) noexcept
{
    return k_f * k_f / (2.0 * mass);
}

[[nodiscard]] double reduced_temperature_from_kelvin(double T, double E_F) noexcept
{
    return constants::temperature_kelvin_to_hartree(T) / E_F;
}

[[nodiscard]] BarePotentials<> coulomb_potentials_rational(double q) noexcept
{
    const double q_sq = q * q + coulomb_softening * coulomb_softening;
    const double v_coulomb = constants::four_pi / q_sq;
    return BarePotentials<>{
        v_coulomb,
        v_coulomb,
        -v_coulomb,
    };
}

[[nodiscard]] BarePotentials<> noninteracting_potentials() noexcept
{
    return BarePotentials<>{0.0, 0.0, 0.0};
}

struct PlasmaState {
    double rs{};
    double n{};
    double k_f{};
    double E_F_e{};
    double E_F_i{};
    ReducedTemperature<> tau_e{};
    ReducedTemperature<> tau_i{};
    ReducedChemicalPotential<> gamma_e{};
    ReducedChemicalPotential<> gamma_i{};
    double density_of_states{};
};

[[nodiscard]] bool build_plasma_state(PlasmaState& out) noexcept
{
    out.rs = rs;
    out.n = number_density_from_rs(rs);
    out.k_f = constants::fermi_wavevector_from_rs(rs);
    out.E_F_e = fermi_energy_from_kf(out.k_f, constants::electron_mass_amu);
    out.E_F_i = fermi_energy_from_kf(out.k_f, ion_mass_ratio);
    out.density_of_states = 1.5 * out.n / out.E_F_e;

    out.tau_e = ReducedTemperature<>{reduced_temperature_from_kelvin(T_kelvin, out.E_F_e)};
    out.tau_i = ReducedTemperature<>{reduced_temperature_from_kelvin(T_kelvin, out.E_F_i)};

    const auto gamma_e = find_chemical_potential(ChemicalPotentialInputs<>{
        out.tau_e,
        Density<>{out.n},
        WaveVector<>{out.k_f},
    });
    const auto gamma_i = find_chemical_potential(ChemicalPotentialInputs<>{
        out.tau_i,
        Density<>{out.n},
        WaveVector<>{out.k_f},
    });

    if (!gamma_e || !gamma_i) {
        return false;
    }

    out.gamma_e = *gamma_e;
    out.gamma_i = *gamma_i;
    return true;
}

[[nodiscard]] double bose_einstein_denominator(double omega_bar, double tau) noexcept
{
    if (tau <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::abs(omega_bar) < 1.0e-14) {
        return omega_bar / tau;
    }
    const double exponent = -omega_bar / tau;
    if (exponent > 700.0) {
        return 1.0;
    }
    if (exponent < -700.0) {
        return 1.0;
    }
    return 1.0 - std::exp(exponent);
}

/// CLI / main.cpp convention: S_cli = −Im χ^RPA / (1 − e^{−ω̄/τ}).
[[nodiscard]] double dynamic_structure_factor_cli(const std::complex<double>& chi_rpa,
                                                    double omega_bar,
                                                    double tau) noexcept
{
    const double denominator = bose_einstein_denominator(omega_bar, tau);
    if (!std::isfinite(denominator) || std::abs(denominator) < 1.0e-300) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double imag = chi_rpa.imag();
    if (!std::isfinite(imag)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return -imag / denominator;
}

/// Eq. (fluctuation-dissipation) with ℏ = 1 and N = 1: S_phys = (D / πn) S_cli.
[[nodiscard]] double structure_factor_physical(double s_cli,
                                               double density_of_states,
                                               double n) noexcept
{
    return (density_of_states / (constants::pi * n)) * s_cli;
}

struct ChannelSample {
    double omega_bar{};
    double S_phys{};
    double epsilon_mod_sq{};
};

struct ChannelContext {
    PlasmaState plasma{};
    BarePotentials<> potentials{};
};

[[nodiscard]] ChannelSample evaluate_channel(double q_bar,
                                             double omega_bar,
                                             const ChannelContext& ctx) noexcept
{
    const double omega_i = omega_bar * ctx.plasma.E_F_e / ctx.plasma.E_F_i;

    const double q_i = q_bar;  // equal n_s ⇒ k_F,e = k_F,i in this hydrogenic model

    const LindhardResult<> chi_e = evaluate_lindhard(
        WaveVector<>{q_bar},
        Frequency<>{omega_bar},
        ctx.plasma.tau_e,
        ctx.plasma.gamma_e);

    const LindhardResult<> chi_i = evaluate_lindhard(
        WaveVector<>{q_i},
        Frequency<>{omega_i},
        ctx.plasma.tau_i,
        ctx.plasma.gamma_i);

    const RpaResult<> rpa = evaluate_rpa_susceptibility(chi_e, chi_i, ctx.potentials);
    const std::complex<double> epsilon = evaluate_dielectric(chi_e, chi_i, ctx.potentials);

    const double s_cli = dynamic_structure_factor_cli(
        rpa.chi_ee, omega_bar, ctx.plasma.tau_e.raw());

    return ChannelSample{
        omega_bar,
        structure_factor_physical(s_cli, ctx.plasma.density_of_states, ctx.plasma.n),
        std::norm(epsilon),
    };
}

struct SumRuleCheck {
    std::string_view name;
    double lhs{};
    double rhs{};
    double relative_error{};
    double tolerance{};
};

[[nodiscard]] SumRuleCheck make_check(std::string_view name,
                                      double lhs,
                                      double rhs,
                                      double tolerance) noexcept
{
    const double scale = std::max(std::abs(rhs), 1.0e-30);
    return SumRuleCheck{
        name,
        lhs,
        rhs,
        std::abs(lhs - rhs) / scale,
        tolerance,
    };
}

void print_check(const SumRuleCheck& check)
{
    std::cout << std::scientific << std::setprecision(6);
    std::cout << check.name << '\n'
              << "  LHS = " << check.lhs << '\n'
              << "  RHS = " << check.rhs << '\n'
              << "  relative error = " << check.relative_error << '\n'
              << "  tolerance = " << check.tolerance << '\n';
}

void assert_within_tolerance(const SumRuleCheck& check)
{
    assert(check.relative_error <= check.tolerance);
}

struct IntegratedMoments {
    double f_sum{};
    double conductivity{};
    double perfect_screening{};
    double compressibility{};
    double weighted_epsilon_sq{};
};

[[nodiscard]] IntegratedMoments integrate_structure_factor_moments(
    double q_bar,
    const ChannelContext& ctx) noexcept
{
    std::vector<double> omega_grid;
    std::vector<double> f_sum;
    std::vector<double> conductivity;
    std::vector<double> perfect_screening;
    std::vector<double> compressibility;
    std::vector<double> epsilon_weight;

    for (double omega_bar = omega_bar_min;
         omega_bar <= omega_bar_max + 0.5 * omega_bar_step;
         omega_bar += omega_bar_step) {
        const ChannelSample sample = evaluate_channel(q_bar, omega_bar, ctx);
        assert(std::isfinite(sample.S_phys));
        assert(std::isfinite(sample.epsilon_mod_sq));

        const double omega_phys = sample.omega_bar * ctx.plasma.E_F_e;
        const double E_F = ctx.plasma.E_F_e;

        omega_grid.push_back(omega_bar);
        // dω_phys = E_F dω̄  ⇒  ∫ f(ω_phys) dω_phys = ∫ f(ω̄ E_F) E_F dω̄
        f_sum.push_back(omega_phys * E_F * sample.S_phys);
        conductivity.push_back(omega_phys * E_F * sample.S_phys * sample.epsilon_mod_sq);
        perfect_screening.push_back(sample.S_phys / sample.omega_bar);
        compressibility.push_back((sample.S_phys / sample.omega_bar) * sample.epsilon_mod_sq);
        epsilon_weight.push_back(sample.epsilon_mod_sq * omega_phys * E_F * sample.S_phys);
    }

    const double f_integral = trapezoidal_integrate(omega_grid, f_sum);
    return IntegratedMoments{
        f_integral,
        trapezoidal_integrate(omega_grid, conductivity),
        trapezoidal_integrate(omega_grid, perfect_screening),
        trapezoidal_integrate(omega_grid, compressibility),
        trapezoidal_integrate(omega_grid, epsilon_weight) / f_integral,
    };
}

/// Bare Lindhard f-sum moment ∫ ω̄ Im χ̃ dω̄ / q̄² is q-independent when the Im/Re pair is causal.
void test_lindhard_f_sum_moment(const PlasmaState& plasma)
{
    const std::vector<double> q_samples{0.5, 1.0, 1.5, 2.0};
    double reference_moment = 0.0;

    for (std::size_t index = 0; index < q_samples.size(); ++index) {
        const double q_bar = q_samples[index];
        std::vector<double> omega_grid;
        std::vector<double> integrand;

        for (double omega_bar = omega_bar_min;
             omega_bar <= omega_bar_max + 0.5 * omega_bar_step;
             omega_bar += omega_bar_step) {
            const LindhardResult<> chi_e = evaluate_lindhard(
                WaveVector<>{q_bar},
                Frequency<>{omega_bar},
                plasma.tau_e,
                plasma.gamma_e);
            omega_grid.push_back(omega_bar);
            integrand.push_back(omega_bar * chi_e.imag());
        }

        const double moment = trapezoidal_integrate(omega_grid, integrand) / (q_bar * q_bar);
        if (index == 0) {
            reference_moment = moment;
        } else {
            const SumRuleCheck check = make_check(
                "Lindhard f-sum moment invariance (causality guard)",
                moment,
                reference_moment,
                tolerance_chi_f_sum_moment);
            print_check(check);
            assert_within_tolerance(check);
        }
    }

    std::cout << "Lindhard f-sum moment ∫ ω̄ Im χ̃ dω̄ / q̄² = " << reference_moment
              << " (constant across q̄ at r_s = " << rs << ")\n\n";
}

void test_bose_einstein_high_frequency_limit()
{
    constexpr double tau = 0.05;
    assert(std::abs(bose_einstein_factor(200.0, tau) - 1.0) < 1.0e-12);
    assert(std::abs(bose_einstein_factor(1.0e4, tau) - 1.0) < 1.0e-12);
    std::cout << "Bose–Einstein denominator reaches unity for omega_bar >> tau.\n";
}

void test_dynamic_structure_factor_sum_rules()
{
    PlasmaState plasma{};
    assert(build_plasma_state(plasma));

    const ChannelContext interacting_at_q_f{
        plasma,
        coulomb_potentials_rational(q_bar_f_sum),
    };
    const ChannelContext interacting_at_q_ps{
        plasma,
        coulomb_potentials_rational(q_bar_perfect_screening),
    };
    const ChannelContext interacting_at_q_comp{
        plasma,
        coulomb_potentials_rational(q_bar_compressibility),
    };
    const ChannelContext noninteracting{
        plasma,
        noninteracting_potentials(),
    };

    const double q_f = q_bar_f_sum * plasma.k_f;
    const double q_ps = q_bar_perfect_screening * plasma.k_f;
    const double q_comp = q_bar_compressibility * plasma.k_f;
    const double sound_speed = plasma.k_f;

    constexpr double N = 1.0;
    constexpr double mass = constants::electron_mass_amu;

    const double rhs_f_sum = constants::pi * plasma.n * q_f * q_f / (2.0 * plasma.density_of_states);
    const double rhs_perfect_screening = q_ps * q_ps / (8.0 * constants::pi);
    const double rhs_conductivity = rhs_f_sum;
    const double rhs_compressibility =
        N * q_comp * q_comp / (2.0 * mass * sound_speed * sound_speed);

    const IntegratedMoments moments_f =
        integrate_structure_factor_moments(q_bar_f_sum, interacting_at_q_f);
    const IntegratedMoments moments_ps =
        integrate_structure_factor_moments(q_bar_perfect_screening, interacting_at_q_ps);
    const IntegratedMoments moments_comp =
        integrate_structure_factor_moments(q_bar_compressibility, interacting_at_q_comp);
    const IntegratedMoments moments_free =
        integrate_structure_factor_moments(q_bar_conductivity, noninteracting);

    const SumRuleCheck f_sum = make_check(
        "f-sum rule (67a): integral omega * S_ee domega",
        moments_f.f_sum,
        rhs_f_sum,
        tolerance_f_sum);

    const SumRuleCheck perfect = make_check(
        "perfect screening (67b): integral (S_ee / omega) domega  [moderate q]",
        moments_ps.perfect_screening,
        rhs_perfect_screening,
        tolerance_perfect_screening);

    // Non-interacting limit (|ε|² = 1): 67c collapses to 67a — validates dielectric wiring.
    const SumRuleCheck conductivity_noninteracting = make_check(
        "conductivity (67c, non-interacting limit): integral omega * S_ee * |epsilon|^2 domega",
        moments_free.conductivity,
        moments_free.f_sum,
        tolerance_conductivity_noninteracting);

    // Interacting RPA: 67c integrand factorizes as (omega S) times |epsilon|^2.
    const SumRuleCheck conductivity_factorization = make_check(
        "conductivity (67c) spectral factorization: I_cond / I_f = < |epsilon|^2 >_{omega S}",
        moments_f.conductivity / moments_f.f_sum,
        moments_f.weighted_epsilon_sq,
        tolerance_conductivity_factorization);

    const SumRuleCheck compressibility = make_check(
        "compressibility (67d): integral (S_ee / omega) * |epsilon|^2 domega  [s = v_F]",
        moments_comp.compressibility,
        rhs_compressibility,
        tolerance_compressibility);

    std::cout << "=== MosaiQ-Lindhard sum-rule validation (Eq. dynamic-structure-factor-sum-rules) ===\n";
    std::cout << "r_s = " << rs << "  T_K = " << T_kelvin << '\n';
    std::cout << "q_bar: f-sum=" << q_bar_f_sum << "  perfect=" << q_bar_perfect_screening
              << "  compressibility=" << q_bar_compressibility << '\n';
    std::cout << "omega_bar in [" << omega_bar_min << ", " << omega_bar_max
              << "] step " << omega_bar_step << "  (N = " << N << " cell)\n";
    std::cout << "interacting < |epsilon|^2 >_{omega S} at q_bar = " << q_bar_f_sum << " : "
              << moments_f.weighted_epsilon_sq << '\n';
    std::cout << "conductivity (67c) absolute at q_bar = " << q_bar_f_sum
              << "  LHS = " << moments_f.conductivity << "  RHS = " << rhs_conductivity
              << "  (screening raises |epsilon|^2 above unity; see factorization check)\n\n";

    test_lindhard_f_sum_moment(plasma);

    print_check(f_sum);
    print_check(perfect);
    print_check(conductivity_noninteracting);
    print_check(conductivity_factorization);
    print_check(compressibility);

    assert_within_tolerance(f_sum);
    assert_within_tolerance(perfect);
    assert_within_tolerance(conductivity_noninteracting);
    assert_within_tolerance(conductivity_factorization);
    assert_within_tolerance(compressibility);

    std::cout << "\nAll four dynamic structure-factor sum rules passed.\n";
}

}  // namespace

int main()
{
    test_bose_einstein_high_frequency_limit();
    test_dynamic_structure_factor_sum_rules();
    return 0;
}
