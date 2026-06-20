#include "engine/Lindhard.hpp"

#include "physics/Constants.hpp"

#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <ranges>

namespace sophus {

namespace {

template<ScalarPhysical T>
[[nodiscard]] T branch_xi(Frequency<T> omega, WaveVector<T> q, T sign) noexcept
{
    return ((omega.raw() / q.raw()) + sign * q.raw()) / T{2};
}

template<ScalarPhysical T>
[[nodiscard]] T stable_fermi_log(T arg) noexcept
{
    if (arg >= T{0}) {
        return arg + std::log1p(std::exp(-arg));
    }
    return std::log1p(std::exp(arg));
}

template<ScalarPhysical T>
[[nodiscard]] T fermi_log_kernel(T kh, ReducedChemicalPotential<T> gamma, ReducedTemperature<T> tau) noexcept
{
    const T arg = (gamma.raw() - kh * kh) / tau.raw();
    return stable_fermi_log(arg);
}

template<ScalarPhysical T>
[[nodiscard]] T cauchy_principal_weight(T xi, T kh, T pi_over_h) noexcept
{
    const T argument = pi_over_h * (xi - kh);
    if (argument == T{0}) {
        return T{0};
    }
    return (T{1} - std::cos(argument)) / argument;
}

template<ScalarPhysical T>
[[nodiscard]] T accumulate_cauchy_principal_value(T xi,
                                                   ReducedChemicalPotential<T> gamma,
                                                   ReducedTemperature<T> tau,
                                                   T h,
                                                   T pi_over_h,
                                                   long half_nodes) noexcept
{
    const auto indices = std::views::iota(-half_nodes, half_nodes + 1);
    return std::transform_reduce(
        std::ranges::begin(indices),
        std::ranges::end(indices),
        T{0},
        std::plus{},
        [&](long k) {
            const T kh = static_cast<T>(k) * h;
            const T f_kh = fermi_log_kernel(kh, gamma, tau);
            return f_kh * cauchy_principal_weight(xi, kh, pi_over_h);
        });
}

}  // namespace

template<ScalarPhysical T>
T imaginary_lindhard(WaveVector<T> q,
                     Frequency<T> omega,
                     ReducedTemperature<T> tau,
                     ReducedChemicalPotential<T> gamma)
{
    const T q_raw = q.raw();
    const T tau_raw = tau.raw();
    const T gamma_raw = gamma.raw();

    const T nu_p = branch_xi(omega, q, T{1});
    const T nu_n = branch_xi(omega, q, T{-1});

    const long double nu_p_sq = static_cast<long double>(nu_p) * static_cast<long double>(nu_p);
    const long double nu_n_sq = static_cast<long double>(nu_n) * static_cast<long double>(nu_n);
    const long double gamma_ld = static_cast<long double>(gamma_raw);
    const long double tau_ld = static_cast<long double>(tau_raw);

    const long double arg_exp_n = (nu_n_sq - gamma_ld) / tau_ld;
    const long double arg_exp_p = (nu_p_sq - gamma_ld) / tau_ld;

    const long double numerator = 1.0L + std::exp(-arg_exp_n);
    const long double denominator = 1.0L + std::exp(-arg_exp_p);
    const long double the_log = std::log(numerator / denominator);

    T imag = static_cast<T>(tau_ld * static_cast<long double>(the_log));
    imag *= -(constants::pi / q_raw);
    return imag;
}

template<ScalarPhysical T>
T real_lindhard_kk(WaveVector<T> q,
                   Frequency<T> omega,
                   ReducedTemperature<T> tau,
                   ReducedChemicalPotential<T> gamma)
{
    const T q_raw = q.raw();
    const T tau_raw = tau.raw();

    constexpr long half_nodes = static_cast<long>(constants::default_kk_sinc_nodes);
    const T h = std::sqrt(constants::sinc_strip_factor * constants::pi / static_cast<T>(half_nodes));
    const T pi_over_h = constants::pi / h;

    const T xi_n = branch_xi(omega, q, T{-1});
    const T xi_p = branch_xi(omega, q, T{1});

    const T cpv_n = accumulate_cauchy_principal_value(
        xi_n, gamma, tau, h, pi_over_h, half_nodes);
    const T cpv_p = accumulate_cauchy_principal_value(
        xi_p, gamma, tau, h, pi_over_h, half_nodes);

    return (cpv_n - cpv_p) * (constants::pi * tau_raw) / (T{4} * q_raw);
}

template<ScalarPhysical T>
LindhardResult<T> evaluate_lindhard(WaveVector<T> q,
                                    Frequency<T> omega,
                                    ReducedTemperature<T> tau,
                                    ReducedChemicalPotential<T> gamma)
{
    return LindhardResult<T>{
        real_lindhard_kk(q, omega, tau, gamma),
        imaginary_lindhard(q, omega, tau, gamma),
    };
}

template double imaginary_lindhard<double>(WaveVector<double>,
                                           Frequency<double>,
                                           ReducedTemperature<double>,
                                           ReducedChemicalPotential<double>);
template double real_lindhard_kk<double>(WaveVector<double>,
                                       Frequency<double>,
                                       ReducedTemperature<double>,
                                       ReducedChemicalPotential<double>);
template LindhardResult<double> evaluate_lindhard<double>(WaveVector<double>,
                                                          Frequency<double>,
                                                          ReducedTemperature<double>,
                                                          ReducedChemicalPotential<double>);

}  // namespace sophus
