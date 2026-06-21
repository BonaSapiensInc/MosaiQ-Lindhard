#include "engine/PlasmonPoleExtractor.hpp"

#include "physics/Constants.hpp"

#include <cmath>

namespace sophus::detail {

[[nodiscard]] double bohm_gross_frequency_impl(double q_bar,
                                               double tau,
                                               double rs,
                                               double mass) noexcept
{
    const double density = 3.0 / (4.0 * constants::pi * rs * rs * rs);
    const double k_f = constants::fermi_wavevector_from_rs(rs);
    const double E_F = k_f * k_f / (2.0 * mass);
    const double omega_p = std::sqrt(4.0 * constants::pi * density / mass);
    const double omega_p_bar = omega_p / E_F;

    const double omega_sq = omega_p_bar * omega_p_bar + 6.0 * tau * q_bar * q_bar;
    return std::sqrt(std::max(omega_sq, 0.0));
}

}  // namespace sophus::detail

namespace sophus {

template<ScalarPhysical T>
T PlasmonPoleExtractor::bohm_gross_frequency(WaveVector<T> q,
                                             ReducedTemperature<T> tau,
                                             T rs,
                                             T mass) noexcept
{
    return static_cast<T>(detail::bohm_gross_frequency_impl(
        static_cast<double>(q.raw()),
        static_cast<double>(tau.raw()),
        static_cast<double>(rs),
        static_cast<double>(mass)));
}

template double PlasmonPoleExtractor::bohm_gross_frequency<double>(
    WaveVector<double>,
    ReducedTemperature<double>,
    double,
    double) noexcept;

}  // namespace sophus
