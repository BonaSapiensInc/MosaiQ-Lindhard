#include "core/SincQuadrature.hpp"

#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>
#include <ranges>

namespace sophus {

template<ScalarPhysical T>
T asinh_sinc_abscissa(T kh) noexcept
{
    return std::log(std::exp(kh) + std::sqrt(T{1} + std::exp(T{2} * kh)));
}

template<ScalarPhysical T>
T sinc_weight_factor(T kh) noexcept
{
    return T{1} / std::sqrt(T{1} + std::exp(-T{2} * kh));
}

template<ScalarPhysical T>
GvLindhardIntegrand<T> make_gv_lindhard_integrand()
{
    return [](const SincSample<T>& sample, const GvLindhardSincArgs<T>& args) -> T {
        const T q = args.q.raw();
        const T omega = args.omega.raw();
        const T tau = args.tau.raw();
        const T gamma = args.gamma.raw();

        const T nu_p = ((omega / q) + q) / T{2};
        const T nu_n = ((omega / q) - q) / T{2};

        const T x = sample.x;

        T argn = (x - nu_n) / (x + nu_n);
        argn = std::sqrt(argn * argn);
        T argp = (x - nu_p) / (x + nu_p);
        argp = std::sqrt(argp * argp);

        T flog = std::log(argn) - std::log(argp);
        flog /= (std::exp((x * x - gamma) / tau) + T{1});
        flog *= x;

        return sample.weight * flog;
    };
}

template<ScalarPhysical T>
SincQuadrature<T>::SincQuadrature(SincQuadratureParams<T> params)
    : params_{params}
{
}

template<ScalarPhysical T>
SincQuadratureParams<T> SincQuadrature<T>::params() const noexcept
{
    return params_;
}

template<ScalarPhysical T>
auto SincQuadrature<T>::node_indices() const
{
    return std::views::iota(-static_cast<long>(params_.half_nodes),
                            static_cast<long>(params_.half_nodes) + 1);
}

template<ScalarPhysical T>
SincSample<T> SincQuadrature<T>::sample_at(long k) const
{
    const T h = params_.step();
    const T kh = h * static_cast<T>(k);
    return SincSample<T>{
        k,
        kh,
        asinh_sinc_abscissa(kh),
        sinc_weight_factor(kh),
    };
}

template<ScalarPhysical T>
T SincQuadrature<T>::integrate(GvLindhardIntegrand<T> integrand,
                               const GvLindhardSincArgs<T>& args) const
{
    const auto indices = node_indices();
    return std::transform_reduce(
        std::ranges::begin(indices),
        std::ranges::end(indices),
        T{0},
        std::plus{},
        [&](long k) { return integrand(sample_at(k), args); });
}

template<ScalarPhysical T>
T SincQuadrature<T>::evaluate_gv_real_lindhard(const GvLindhardSincArgs<T>& args) const
{
    const T h = params_.step();
    const T sum = integrate(make_gv_lindhard_integrand<T>(), args);
    return sum * h / (-args.q.raw());
}

template<ScalarPhysical T>
std::vector<SincSample<T>> SincQuadrature<T>::all_samples() const
{
    std::vector<SincSample<T>> samples;
    samples.reserve(params_.node_count());
    for (long k : node_indices()) {
        samples.push_back(sample_at(k));
    }
    return samples;
}

template double asinh_sinc_abscissa<double>(double) noexcept;
template double sinc_weight_factor<double>(double) noexcept;
template GvLindhardIntegrand<double> make_gv_lindhard_integrand<double>();

template class SincQuadrature<double>;

}  // namespace sophus
