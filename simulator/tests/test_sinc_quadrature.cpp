#include "core/Concepts.hpp"
#include "core/SincQuadrature.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

void assert_finite(const char* label, double value)
{
    assert(!std::isnan(value));
    assert(!std::isinf(value));
    assert(std::isfinite(value));
    std::cout << label << " = " << value << '\n';
}

void test_sample_pipeline()
{
    mosaiq::SincQuadrature<> quad;
    const auto samples = quad.all_samples();

    assert(samples.size() == mosaiq::constants::default_gv_sinc_nodes * 2 + 1);

    for (const auto& sample : samples) {
        assert(std::isfinite(sample.kh));
        assert(std::isfinite(sample.x));
        assert(std::isfinite(sample.weight));
        assert(sample.weight > 0.0);
        assert(sample.x >= 0.0);
    }
}

void test_gv_real_lindhard_finite_domain()
{
    mosaiq::GvLindhardSincArgs<> args{
        mosaiq::WaveVector<>{1.0},
        mosaiq::Frequency<>{0.5},
        mosaiq::ReducedTemperature<>{0.05},
        mosaiq::ReducedChemicalPotential<>{1.0},
    };

    mosaiq::SincQuadrature<> quad;
    const double result = quad.evaluate_gv_real_lindhard(args);

    assert_finite("GV Re chi^L (rational units, q=1, w=0.5, tau=0.05, gamma=1)", result);

    // Physical guardrail: rational-unit real part should remain O(1), not explode.
    assert(std::abs(result) < 10.0);
    assert(result < 0.0);  // legacy prefactor divides by -q with positive q
}

void test_gv_real_lindhard_low_temperature()
{
    mosaiq::GvLindhardSincArgs<> args{
        mosaiq::WaveVector<>{0.8},
        mosaiq::Frequency<>{0.2},
        mosaiq::ReducedTemperature<>{0.01},
        mosaiq::ReducedChemicalPotential<>{1.0},
    };

    mosaiq::SincQuadrature<> quad;
    const double result = quad.evaluate_gv_real_lindhard(args);

    assert_finite("GV Re chi^L (low tau)", result);
    assert(std::abs(result) < 10.0);
}

void test_integrate_matches_evaluate()
{
    mosaiq::GvLindhardSincArgs<> args{
        mosaiq::WaveVector<>{1.2},
        mosaiq::Frequency<>{0.7},
        mosaiq::ReducedTemperature<>{0.08},
        mosaiq::ReducedChemicalPotential<>{0.9},
    };

    mosaiq::SincQuadrature<> quad;
    const double h = quad.params().step();
    const double integrated =
        quad.integrate(mosaiq::make_gv_lindhard_integrand<double>(), args) * h / (-args.q.raw());
    const double evaluated = quad.evaluate_gv_real_lindhard(args);

    assert_finite("integrate pipeline", integrated);
    assert_finite("evaluate_gv_real_lindhard", evaluated);
    assert(std::abs(integrated - evaluated) < 1.0e-12);
}

void test_cauchy_principal_value_kernel_lhopital()
{
    constexpr double pi_over_h = 42.0;
    const double at_resonance = mosaiq::cauchy_principal_value_kernel(0.0, pi_over_h);
    const double near_resonance = mosaiq::cauchy_principal_value_kernel(1.0e-14, pi_over_h);

    assert(at_resonance == 0.0);
    assert(near_resonance == 0.0);
    assert(std::isfinite(mosaiq::cauchy_principal_value_kernel(0.25, pi_over_h)));
}

}  // namespace

int main()
{
    test_cauchy_principal_value_kernel_lhopital();
    test_sample_pipeline();
    test_gv_real_lindhard_finite_domain();
    test_gv_real_lindhard_low_temperature();
    test_integrate_matches_evaluate();

    std::cout << "All SincQuadrature tests passed.\n";
    return 0;
}
