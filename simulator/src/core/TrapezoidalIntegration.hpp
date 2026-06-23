#pragma once

#include <cassert>
#include <cstddef>
#include <vector>

namespace mosaiq {

[[nodiscard]] inline double trapezoidal_integrate(const std::vector<double>& x,
                                                  const std::vector<double>& y) noexcept
{
    assert(x.size() == y.size());
    assert(x.size() >= 2);

    double integral = 0.0;
    for (std::size_t i = 1; i < x.size(); ++i) {
        const double dx = x[i] - x[i - 1];
        integral += 0.5 * (y[i] + y[i - 1]) * dx;
    }
    return integral;
}

}  // namespace mosaiq
