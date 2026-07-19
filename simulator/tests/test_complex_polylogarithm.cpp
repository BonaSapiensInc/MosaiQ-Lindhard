/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "core/ComplexPolylogarithm.hpp"
#include "core/RiemannZetaBorwein.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

namespace {

using namespace mosaiq;

void assert_close_c(const char* label,
                    std::complex<double> got,
                    std::complex<double> expected,
                    double tol)
{
    const double err = std::abs(got - expected);
    assert(std::isfinite(got.real()) && std::isfinite(got.imag()));
    assert(err <= tol);
    std::cout << label << " = (" << got.real() << ", " << got.imag() << ") err=" << err << '\n';
}

void test_zero_and_domain()
{
    const auto z0 = evaluate_complex_polylog(2.0, {0.0, 0.0});
    assert(z0.has_value());
    assert_close_c("Li_2(0)", *z0, {0.0, 0.0}, 1.0e-15);

    assert(!evaluate_complex_polylog(1.0, {0.5, 0.0}).has_value());
    assert(!evaluate_complex_polylog(0.5, {0.5, 0.0}).has_value());
    std::cout << "Domain rejection OK.\n";
}

void test_series_small_z()
{
    // Li_2(0.5) ≈ 0.5822405264650125
    const auto li2 = evaluate_complex_polylog(2.0, {0.5, 0.0});
    assert(li2.has_value());
    assert_close_c("Li_2(0.5)", *li2, {0.5822405264650125, 0.0}, 1.0e-12);
}

void test_unit_boundary_equals_zeta()
{
    const auto li3 = evaluate_complex_polylog(3.0, {1.0, 0.0});
    const auto z3 = riemann_zeta_borwein(3.0);
    assert(li3.has_value() && z3.has_value());
    assert_close_c("Li_3(1)=zeta(3)", *li3, {*z3, 0.0}, 1.0e-14);
}

void test_heavy_path_not_yet()
{
    // Outside series_radius: analytic continuation stub must refuse.
    const auto out = evaluate_complex_polylog(2.0, {0.9, 0.0});
    assert(!out.has_value());
    std::cout << "Heavy-path stub correctly returns nullopt for |z|>=0.75.\n";
}

}  // namespace

int main()
{
    test_zero_and_domain();
    test_series_small_z();
    test_unit_boundary_equals_zeta();
    test_heavy_path_not_yet();
    std::cout << "All ComplexPolylogarithm skeleton tests passed.\n";
    return 0;
}
