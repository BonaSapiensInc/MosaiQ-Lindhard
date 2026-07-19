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

    assert(!evaluate_complex_polylog(0.0, {0.5, 0.0}).has_value());
    assert(!evaluate_complex_polylog(-1.0, {0.5, 0.0}).has_value());
    // z = 1 requires s > 1 (ζ pole / harmonic divergence for s ≤ 1).
    assert(!evaluate_complex_polylog(1.0, {1.0, 0.0}).has_value());
    assert(!evaluate_complex_polylog(0.5, {1.0, 0.0}).has_value());
    // Branch cut [1, ∞): refuse real z > 1 (exact 1 is ζ(s) for s > 1).
    assert(!evaluate_complex_polylog(2.0, {1.5, 0.0}).has_value());

    // 0 < s ≤ 1 is allowed inside the series disk (Pathway A: s = f).
    const auto li_half = evaluate_complex_polylog(0.5, {0.5, 0.0});
    assert(li_half.has_value());
    assert(std::isfinite(li_half->real()) && std::isfinite(li_half->imag()));
    // Li_1(1/2) = −ln(1/2)
    const auto li1 = evaluate_complex_polylog(1.0, {0.5, 0.0});
    assert(li1.has_value());
    assert_close_c("Li_1(0.5)", *li1, {-std::log(0.5), 0.0}, 1.0e-12);
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

void test_heavy_path_inside_unit_disk()
{
    // |z| = 0.9 ≥ series_radius → sinc Mellin path; high-precision reference.
    const auto li2 = evaluate_complex_polylog(2.0, {0.9, 0.0});
    assert(li2.has_value());
    assert_close_c("Li_2(0.9)", *li2, {1.2997147230049588, 0.0}, 1.0e-9);
}

void test_heavy_path_negative_and_exterior()
{
    const auto lim1 = evaluate_complex_polylog(2.0, {-1.0, 0.0});
    assert(lim1.has_value());
    assert_close_c("Li_2(-1)", *lim1, {-0.8224670334241132, 0.0}, 1.0e-9);

    const auto lim2 = evaluate_complex_polylog(2.0, {-2.0, 0.0});
    assert(lim2.has_value());
    assert_close_c("Li_2(-2)", *lim2, {-1.4367463668836808, 0.0}, 1.0e-9);

    const auto lc = evaluate_complex_polylog(2.0, {0.8, 0.3});
    assert(lc.has_value());
    assert_close_c("Li_2(0.8+0.3i)", *lc, {0.9500542362497607, 0.5413399911729875}, 1.0e-9);

    const auto le = evaluate_complex_polylog(2.0, {1.5, 0.2});
    assert(le.has_value());
    assert_close_c("Li_2(1.5+0.2i)", *le, {1.9895165024320542, 1.387727336741074}, 1.0e-7);
}

}  // namespace

int main()
{
    test_zero_and_domain();
    test_series_small_z();
    test_unit_boundary_equals_zeta();
    test_heavy_path_inside_unit_disk();
    test_heavy_path_negative_and_exterior();
    std::cout << "All ComplexPolylogarithm tests passed.\n";
    return 0;
}
