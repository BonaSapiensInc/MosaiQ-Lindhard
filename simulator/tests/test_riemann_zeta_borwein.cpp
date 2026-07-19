/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "core/RiemannZetaBorwein.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

using namespace mosaiq;

void assert_close(const char* label, double got, double expected, double tol)
{
    const double err = std::abs(got - expected);
    const double rel = (expected != 0.0) ? err / std::abs(expected) : err;
    assert(std::isfinite(got));
    assert(rel <= tol || err <= tol);
    std::cout << label << " = " << got << " (expected " << expected << ", rel=" << rel << ")\n";
}

void test_even_integers()
{
    const BorweinPolicy policy{};
    const auto z2 = riemann_zeta_borwein(2.0, policy);
    const auto z4 = riemann_zeta_borwein(4.0, policy);
    assert(z2.has_value());
    assert(z4.has_value());
    assert_close("zeta(2)", *z2, riemann_zeta_even_integer<double>(1), 1.0e-14);
    assert_close("zeta(4)", *z4, riemann_zeta_even_integer<double>(2), 1.0e-14);
}

void test_apery()
{
    // Apery's constant ζ(3) ≈ 1.2020569031595942
    const auto z3 = riemann_zeta_borwein(3.0);
    assert(z3.has_value());
    assert_close("zeta(3)", *z3, 1.2020569031595942, 1.0e-14);
}

void test_non_integer()
{
    // ζ(2.5) reference (high-precision truncated)
    const auto z = riemann_zeta_borwein(2.5);
    assert(z.has_value());
    assert_close("zeta(2.5)", *z, 1.3414872572509171, 1.0e-12);
}

void test_bit_identity()
{
    const auto a = riemann_zeta_borwein(3.0);
    const auto b = riemann_zeta_borwein(3.0);
    assert(a.has_value() && b.has_value());
    assert(*a == *b);
    std::cout << "Bit-identical repeated zeta(3) OK.\n";
}

void test_domain_rejection()
{
    assert(!riemann_zeta_borwein(1.0).has_value());
    assert(!riemann_zeta_borwein(0.5).has_value());
    assert(!riemann_zeta_borwein(-1.0).has_value());
    assert(!riemann_zeta_borwein(2.0, BorweinPolicy{.truncation_order = 0}).has_value());
    std::cout << "Domain rejection OK.\n";
}

}  // namespace

int main()
{
    test_even_integers();
    test_apery();
    test_non_integer();
    test_bit_identity();
    test_domain_rejection();
    std::cout << "All RiemannZetaBorwein tests passed.\n";
    return 0;
}
