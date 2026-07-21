/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#include "engine/Berezinian.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>

using mosaiq::ElectronPhononSuperMatrix;
using mosaiq::computeBerezinian;
using Z = std::complex<double>;

namespace {

/// Analytic N_e=2, N_p=1 fixture:
///   A = diag(2, 3), B = [1; 1], C = [1, 1], D = [2]
///   Ã = A − B (C/D) = [[1.5, −0.5], [−0.5, 2.5]], det(Ã)=3.5
///   sdet = det(Ã)/D = 1.75
ElectronPhononSuperMatrix<Z> analytic_supermatrix()
{
    ElectronPhononSuperMatrix<Z> M(2, 1);
    M.A_at(0, 0) = 2.0;
    M.A_at(1, 1) = 3.0;
    M.B_at(0, 0) = 1.0;
    M.B_at(1, 0) = 1.0;
    M.C_at(0, 0) = 1.0;
    M.C_at(0, 1) = 1.0;
    M.D_at(0, 0) = 2.0;
    return M;
}

}  // namespace

TEST(BerezinianTest, AnalyticTwoByOneMatchesHandCalculation)
{
    const auto M = analytic_supermatrix();

    // Hand: X = C/D = [0.5, 0.5]; Ã = A − B X; sdet = det(Ã)/D
    const double d = 2.0;
    const double a00 = 2.0 - (1.0 * 1.0) / d;  // 1.5
    const double a01 = 0.0 - (1.0 * 1.0) / d;  // -0.5
    const double a10 = 0.0 - (1.0 * 1.0) / d;  // -0.5
    const double a11 = 3.0 - (1.0 * 1.0) / d;  // 2.5
    const Z expected = (a00 * a11 - a01 * a10) / d;

    const Z got = computeBerezinian(M, 0.0);
    EXPECT_TRUE(std::isfinite(got.real()));
    EXPECT_TRUE(std::isfinite(got.imag()));
    EXPECT_NEAR(got.real(), expected.real(), 1e-13);
    EXPECT_NEAR(got.imag(), expected.imag(), 1e-13);
    EXPECT_NEAR(got.real(), 1.75, 1e-13);
}

TEST(BerezinianTest, RetardedEtaRegularizesSingularD)
{
    ElectronPhononSuperMatrix<Z> M(2, 1);
    M.A_at(0, 0) = 1.0;
    M.A_at(1, 1) = 1.0;
    M.B_at(0, 0) = 1.0;
    M.B_at(1, 0) = 0.0;
    M.C_at(0, 0) = 1.0;
    M.C_at(0, 1) = 0.0;
    M.D_at(0, 0) = 0.0;  // singular without η

    constexpr double eta = 1e-3;
    const Z got = computeBerezinian(M, eta);

    EXPECT_FALSE(std::isnan(got.real()));
    EXPECT_FALSE(std::isnan(got.imag()));
    EXPECT_TRUE(std::isfinite(got.real()));
    EXPECT_TRUE(std::isfinite(got.imag()));

    // D_reg = iη ⇒ X = C/D_reg = −i/η on the (0,0) channel;
    // Ã = [[1 − (−i/η), 0], [0, 1]] = [[1 + i/η, 0], [0, 1]]
    // sdet = det(Ã)/det(D_reg) = (1 + i/η) / (iη)
    const Z expected = Z{1.0, 1.0 / eta} / Z{0.0, eta};
    const double scale = std::max(1.0, std::abs(expected));
    EXPECT_NEAR(got.real(), expected.real(), 1e-12 * scale);
    EXPECT_NEAR(got.imag(), expected.imag(), 1e-12 * scale);
}
