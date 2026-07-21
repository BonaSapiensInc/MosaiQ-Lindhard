/* ==========================================================================
 * MOSAIQ-LINDHARD ENGINE
 * Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
 * * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
 * Free for non-commercial personal and academic research use only.
 * Commercial, governmental, and public institutional use requires a
 * separate paid license. See LICENSE file for details.
 * * Contact: kim.ingee@bonasapiens.com
 * ========================================================================== */

#pragma once

#include "core/Concepts.hpp"
#include "core/LinalgUtils.hpp"
#include "core/SuperMatrix.hpp"

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace mosaiq {

/// Schur complement Ã = A − B X with D_reg X = C, D_reg = D + iη I.
/// Doctrine: docs/BEREZINIAN_ARCHITECTURE.md §2–§3.
template<ScalarPhysical Real = double>
struct SchurComplementResult {
    std::size_t n_e{0};
    std::size_t n_p{0};
    Real eta{0};
    /// det(D_reg) extracted from the successful LU/LDLT factor (no extra O(N³)).
    std::complex<Real> det_D{1};
    /// Row-major Ã (N_e × N_e).
    std::vector<std::complex<Real>> A_tilde{};
    /// Row-major X solving D_reg X = C (N_p × N_e); retained for audits.
    std::vector<std::complex<Real>> X{};

    [[nodiscard]] std::span<const std::complex<Real>> schur() const noexcept
    {
        return A_tilde;
    }
};

/// Apply retarded regularization D → D + iη I on a square row-major block.
template<ScalarPhysical Real>
void apply_retarded_eta(std::span<std::complex<Real>> D, std::size_t n, Real eta) noexcept
{
    assert(D.size() == n * n);
    const std::complex<Real> i_eta{Real{0}, eta};
    for (std::size_t i = 0; i < n; ++i) {
        D[i * n + i] += i_eta;
    }
}

/// Compute the fermion-block Schur complement via the DX = C strategy.
///
/// Steps:
///   1. Copy D_pp → D_reg and load iη on the diagonal (causality).
///   2. Factor D_reg (LDLT preferred; LU fallback); harvest det(D_reg) from the factor.
///   3. Solve D_reg X = C on the same factor; Ã = A_ee − B_ep X.
///
/// Never forms D^{-1} explicitly.
template<ScalarPhysical Real = double>
[[nodiscard]] SchurComplementResult<Real> computeSchurComplement(
    const ElectronPhononSuperMatrix<std::complex<Real>>& M,
    Real eta,
    bool prefer_ldlt = true)
{
    if (eta < Real{0}) {
        throw std::invalid_argument("computeSchurComplement: eta must be non-negative");
    }

    const std::size_t n_e = M.n_e();
    const std::size_t n_p = M.n_p();

    SchurComplementResult<Real> out;
    out.n_e = n_e;
    out.n_p = n_p;
    out.eta = eta;
    out.A_tilde.assign(n_e * n_e, std::complex<Real>{});
    out.X.assign(n_p * n_e, std::complex<Real>{});

    using Z = std::complex<Real>;

    // --- workspace copies (factorization mutates D; RHS becomes X) ---
    std::vector<Z> D_work(M.D().begin(), M.D().end());
    apply_retarded_eta<Real>(D_work, n_p, eta);

    std::copy(M.C().begin(), M.C().end(), out.X.begin());

    auto try_solve = [&](bool use_ldlt) -> bool {
        std::vector<Z> D_try = D_work;
        std::vector<Z> X_try = out.X;
        auto D_view = linalg::make_view(std::span<Z>{D_try}, n_p, n_p);
        auto X_view = linalg::make_view(std::span<Z>{X_try}, n_p, n_e);

        if (use_ldlt) {
            if (!linalg::ldlt_factorize(D_view)) {
                return false;
            }
            out.det_D = linalg::ldlt_determinant(D_view);
            linalg::ldlt_solve(D_view, X_view);
        } else {
            std::vector<std::size_t> pivots(n_p);
            if (!linalg::lu_factorize(D_view, pivots)) {
                return false;
            }
            out.det_D = linalg::lu_determinant(D_view, pivots);
            linalg::lu_solve(D_view, pivots, X_view);
        }

        out.X.swap(X_try);
        return true;
    };

    bool solved = false;
    if (prefer_ldlt) {
        solved = try_solve(true);
    }
    if (!solved) {
        solved = try_solve(false);
    }
    if (!solved) {
        throw std::runtime_error(
            "computeSchurComplement: D_reg is singular after iη regularization");
    }

    linalg::schur_from_BX<Z>(n_e, n_p, M.A(), M.B(), out.X, out.A_tilde);
    return out;
}

/// Convenience overload: promote a real-valued super-matrix to complex, then Schur.
template<ScalarPhysical Real = double>
[[nodiscard]] SchurComplementResult<Real> computeSchurComplement(
    const ElectronPhononSuperMatrix<Real>& M_real,
    Real eta,
    bool prefer_ldlt = true)
{
    ElectronPhononSuperMatrix<std::complex<Real>> M(M_real.n_e(), M_real.n_p());
    auto promote = [](std::span<const Real> src, std::span<std::complex<Real>> dst) {
        for (std::size_t i = 0; i < src.size(); ++i) {
            dst[i] = std::complex<Real>{src[i], Real{0}};
        }
    };
    promote(M_real.A(), M.A());
    promote(M_real.B(), M.B());
    promote(M_real.C(), M.C());
    promote(M_real.D(), M.D());
    return computeSchurComplement<Real>(M, eta, prefer_ldlt);
}

/// Berezinian (superdeterminant): sdet(M) = det(Ã) / det(D_reg), Ã = A − B D_reg^{-1} C.
template<ScalarPhysical Real = double>
[[nodiscard]] std::complex<Real> computeBerezinian(
    const ElectronPhononSuperMatrix<std::complex<Real>>& M,
    Real eta)
{
    using Z = std::complex<Real>;

    const auto schur = computeSchurComplement<Real>(M, eta);
    if (schur.det_D == Z{}) {
        throw std::runtime_error("computeBerezinian: det(D_reg) vanished");
    }

    std::vector<Z> A_work = schur.A_tilde;
    auto A_view = linalg::make_view(std::span<Z>{A_work}, schur.n_e, schur.n_e);
    std::vector<std::size_t> pivots(schur.n_e);
    if (!linalg::lu_factorize(A_view, pivots)) {
        throw std::runtime_error("computeBerezinian: Schur block Ã is singular");
    }

    const Z det_A = linalg::lu_determinant(A_view, pivots);
    return det_A / schur.det_D;
}

}  // namespace mosaiq
