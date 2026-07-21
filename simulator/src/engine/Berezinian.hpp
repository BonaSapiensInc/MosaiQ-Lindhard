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
///   2. Copy C_pe → X and solve D_reg X = C (LDLT preferred; LU fallback).
///   3. Ã = A_ee − B_ep X.
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
        const bool ok =
            linalg::solve_DX_eq_C<Z>(D_try, n_p, X_try, n_e, use_ldlt);
        if (ok) {
            out.X.swap(X_try);
        }
        return ok;
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

}  // namespace mosaiq
