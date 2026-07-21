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

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mosaiq {
namespace linalg {

/// Zero-dependency dense linear algebra for Berezinian Schur complements.
/// Doctrine: docs/BEREZINIAN_ARCHITECTURE.md §3.1 — never form D^{-1} explicitly;
/// solve D X = C then Ã = A − B X.

inline constexpr double k_pivot_epsilon = 1e-14;

template<typename T>
[[nodiscard]] constexpr T abs_like(T x) noexcept
{
    using std::abs;
    return abs(x);
}

template<typename T>
[[nodiscard]] constexpr T abs_like(std::complex<T> z) noexcept
{
    return std::abs(z);
}

/// Row-major dense matrix view (non-owning).
template<typename T>
struct MatrixView {
    std::span<T> data{};
    std::size_t rows{0};
    std::size_t cols{0};

    [[nodiscard]] T& operator()(std::size_t i, std::size_t j) noexcept
    {
        assert(i < rows && j < cols);
        return data[i * cols + j];
    }
    [[nodiscard]] const T& operator()(std::size_t i, std::size_t j) const noexcept
    {
        assert(i < rows && j < cols);
        return data[i * cols + j];
    }
};

template<typename T>
[[nodiscard]] constexpr MatrixView<T> make_view(std::span<T> data,
                                                std::size_t rows,
                                                std::size_t cols) noexcept
{
    assert(data.size() == rows * cols);
    return MatrixView<T>{data, rows, cols};
}

/// C ← α A B + β C  (all row-major). Dimensions: A(m×k), B(k×n), C(m×n).
template<typename T>
void gemm(std::size_t m,
          std::size_t n,
          std::size_t k,
          T alpha,
          std::span<const T> A,
          std::span<const T> B,
          T beta,
          std::span<T> C) noexcept
{
    assert(A.size() == m * k);
    assert(B.size() == k * n);
    assert(C.size() == m * n);

    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            T sum{};
            for (std::size_t p = 0; p < k; ++p) {
                sum += A[i * k + p] * B[p * n + j];
            }
            C[i * n + j] = alpha * sum + beta * C[i * n + j];
        }
    }
}

/// In-place LU with partial pivoting: A → P^{-1} L U (unit lower L).
/// `pivots[i]` = row exchanged with i at step i.
template<typename T>
[[nodiscard]] bool lu_factorize(MatrixView<T> A, std::span<std::size_t> pivots)
{
    const std::size_t n = A.rows;
    if (A.rows != A.cols || pivots.size() != n) {
        throw std::invalid_argument("lu_factorize: square matrix and pivot length required");
    }

    for (std::size_t i = 0; i < n; ++i) {
        pivots[i] = i;
    }

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        auto max_val = abs_like(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            const auto candidate = abs_like(A(i, k));
            if (candidate > max_val) {
                max_val = candidate;
                pivot = i;
            }
        }
        if (max_val < static_cast<decltype(max_val)>(k_pivot_epsilon)) {
            return false;
        }
        if (pivot != k) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(A(k, j), A(pivot, j));
            }
            std::swap(pivots[k], pivots[pivot]);
        }

        const T akk = A(k, k);
        for (std::size_t i = k + 1; i < n; ++i) {
            A(i, k) /= akk;
            const T lik = A(i, k);
            for (std::size_t j = k + 1; j < n; ++j) {
                A(i, j) -= lik * A(k, j);
            }
        }
    }
    return true;
}

/// Solve A X = B after `lu_factorize(A, pivots)`. B is n×nrhs, overwritten with X.
template<typename T>
void lu_solve(MatrixView<T> LU, std::span<const std::size_t> pivots, MatrixView<T> B)
{
    const std::size_t n = LU.rows;
    const std::size_t nrhs = B.cols;
    if (LU.rows != LU.cols || B.rows != n || pivots.size() != n) {
        throw std::invalid_argument("lu_solve: dimension mismatch");
    }

    // Apply the same row permutation recorded during factorization.
    std::vector<T> scratch(static_cast<std::size_t>(n * nrhs));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < nrhs; ++j) {
            scratch[i * nrhs + j] = B(pivots[i], j);
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < nrhs; ++j) {
            B(i, j) = scratch[i * nrhs + j];
        }
    }

    // Forward substitution: L Y = P B  (unit diagonal).
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < nrhs; ++j) {
            T sum = B(i, j);
            for (std::size_t k = 0; k < i; ++k) {
                sum -= LU(i, k) * B(k, j);
            }
            B(i, j) = sum;
        }
    }

    // Back substitution: U X = Y.
    for (std::size_t i = n; i-- > 0;) {
        for (std::size_t j = 0; j < nrhs; ++j) {
            T sum = B(i, j);
            for (std::size_t k = i + 1; k < n; ++k) {
                sum -= LU(i, k) * B(k, j);
            }
            B(i, j) = sum / LU(i, i);
        }
    }
}

/// In-place LDLT for real-symmetric / complex-symmetric matrices (no pivoting).
/// Overwrites A with L (unit lower, strict) below diagonal and D on diagonal.
/// Returns false if a diagonal pivot is singular.
template<typename T>
[[nodiscard]] bool ldlt_factorize(MatrixView<T> A)
{
    const std::size_t n = A.rows;
    if (A.rows != A.cols) {
        throw std::invalid_argument("ldlt_factorize: square matrix required");
    }

    for (std::size_t j = 0; j < n; ++j) {
        // d_j = a_jj − Σ_{k<j} l_jk² d_k
        T dj = A(j, j);
        for (std::size_t k = 0; k < j; ++k) {
            const T ljk = A(j, k);
            dj -= ljk * ljk * A(k, k);
        }
        if (abs_like(dj) < static_cast<decltype(abs_like(dj))>(k_pivot_epsilon)) {
            return false;
        }
        A(j, j) = dj;

        for (std::size_t i = j + 1; i < n; ++i) {
            T lij = A(i, j);
            for (std::size_t k = 0; k < j; ++k) {
                lij -= A(i, k) * A(j, k) * A(k, k);
            }
            A(i, j) = lij / dj;
            A(j, i) = A(i, j);  // keep symmetric mirror for convenience
        }
    }
    return true;
}

/// Solve A X = B after `ldlt_factorize(A)`. B is n×nrhs, overwritten with X.
template<typename T>
void ldlt_solve(MatrixView<T> LD, MatrixView<T> B)
{
    const std::size_t n = LD.rows;
    const std::size_t nrhs = B.cols;
    if (LD.rows != LD.cols || B.rows != n) {
        throw std::invalid_argument("ldlt_solve: dimension mismatch");
    }

    // L Y = B
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < nrhs; ++j) {
            T sum = B(i, j);
            for (std::size_t k = 0; k < i; ++k) {
                sum -= LD(i, k) * B(k, j);
            }
            B(i, j) = sum;
        }
    }

    // D Z = Y
    for (std::size_t i = 0; i < n; ++i) {
        const T di = LD(i, i);
        for (std::size_t j = 0; j < nrhs; ++j) {
            B(i, j) /= di;
        }
    }

    // Lᵀ X = Z
    for (std::size_t i = n; i-- > 0;) {
        for (std::size_t j = 0; j < nrhs; ++j) {
            T sum = B(i, j);
            for (std::size_t k = i + 1; k < n; ++k) {
                sum -= LD(k, i) * B(k, j);
            }
            B(i, j) = sum;
        }
    }
}

/// Solve D X = RHS for square D (n×n) and RHS (n×nrhs).
/// Prefers LDLT when `prefer_ldlt` is true; on LDLT failure returns false without
/// attempting LU on the mutated factor (caller must retry with a fresh D copy and
/// `prefer_ldlt == false`). When `prefer_ldlt` is false, uses LU with partial pivoting.
template<typename T>
[[nodiscard]] bool solve_DX_eq_C(std::span<T> D_workspace,
                                 std::size_t n,
                                 std::span<T> rhs_workspace,
                                 std::size_t nrhs,
                                 bool prefer_ldlt = true)
{
    if (D_workspace.size() != n * n || rhs_workspace.size() != n * nrhs) {
        throw std::invalid_argument("solve_DX_eq_C: size mismatch");
    }

    auto D = make_view(D_workspace, n, n);
    auto X = make_view(rhs_workspace, n, nrhs);

    if (prefer_ldlt) {
        if (!ldlt_factorize(D)) {
            return false;
        }
        ldlt_solve(D, X);
        return true;
    }

    std::vector<std::size_t> pivots(n);
    if (!lu_factorize(D, pivots)) {
        return false;
    }
    lu_solve(D, pivots, X);
    return true;
}

/// Ã = A − B X  with A,Ã (m×m), B (m×n), X (n×m).
template<typename T>
void schur_from_BX(std::size_t m,
                   std::size_t n,
                   std::span<const T> A,
                   std::span<const T> B,
                   std::span<const T> X,
                   std::span<T> A_tilde)
{
    if (A.size() != m * m || B.size() != m * n || X.size() != n * m ||
        A_tilde.size() != m * m) {
        throw std::invalid_argument("schur_from_BX: size mismatch");
    }
    for (std::size_t i = 0; i < m * m; ++i) {
        A_tilde[i] = A[i];
    }
    gemm(m, m, n, T{-1}, B, X, T{1}, A_tilde);
}

}  // namespace linalg
}  // namespace mosaiq
