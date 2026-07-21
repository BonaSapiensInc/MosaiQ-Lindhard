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
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace mosaiq {

/// Contiguous electron–phonon graded block layout (zero-copy `std::span` views):
///
///   [ A_ee (N_e×N_e) | B_ep (N_e×N_p) | C_pe (N_p×N_e) | D_pp (N_p×N_p) ]
///
/// Storage is row-major within each block. N_e and N_p may differ.
/// Element type `T` may be real or `std::complex<Real>` (Keldysh / iη path).
/// Doctrine: docs/BEREZINIAN_ARCHITECTURE.md §3.3.
template<typename T = double>
class ElectronPhononSuperMatrix {
public:
    ElectronPhononSuperMatrix() = default;

    ElectronPhononSuperMatrix(std::size_t n_e, std::size_t n_p)
        : n_e_{n_e}
        , n_p_{n_p}
        , storage_(total_entries(n_e, n_p), T{})
    {
        if (n_e == 0 || n_p == 0) {
            throw std::invalid_argument(
                "ElectronPhononSuperMatrix: N_e and N_p must be positive");
        }
    }

    [[nodiscard]] std::size_t n_e() const noexcept { return n_e_; }
    [[nodiscard]] std::size_t n_p() const noexcept { return n_p_; }

    [[nodiscard]] static constexpr std::size_t total_entries(std::size_t n_e,
                                                            std::size_t n_p) noexcept
    {
        return n_e * n_e + n_e * n_p + n_p * n_e + n_p * n_p;
    }

    /// Fermion block A_ee — contiguous view, no deep copy.
    [[nodiscard]] std::span<T> A() noexcept
    {
        return {storage_.data() + offset_A(), n_e_ * n_e_};
    }
    [[nodiscard]] std::span<const T> A() const noexcept
    {
        return {storage_.data() + offset_A(), n_e_ * n_e_};
    }

    /// Electron→phonon cross block B_ep.
    [[nodiscard]] std::span<T> B() noexcept
    {
        return {storage_.data() + offset_B(), n_e_ * n_p_};
    }
    [[nodiscard]] std::span<const T> B() const noexcept
    {
        return {storage_.data() + offset_B(), n_e_ * n_p_};
    }

    /// Phonon→electron cross block C_pe.
    [[nodiscard]] std::span<T> C() noexcept
    {
        return {storage_.data() + offset_C(), n_p_ * n_e_};
    }
    [[nodiscard]] std::span<const T> C() const noexcept
    {
        return {storage_.data() + offset_C(), n_p_ * n_e_};
    }

    /// Boson / phonon dynamical block D_pp.
    [[nodiscard]] std::span<T> D() noexcept
    {
        return {storage_.data() + offset_D(), n_p_ * n_p_};
    }
    [[nodiscard]] std::span<const T> D() const noexcept
    {
        return {storage_.data() + offset_D(), n_p_ * n_p_};
    }

    /// Full contiguous backing store (A|B|C|D).
    [[nodiscard]] std::span<T> data() noexcept { return storage_; }
    [[nodiscard]] std::span<const T> data() const noexcept { return storage_; }

    /// Element access within a block (row-major).
    [[nodiscard]] T& A_at(std::size_t i, std::size_t j) noexcept
    {
        assert(i < n_e_ && j < n_e_);
        return storage_[offset_A() + i * n_e_ + j];
    }
    [[nodiscard]] const T& A_at(std::size_t i, std::size_t j) const noexcept
    {
        assert(i < n_e_ && j < n_e_);
        return storage_[offset_A() + i * n_e_ + j];
    }

    [[nodiscard]] T& B_at(std::size_t i, std::size_t j) noexcept
    {
        assert(i < n_e_ && j < n_p_);
        return storage_[offset_B() + i * n_p_ + j];
    }
    [[nodiscard]] const T& B_at(std::size_t i, std::size_t j) const noexcept
    {
        assert(i < n_e_ && j < n_p_);
        return storage_[offset_B() + i * n_p_ + j];
    }

    [[nodiscard]] T& C_at(std::size_t i, std::size_t j) noexcept
    {
        assert(i < n_p_ && j < n_e_);
        return storage_[offset_C() + i * n_e_ + j];
    }
    [[nodiscard]] const T& C_at(std::size_t i, std::size_t j) const noexcept
    {
        assert(i < n_p_ && j < n_e_);
        return storage_[offset_C() + i * n_e_ + j];
    }

    [[nodiscard]] T& D_at(std::size_t i, std::size_t j) noexcept
    {
        assert(i < n_p_ && j < n_p_);
        return storage_[offset_D() + i * n_p_ + j];
    }
    [[nodiscard]] const T& D_at(std::size_t i, std::size_t j) const noexcept
    {
        assert(i < n_p_ && j < n_p_);
        return storage_[offset_D() + i * n_p_ + j];
    }

private:
    [[nodiscard]] constexpr std::size_t offset_A() const noexcept { return 0; }
    [[nodiscard]] constexpr std::size_t offset_B() const noexcept { return n_e_ * n_e_; }
    [[nodiscard]] constexpr std::size_t offset_C() const noexcept
    {
        return n_e_ * n_e_ + n_e_ * n_p_;
    }
    [[nodiscard]] constexpr std::size_t offset_D() const noexcept
    {
        return n_e_ * n_e_ + n_e_ * n_p_ + n_p_ * n_e_;
    }

    std::size_t n_e_{0};
    std::size_t n_p_{0};
    std::vector<T> storage_{};
};

}  // namespace mosaiq
