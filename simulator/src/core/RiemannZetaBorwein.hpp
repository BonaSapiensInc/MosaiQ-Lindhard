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
#include "physics/Constants.hpp"

#include <cstddef>
#include <optional>

namespace mosaiq {

/// Policy for deterministic Borwein evaluation of ζ(s). Truncation order is pinned.
struct BorweinPolicy {
    std::size_t truncation_order{constants::default_borwein_truncation_order};
    double absolute_tolerance{constants::default_borwein_absolute_tolerance};
};

/// Evaluate ζ(s) for real s > 1 by Borwein's convergent alternating series.
/// Returns nullopt if s ≤ 1, non-finite, or truncation_order == 0.
[[nodiscard]] std::optional<double> riemann_zeta_borwein(double s,
                                                         BorweinPolicy policy = {});

/// Closed forms ζ(2k) for k = 1,2 (regression anchors). Returns 0 for unsupported k.
template<ScalarPhysical T = double>
[[nodiscard]] constexpr T riemann_zeta_even_integer(std::size_t k) noexcept
{
    if (k == 1) {
        return static_cast<T>(constants::pi * constants::pi / 6.0);
    }
    if (k == 2) {
        const T pi2 = static_cast<T>(constants::pi * constants::pi);
        return pi2 * pi2 / static_cast<T>(90);
    }
    return T{};
}

}  // namespace mosaiq
