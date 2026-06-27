#!/usr/bin/env python3

# ==========================================================================
# MOSAIQ-LINDHARD ENGINE
# Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
# * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
# Free for non-commercial personal and academic research use only.
# Commercial, governmental, and public institutional use requires a
# separate paid license. See LICENSE file for details.
# * Contact: kim.ingee@bonasapiens.com
# ==========================================================================

"""Figure 7: thermal anatomy of bare electron Lindhard blurring (2x2 companion panel)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, TwoSlopeNorm

from plot_common import OUTPUT_DIR, save_figure
from plot_lindhard_t0_2d import calc_t0_lindhard
from plot_Sqw import (
    CONTOUR_ISOLINE_ALPHA,
    CONTOUR_ISOLINE_COLOR,
    CONTOUR_ISOLINE_LINEWIDTH,
    ELECTRON_DISPLAY_BOUNDS,
    build_contour_grid,
    configure_seaborn,
    filter_contour_line_levels,
    format_contour_display_axes,
    frame_contour_axis,
    interpolate_to_display_mesh,
    load_gridded,
    log_contour_levels,
    prepare_plot_grid,
    rainbow_log_cmap,
)

RS = 1.0
Q_SLICE = 1.0
ELECTRON_Q_MAX = 4.0
IM_CHI_E_COL = 3
S_EE_COL = 3

TEMPERATURES_K: tuple[int, ...] = (11, 1000, 10000)
LADDER_PANELS: tuple[tuple[int, str], ...] = (
    (1000, "(a)"),
    (11, "(b)"),
)
SLICE_TEMPERATURES_K: tuple[int, ...] = (11, 1000, 10000)

SIMULATOR = Path(__file__).resolve().parent.parent / "simulator" / "build" / "mosaiq_simulator"


def lindhard_dat_path(t_kelvin: int) -> Path:
    if t_kelvin == 10000:
        primary = OUTPUT_DIR / "output_lindhard_base.dat"
        tagged = OUTPUT_DIR / "output_lindhard_base_T10000.dat"
        return tagged if tagged.is_file() else primary
    return OUTPUT_DIR / f"output_lindhard_base_T{t_kelvin}.dat"


def structure_factor_dat_path() -> Path:
    tagged = OUTPUT_DIR / "output_structure_factor_T10000.dat"
    primary = OUTPUT_DIR / "output_structure_factor.dat"
    return tagged if tagged.is_file() else primary


def ensure_lindhard_data(t_kelvin: int) -> Path:
    path = lindhard_dat_path(t_kelvin)
    if path.is_file():
        return path
    if not SIMULATOR.is_file():
        raise SystemExit(
            f"Missing {path} and simulator not built. Run:\n"
            "  cmake -S simulator -B simulator/build -DCMAKE_BUILD_TYPE=Release\n"
            "  cmake --build simulator/build --target mosaiq_simulator\n"
            f"  {SIMULATOR} {RS} {t_kelvin}"
        )
    subprocess.run([str(SIMULATOR), str(RS), str(float(t_kelvin))], check=True)
    generated = OUTPUT_DIR / "output_lindhard_base.dat"
    if t_kelvin != 10000:
        path.write_bytes(generated.read_bytes())
    else:
        path = generated
    return path


def load_im_chi_electron(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    rows = data[data[:, 0] <= ELECTRON_Q_MAX + 1.0e-9]
    q_lo, q_hi, w_lo, w_hi = ELECTRON_DISPLAY_BOUNDS
    q_grid, w_grid, value_grid = build_contour_grid(
        rows,
        IM_CHI_E_COL,
        q_min=max(q_lo, 0.01),
        q_max=q_hi,
        w_max=w_hi,
        n_omega=400,
        w_min=max(w_lo, 0.01),
    )
    magnitude = np.abs(value_grid)
    return interpolate_to_display_mesh(q_grid, w_grid, magnitude, q_lo, q_hi, w_lo, w_hi)


def spe_boundaries(q_bar: float) -> tuple[float, float]:
    """Analytic T=0 single-particle edges at fixed q (positive omega axis)."""
    lower = abs(q_bar * q_bar - 2.0 * q_bar)
    upper = q_bar * q_bar + 2.0 * q_bar
    return lower, upper


def nearest_q_index(q_axis: np.ndarray, q_target: float) -> int:
    return int(np.argmin(np.abs(q_axis - q_target)))


def slice_at_q(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    grid: np.ndarray,
    q_target: float,
) -> tuple[np.ndarray, np.ndarray]:
    q_axis = q_grid[0, :]
    w_axis = w_grid[:, 0]
    idx = nearest_q_index(q_axis, q_target)
    values = np.abs(grid[:, idx])
    return w_axis, values


def draw_log_contour_panel(
    ax: plt.Axes,
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    z_grid: np.ndarray,
    *,
    title: str,
    shared_norm: LogNorm | None = None,
) -> LogNorm | None:
    plot_grid, norm = prepare_plot_grid(z_grid)
    color_norm = shared_norm if shared_norm is not None else norm
    if color_norm is None:
        ax.text(0.5, 0.5, "No positive spectral weight", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title, fontweight="bold")
        return None

    cmap = rainbow_log_cmap()
    fill_levels, line_levels = log_contour_levels(color_norm, n_filled=200, n_lines=72)
    line_levels = filter_contour_line_levels(line_levels, color_norm)

    ax.contourf(
        q_grid,
        w_grid,
        plot_grid,
        levels=fill_levels,
        cmap=cmap,
        norm=color_norm,
        extend="min",
    )
    ax.contour(
        q_grid,
        w_grid,
        plot_grid,
        levels=line_levels,
        norm=color_norm,
        colors=CONTOUR_ISOLINE_COLOR,
        linewidths=CONTOUR_ISOLINE_LINEWIDTH,
        alpha=CONTOUR_ISOLINE_ALPHA,
    )
    frame_contour_axis(ax)
    format_contour_display_axes(ax, *ELECTRON_DISPLAY_BOUNDS)
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel(r"$\bar{q}$")
    ax.set_ylabel(r"$\bar{\omega}$")
    return color_norm


def ladder_shared_norm(temperatures: tuple[int, ...]) -> LogNorm:
    maxima: list[float] = []
    minima: list[float] = []
    for t_kelvin in temperatures:
        path = ensure_lindhard_data(t_kelvin)
        _, _, z_grid = load_im_chi_electron(path)
        positive = z_grid[np.isfinite(z_grid) & (z_grid > 0.0)]
        if positive.size:
            maxima.append(float(np.max(positive)))
            minima.append(float(np.min(positive)))
    if not maxima:
        raise SystemExit("No positive bare Im chi_e data for thermal ladder panels.")
    vmax = max(maxima)
    vmin = max(min(minima), vmax * 1.0e-6, 1.0e-12)
    return LogNorm(vmin=vmin, vmax=vmax)


def render_slice_panel(ax: plt.Axes) -> None:
    w_lo, w_hi = ELECTRON_DISPLAY_BOUNDS[2], ELECTRON_DISPLAY_BOUNDS[3]
    w_fine = np.linspace(max(w_lo, 0.01), w_hi, 600)
    q_line = np.full_like(w_fine, Q_SLICE)
    _, im_t0 = calc_t0_lindhard(q_line, w_fine)
    ax.plot(
        w_fine,
        np.abs(im_t0),
        color="black",
        linewidth=2.0,
        label=r"$T=0$ analytic",
        zorder=5,
    )

    colors = {11: "#1f3b73", 1000: "#8b4513", 10000: "#b22222"}
    for t_kelvin in SLICE_TEMPERATURES_K:
        path = ensure_lindhard_data(t_kelvin)
        q_grid, w_grid, z_grid = load_im_chi_electron(path)
        w_axis, values = slice_at_q(q_grid, w_grid, z_grid, Q_SLICE)
        ax.plot(
            w_axis,
            values,
            color=colors[t_kelvin],
            linewidth=1.6,
            label=rf"$T={t_kelvin:,}\,\mathrm{{K}}$",
        )

    omega_lo, omega_hi = spe_boundaries(Q_SLICE)
    ax.axvline(omega_lo, color="black", linestyle="--", linewidth=1.2, alpha=0.75)
    ax.axvline(omega_hi, color="black", linestyle="--", linewidth=1.2, alpha=0.75)

    ax.set_yscale("log")
    ax.set_xlim(w_lo, w_hi)
    ax.set_xlabel(r"$\bar{\omega}$ at $\bar{q}=1.0$")
    ax.set_ylabel(r"$|\Im \tilde{\chi}_e^L|$")
    ax.set_title(r"(c) Thermal Fermi-step broadening at $\bar{q}=1.0$", fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ax.grid(True, alpha=0.25, linewidth=0.6)

    y_top = ax.get_ylim()[1]
    ax.text(omega_lo + 0.04, y_top * 0.55, rf"$\bar{{\omega}}={omega_lo:.1f}$", fontsize=9, color="black")
    ax.text(omega_hi - 0.55, y_top * 0.55, rf"$\bar{{\omega}}={omega_hi:.1f}$", fontsize=9, color="black")


def render_bare_rpa_panel(ax: plt.Axes) -> None:
    chi_path = ensure_lindhard_data(10000)
    sf_path = structure_factor_dat_path()
    if not sf_path.is_file():
        raise SystemExit(
            f"Missing {sf_path}. Run: {SIMULATOR} {RS} 10000 "
            "and save output/output_structure_factor_T10000.dat"
        )

    q_chi, w_chi, chi_grid = load_im_chi_electron(chi_path)
    log_bare = np.log10(np.maximum(chi_grid, 1.0e-30))

    sf_data = np.loadtxt(sf_path, comments="#")
    sf_rows = sf_data[sf_data[:, 0] <= ELECTRON_Q_MAX + 1.0e-9]
    q_lo, q_hi, w_lo, w_hi = ELECTRON_DISPLAY_BOUNDS
    q_sf, w_sf, s_grid = build_contour_grid(
        sf_rows,
        S_EE_COL,
        q_min=max(q_lo, 0.01),
        q_max=q_hi,
        w_max=w_hi,
        n_omega=400,
        w_min=max(w_lo, 0.01),
    )
    _, _, s_display = interpolate_to_display_mesh(
        q_sf, w_sf, np.abs(s_grid), q_lo, q_hi, w_lo, w_hi
    )
    log_rpa = np.log10(np.maximum(s_display, 1.0e-30))
    contrast = log_rpa - log_bare

    finite = contrast[np.isfinite(contrast)]
    if finite.size == 0:
        raise SystemExit("Bare-vs-RPA contrast map is empty.")

    vmax = float(np.percentile(finite, 98.0))
    vmin = float(np.percentile(finite, 2.0))
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=max(vmax, 0.5))

    cmap = plt.get_cmap("RdBu_r").copy()
    frame_contour_axis(ax)
    mesh = ax.pcolormesh(q_chi, w_chi, contrast, cmap=cmap, norm=norm, shading="auto")
    format_contour_display_axes(ax, *ELECTRON_DISPLAY_BOUNDS)
    ax.set_xlabel(r"$\bar{q}$")
    ax.set_ylabel(r"$\bar{\omega}$")
    ax.set_title(
        r"(d) $\log_{10}|S_{ee}| - \log_{10}|\Im \tilde{\chi}_e^L|$ at $T=10{,}000\,\mathrm{K}$",
        fontweight="bold",
    )
    cbar = ax.figure.colorbar(mesh, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"Interaction contrast [decades]")


def main() -> None:
    configure_seaborn()
    shared_norm = ladder_shared_norm(tuple(t for t, _ in LADDER_PANELS))

    fig, axes = plt.subplots(2, 2, figsize=(14, 11), layout="constrained")

    for ax, (t_kelvin, panel_tag) in zip(axes.flat[:2], LADDER_PANELS, strict=True):
        path = ensure_lindhard_data(t_kelvin)
        q_grid, w_grid, z_grid = load_im_chi_electron(path)
        draw_log_contour_panel(
            ax,
            q_grid,
            w_grid,
            z_grid,
            title=rf"{panel_tag} $|\Im \tilde{{\chi}}_e^L|$, $T={t_kelvin:,}\,\mathrm{{K}}$",
            shared_norm=shared_norm,
        )

    render_slice_panel(axes[1, 0])
    render_bare_rpa_panel(axes[1, 1])

    fig.suptitle(
        rf"Thermal anatomy of bare electron Lindhard blurring ($r_s={RS:.1f}$)",
        fontsize=14,
        fontweight="bold",
    )

    saved = save_figure(fig, "thermal_anatomy")
    plt.close(fig)
    print(f"Saved {saved}")


if __name__ == "__main__":
    main()
