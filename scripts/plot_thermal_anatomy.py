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

"""Figure 7: thermal anatomy of bare electron Lindhard real-part sharpening (four panel PDFs)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm

from plot_common import OUTPUT_DIR, save_figure
from plot_lindhard_t0_2d import calc_t0_lindhard
from plot_Sqw import (
    ELECTRON_DISPLAY_BOUNDS,
    build_contour_grid,
    configure_seaborn,
    draw_contour_isolines,
    fill_display_border,
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
RE_CHI_E_COL = 2
S_EE_COL = 3

COLORBAR_FRACTION = 0.046
COLORBAR_PAD = 0.04
LADDER_CBAR_LABEL = r"$|\Re \tilde{\chi}_e^L|$"
RATIO_CBAR_LABEL = r"$S_{ee}\,/\,|\Re \tilde{\chi}_e^L|$"
SINGLE_PANEL_FIGSIZE = (10, 7)
COLD_LADDER_T_K = 70
COLD_LADDER_Q_POINTS = 79

LADDER_PANELS: tuple[tuple[int, str], ...] = (
    (1000, "thermal_anatomy_a"),
    (COLD_LADDER_T_K, "thermal_anatomy_b"),
)

PANEL_C_STEM = "thermal_anatomy_c"
PANEL_D_STEM = "thermal_anatomy_d"

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


def load_re_chi_electron(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Match Figure 3 pipeline; fall back to contour-grid interpolation for sparse exports."""
    data = np.loadtxt(path, comments="#")
    rows = data[data[:, 0] <= ELECTRON_Q_MAX + 1.0e-9]
    q_lo, q_hi, w_lo, w_hi = ELECTRON_DISPLAY_BOUNDS
    try:
        q_grid, w_grid, value_grid = load_gridded(rows, RE_CHI_E_COL)
    except ValueError:
        q_data_min = float(np.min(rows[:, 0]))
        positive_w = rows[:, 1][rows[:, 1] > 0.0]
        w_data_min = float(np.min(positive_w)) if positive_w.size else 1.0e-2
        q_grid, w_grid, value_grid = build_contour_grid(
            rows,
            RE_CHI_E_COL,
            q_min=0.0,
            q_max=q_hi,
            w_max=w_hi,
            n_omega=400,
            w_min=w_data_min,
            n_q=COLD_LADDER_Q_POINTS,
        )
    magnitude = np.abs(value_grid)
    return interpolate_to_display_mesh(q_grid, w_grid, magnitude, q_lo, q_hi, w_lo, w_hi)


def spe_boundaries(q_bar: float) -> tuple[float, float]:
    """Analytic T=0 single-particle edges at fixed q (positive omega axis)."""
    lower = abs(q_bar * q_bar - 2.0 * q_bar)
    upper = q_bar * q_bar + 2.0 * q_bar
    return lower, upper


def slice_re_chi_at_q(path: Path, q_target: float) -> tuple[np.ndarray, np.ndarray]:
    """Extract |Re chi_e| along fixed q from raw simulator samples (2D q–omega interpolation)."""
    rows = np.loadtxt(path, comments="#")
    rows = rows[rows[:, 0] <= ELECTRON_Q_MAX + 1.0e-9]
    q_pts = rows[:, 0]
    w_pts = rows[:, 1]
    v_pts = np.abs(rows[:, RE_CHI_E_COL])

    w_lo, w_hi = ELECTRON_DISPLAY_BOUNDS[2], ELECTRON_DISPLAY_BOUNDS[3]
    w_unique = np.unique(w_pts)
    if w_unique.size >= 2:
        dw = float(np.median(np.diff(np.sort(w_unique))))
    else:
        dw = 0.02
    half_band = max(dw * 0.75, 0.01)

    w_fine = np.linspace(max(w_lo, 0.01), w_hi, 600)
    values = np.full_like(w_fine, np.nan)
    for index, omega in enumerate(w_fine):
        band = np.abs(w_pts - omega) <= half_band
        if not np.any(band):
            continue
        q_band = q_pts[band]
        v_band = v_pts[band]
        order = np.argsort(q_band)
        q_sorted = q_band[order]
        v_sorted = v_band[order]
        q_unique, unique_idx = np.unique(q_sorted, return_index=True)
        v_unique = v_sorted[unique_idx]
        if q_target < q_unique[0] - 1.0e-9 or q_target > q_unique[-1] + 1.0e-9:
            continue
        if q_unique.size == 1:
            if abs(q_unique[0] - q_target) < 0.06:
                values[index] = v_unique[0]
            continue
        values[index] = np.interp(q_target, q_unique, v_unique)
    return w_fine, values


def plot_log_segments(
    ax: plt.Axes,
    w_axis: np.ndarray,
    values: np.ndarray,
    *,
    zorder: int = 3,
    **plot_kwargs,
) -> None:
    """Plot only finite, positive samples; break lines at NaN gaps for log axes."""
    mask = np.isfinite(values) & (values > 0.0)
    if not np.any(mask):
        return

    start: int | None = None
    for index, active in enumerate(mask):
        if active and start is None:
            start = index
        elif not active and start is not None:
            ax.plot(w_axis[start:index], values[start:index], zorder=zorder, **plot_kwargs)
            start = None
    if start is not None:
        ax.plot(w_axis[start:], values[start:], zorder=zorder, **plot_kwargs)


def render_log_rainbow_field(
    ax: plt.Axes,
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    z_grid: np.ndarray,
    color_norm: LogNorm,
    *,
    title: str,
    draw_isolines: bool = True,
    fill_border: bool = False,
) -> ScalarMappable:
    """Edge-to-edge log-rainbow field (contourf), matching Figure 3 styling."""
    field = fill_display_border(z_grid) if fill_border else z_grid
    plot_grid, _ = prepare_plot_grid(field)
    q_lo, q_hi, w_lo, w_hi = ELECTRON_DISPLAY_BOUNDS
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
    if draw_isolines:
        draw_contour_isolines(
            ax,
            q_grid,
            w_grid,
            plot_grid,
            color_norm,
            line_levels,
        )
    frame_contour_axis(ax)
    format_contour_display_axes(ax, q_lo, q_hi, w_lo, w_hi)
    ax.set_aspect("auto")
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel(r"$\bar{q}$")
    ax.set_ylabel(r"$\bar{\omega}$")
    scalar_map = ScalarMappable(cmap=cmap, norm=color_norm)
    scalar_map.set_array([])
    return scalar_map


def attach_log_colorbar(
    fig: plt.Figure,
    axes: plt.Axes | np.ndarray,
    mappable: ScalarMappable,
    label: str,
) -> plt.colorbar.Colorbar:
    return fig.colorbar(
        mappable,
        ax=axes,
        fraction=COLORBAR_FRACTION,
        pad=COLORBAR_PAD,
        label=label,
    )


def ladder_shared_norm(temperatures: tuple[int, ...]) -> LogNorm:
    maxima: list[float] = []
    minima: list[float] = []
    for t_kelvin in temperatures:
        path = ensure_lindhard_data(t_kelvin)
        _, _, z_grid = load_re_chi_electron(path)
        positive = z_grid[np.isfinite(z_grid) & (z_grid > 0.0)]
        if positive.size:
            maxima.append(float(np.max(positive)))
            minima.append(float(np.min(positive)))
    if not maxima:
        raise SystemExit("No positive bare Re chi_e data for thermal ladder panels.")
    vmax = max(maxima)
    vmin = max(min(minima), vmax * 1.0e-6, 1.0e-12)
    return LogNorm(vmin=vmin, vmax=vmax)


def render_slice_panel(ax: plt.Axes) -> None:
    w_lo, w_hi = ELECTRON_DISPLAY_BOUNDS[2], ELECTRON_DISPLAY_BOUNDS[3]
    w_fine = np.linspace(max(w_lo, 0.01), w_hi, 600)
    q_line = np.full_like(w_fine, Q_SLICE)
    re_t0, _ = calc_t0_lindhard(q_line, w_fine)
    omega_lo, omega_hi = spe_boundaries(Q_SLICE)
    ax.axvline(omega_lo, color="0.35", linestyle=":", linewidth=1.0, alpha=0.7, zorder=1)
    ax.axvline(omega_hi, color="0.35", linestyle=":", linewidth=1.0, alpha=0.7, zorder=1)

    finite_styles: tuple[tuple[int, str, str, float], ...] = (
        (10000, "#b22222", "dashdot", 1.7),
        (1000, "#8b4513", "dashed", 1.7),
        (COLD_LADDER_T_K, "#1f3b73", "solid", 1.8),
    )
    for t_kelvin, color, linestyle, linewidth in finite_styles:
        path = ensure_lindhard_data(t_kelvin)
        w_axis, values = slice_re_chi_at_q(path, Q_SLICE)
        plot_log_segments(
            ax,
            w_axis,
            values,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=rf"$T={t_kelvin:,}\,\mathrm{{K}}$",
            zorder=4,
        )

    plot_log_segments(
        ax,
        w_fine,
        np.abs(re_t0),
        color="black",
        linestyle="solid",
        linewidth=2.2,
        label=r"$T=0$ analytic",
        zorder=6,
    )

    ax.set_yscale("log")
    ax.set_xlim(w_lo, w_hi)
    ax.margins(x=0)
    ax.set_xlabel(r"$\bar{\omega}$ at $\bar{q}=1.0$")
    ax.set_ylabel(r"$|\Re \tilde{\chi}_e^L|$")
    ax.set_title(r"Real-part sharpening at $\bar{q}=1.0$", fontweight="bold")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.95, edgecolor="0.75")
    ax.grid(True, alpha=0.25, linewidth=0.6)


def load_rpa_bare_ratio() -> tuple[np.ndarray, np.ndarray, np.ndarray, LogNorm]:
    chi_path = ensure_lindhard_data(10000)
    sf_path = structure_factor_dat_path()
    if not sf_path.is_file():
        raise SystemExit(
            f"Missing {sf_path}. Run: {SIMULATOR} {RS} 10000 "
            "and save output/output_structure_factor_T10000.dat"
        )

    q_grid, w_grid, chi_grid = load_re_chi_electron(chi_path)

    sf_data = np.loadtxt(sf_path, comments="#")
    sf_rows = sf_data[sf_data[:, 0] <= ELECTRON_Q_MAX + 1.0e-9]
    q_lo, q_hi, w_lo, w_hi = ELECTRON_DISPLAY_BOUNDS
    q_sf, w_sf, s_grid = load_gridded(sf_rows, S_EE_COL)
    _, _, s_display = interpolate_to_display_mesh(
        q_sf, w_sf, np.abs(s_grid), q_lo, q_hi, w_lo, w_hi
    )

    chi_floor = np.maximum(chi_grid, 1.0e-30)
    ratio = np.where(
        np.isfinite(s_display) & np.isfinite(chi_grid) & (chi_grid > 0.0),
        s_display / chi_floor,
        np.nan,
    )
    _, norm = prepare_plot_grid(ratio)
    if norm is None:
        raise SystemExit("Bare-vs-RPA ratio map is empty.")
    return q_grid, w_grid, ratio, norm


def render_bare_rpa_panel(ax: plt.Axes) -> ScalarMappable:
    q_grid, w_grid, ratio, norm = load_rpa_bare_ratio()
    return render_log_rainbow_field(
        ax,
        q_grid,
        w_grid,
        ratio,
        norm,
        title=(
            r"RPA/bare ratio $S_{ee}/|\Re \tilde{\chi}_e^L|$"
            r" at $T=10{,}000\,\mathrm{K}$"
        ),
        draw_isolines=True,
    )


def save_ladder_panel(
    t_kelvin: int,
    stem: str,
    shared_norm: LogNorm,
) -> Path:
    path = ensure_lindhard_data(t_kelvin)
    q_grid, w_grid, z_grid = load_re_chi_electron(path)
    fig, ax = plt.subplots(figsize=SINGLE_PANEL_FIGSIZE, layout="constrained")
    mesh = render_log_rainbow_field(
        ax,
        q_grid,
        w_grid,
        z_grid,
        shared_norm,
        title=rf"$|\Re \tilde{{\chi}}_e^L|$, $T={t_kelvin:,}\,\mathrm{{K}}$",
        fill_border=(t_kelvin == COLD_LADDER_T_K),
    )
    attach_log_colorbar(fig, ax, mesh, LADDER_CBAR_LABEL)
    saved = save_figure(fig, stem)
    plt.close(fig)
    return saved


def save_slice_panel() -> Path:
    fig, ax = plt.subplots(figsize=SINGLE_PANEL_FIGSIZE, layout="constrained")
    render_slice_panel(ax)
    saved = save_figure(fig, PANEL_C_STEM)
    plt.close(fig)
    return saved


def save_ratio_panel() -> Path:
    fig, ax = plt.subplots(figsize=SINGLE_PANEL_FIGSIZE, layout="constrained")
    mesh = render_bare_rpa_panel(ax)
    attach_log_colorbar(fig, ax, mesh, RATIO_CBAR_LABEL)
    saved = save_figure(fig, PANEL_D_STEM)
    plt.close(fig)
    return saved


def main() -> None:
    configure_seaborn()
    shared_norm = ladder_shared_norm(tuple(t for t, _ in LADDER_PANELS))

    saved_paths: list[Path] = []
    for t_kelvin, stem in LADDER_PANELS:
        saved_paths.append(save_ladder_panel(t_kelvin, stem, shared_norm))
    saved_paths.append(save_slice_panel())
    saved_paths.append(save_ratio_panel())

    for path in saved_paths:
        print(f"Saved {path}")


if __name__ == "__main__":
    main()
