# ==========================================================================
# MOSAIQ-LINDHARD ENGINE
# Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
# * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
# Free for non-commercial personal and academic research use only.
# Commercial, governmental, and public institutional use requires a
# separate paid license. See LICENSE file for details.
# * Contact: kim.ingee@bonasapiens.com
# ==========================================================================

"""Shared defaults for MosaiQ-Lindhard manuscript figure scripts."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_MISC_DIR = OUTPUT_DIR / "misc"
MANUSCRIPT_FIGURES_DIR = PROJECT_ROOT / "manuscript" / "figures"
DEFAULT_FIG_FORMAT = "pdf"

# PDFs referenced by manuscript/two-fermi.tex (APS submission bundle).
MANUSCRIPT_FIGURE_NAMES: tuple[str, ...] = (
    "chi_e_re",
    "chi_e_im",
    "chi_i_re",
    "chi_i_im",
    "chi_ei_re",
    "chi_ei_im",
    "S_ee_contour",
    "S_ii_contour",
    "S_ei_contour",
    "plasmon_dispersion",
    "zeta_rpa_dispersion",
    "t0_analytic_re_contour",
    "t0_analytic_im_contour",
    "t0_limit_validation",
    "thermal_anatomy_a",
    "thermal_anatomy_b",
    "thermal_anatomy_c",
    "thermal_anatomy_d",
    "Sq_gamma_ee",
    "Sq_gamma_ii",
    "Sq_gamma_ei",
    "lindhard_rs_t_excision",
)

PDF_RCPARAMS: dict[str, object] = {
    "savefig.format": DEFAULT_FIG_FORMAT,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def apply_pdf_rcparams(extra: dict[str, object] | None = None) -> None:
    rc = {**PDF_RCPARAMS}
    if extra:
        rc.update(extra)
    plt.rcParams.update(rc)


def output_path(name: str) -> Path:
    """Resolve a manuscript figure path under output/ using the default PDF format."""
    return OUTPUT_DIR / f"{Path(name).stem}.{DEFAULT_FIG_FORMAT}"


def format_temperature_k(t_kelvin: float) -> str:
    """Matplotlib mathtext temperature tag, e.g. $T = 10\\,000\\,\\mathrm{K}$."""
    if t_kelvin <= 0.0:
        return r"$T = 0$"
    t_int = int(round(t_kelvin))
    if t_int >= 1000:
        grouped = f"{t_int:,}".replace(",", r"\,")
        return rf"$T = {grouped}\,\mathrm{{K}}$"
    return rf"$T = {t_int}\,\mathrm{{K}}$"


def misc_output_path(name: str) -> Path:
    """Resolve an exploratory / non-manuscript figure path under output/misc/."""
    return OUTPUT_MISC_DIR / f"{Path(name).stem}.{DEFAULT_FIG_FORMAT}"


def publish_manuscript_figure(source: Path) -> Path:
    """Copy a generated figure into manuscript/figures/ for LaTeX and APS submission."""
    MANUSCRIPT_FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    destination = MANUSCRIPT_FIGURES_DIR / source.name
    if source.resolve() != destination.resolve():
        destination.write_bytes(source.read_bytes())
    return destination


def sync_manuscript_figures() -> list[Path]:
    """Copy all manuscript figure PDFs from output/ into manuscript/figures/."""
    synced: list[Path] = []
    missing: list[str] = []
    for stem in MANUSCRIPT_FIGURE_NAMES:
        source = OUTPUT_DIR / f"{stem}.{DEFAULT_FIG_FORMAT}"
        if not source.is_file():
            missing.append(source.name)
            continue
        synced.append(publish_manuscript_figure(source))
    if missing:
        missing_list = ", ".join(missing)
        raise FileNotFoundError(
            f"Missing manuscript figure(s) in {OUTPUT_DIR}: {missing_list}"
        )
    return synced


def save_figure(fig: Figure, name: str, **kwargs: object) -> Path:
    path = output_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    save_kwargs: dict[str, object] = {
        "format": DEFAULT_FIG_FORMAT,
        "dpi": PDF_RCPARAMS["savefig.dpi"],
        "bbox_inches": "tight",
    }
    save_kwargs.update(kwargs)
    fig.savefig(path, **save_kwargs)
    publish_manuscript_figure(path)
    return path
