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
DEFAULT_FIG_FORMAT = "pdf"

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


def misc_output_path(name: str) -> Path:
    """Resolve an exploratory / non-manuscript figure path under output/misc/."""
    return OUTPUT_MISC_DIR / f"{Path(name).stem}.{DEFAULT_FIG_FORMAT}"


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
    return path
