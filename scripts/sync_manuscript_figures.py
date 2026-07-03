#!/usr/bin/env python3

"""Copy manuscript figure PDFs from output/ into manuscript/figures/."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_common import sync_manuscript_figures


def main() -> None:
    synced = sync_manuscript_figures()
    for path in synced:
        print(f"Synced {path}")


if __name__ == "__main__":
    main()
